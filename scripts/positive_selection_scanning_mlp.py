"""
script to run MLP for positive selection scanning classification
references: dhruv raghavan, hannah snell
"""

import pandas as pd
import re
import numpy as np
import glob
import os
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler, LabelEncoder
from tensorflow.keras.models import Sequential
from tensorflow.keras.regularizers import l2
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization
from tensorflow.keras.optimizers import Adam
from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix, classification_report

def parse_filename(filename):
    """Extract population, generations, selection coefficient, and simulation number from filename"""
    basename = os.path.basename(filename)
    
    # Remove extension first - handle both .swifr and .swifr_final
    if basename.endswith('.swifr_final'):
        name_without_ext = basename.rsplit('.swifr_final', 1)[0]
    elif basename.endswith('.swifr'):
        name_without_ext = basename.rsplit('.swifr', 1)[0]
    else:
        return None, None, None, None
    
    # Split by dots
    parts = name_without_ext.split('.')
    
    if len(parts) >= 4:
        population = parts[0]
        generations = int(parts[1])
        
        # Handle selection coefficient - it might be split across two parts due to decimal point
        if len(parts) == 5:  # Format: pop.gens.sel1.sel2.simnum (e.g., YRI.300.0.05.22)
            selcoeff = float(f"{parts[2]}.{parts[3]}")
            sim_num = int(parts[4])
        elif len(parts) == 4:  # Format: pop.gens.sel.simnum (e.g., YRI.300.1.22)
            selcoeff = float(parts[2])
            sim_num = int(parts[3])
        else:
            return None, None, None, None
            
        return population, generations, selcoeff, sim_num
    
    return None, None, None, None

def load_data(base_path="data/deep_learning_sims_20251003", data_split="training_data"):

    data_path = os.path.join(base_path, data_split)
    
    # Site and window data paths
    site_path = os.path.join(data_path, "site_sims")
    window_path = os.path.join(data_path, "window_sims")
    
    # Site columns to keep (including XP-EHH column)
    site_features = ["SNP_name", "Physical_Distance", "Map_Distance", "XP-EHH", "DDAF", "nSL", "iHS", "FST"]
    
    # Window columns to keep
    window_features = ["SNP_name", "Theta_Pi", "Theta_W", "Tajima_D", "Fay_Wu_H", "Zeng_E", "Garud_H", "PBS", "NSS"]
    
    combined_data = []
    
    print(f"Processing data from {data_split}...")
    
    if not os.path.exists(site_path) or not os.path.exists(window_path):
        print(f"Warning: Paths not found - {site_path} or {window_path}")
        return pd.DataFrame()
        
    # Get all site and window files - handle both extensions
    site_files = sorted(glob.glob(os.path.join(site_path, "*.swifr")) + 
                       glob.glob(os.path.join(site_path, "*.swifr_final")))
    window_files = sorted(glob.glob(os.path.join(window_path, "*.swifr")) + 
                         glob.glob(os.path.join(window_path, "*.swifr_final")))
    
    print(f"Found {len(site_files)} site files and {len(window_files)} window files")
    
    # Process each site file and find its matching window file
    for site_file in tqdm(site_files, desc="Processing files"):
        # Get corresponding window file (should have same name)
        site_basename = os.path.basename(site_file)
        
        # Try to find matching window file with same extension
        window_file = os.path.join(window_path, site_basename)
        if not os.path.exists(window_file):
            # Try alternative extension
            if site_basename.endswith('.swifr_final'):
                alt_basename = site_basename.replace('.swifr_final', '.swifr')
            else:
                alt_basename = site_basename.replace('.swifr', '.swifr_final')
            window_file = os.path.join(window_path, alt_basename)
        
        if not os.path.exists(window_file):
            print(f"Warning: No matching window file for {site_file}")
            continue
            
        # Parse filename to get metadata
        pop, gens, selcoeff, sim_num = parse_filename(site_file)
        
        if pop is None:
            print(f"Warning: Could not parse filename {site_file}")
            continue
            
        try:
            # Load data - use raw string for regex pattern
            site_df = pd.read_csv(site_file, sep=r'\s+', header=0)
            window_df = pd.read_csv(window_file, sep=r'\s+', header=0)

            # Replace missing values
            site_df.replace(-998, np.nan, inplace=True)
            window_df.replace(-998, np.nan, inplace=True)

            # Ensure we have required columns
            if not all(col in site_df.columns for col in site_features):
                print(f"Warning: Missing site columns in {site_file}")
                print(f"Available columns: {site_df.columns.tolist()}")
                continue
            
            if not all(col in window_df.columns for col in window_features):
                print(f"Warning: Missing window columns in {window_file}")
                print(f"Available columns: {window_df.columns.tolist()}")
                continue
            
            # Process and combine data with proper labeling
            sim_combined = process_simulation_pair_with_labeling(
                site_df, window_df, site_features, window_features, 
                pop, gens, selcoeff, sim_num
            )
            
            if not sim_combined.empty:
                combined_data.append(sim_combined)
            
        except Exception as e:
            print(f"Error processing file {site_file}: {str(e)}")
    
    # Combine all data
    if combined_data:
        final_df = pd.concat(combined_data, ignore_index=True)
        print(f"\nFinal combined dataset shape: {final_df.shape}")
        print(f"Class distribution:\n{final_df['label'].value_counts()}")
        
        return final_df
    else:
        print("No data was successfully loaded!")
        return pd.DataFrame()

def process_simulation_pair_with_labeling(site_df, window_df, site_features, window_features, 
                                        population, generations, selcoeff, sim_num):
    """
    Process a single pair of site and window files with proper labeling based on position rules
    """
    # Filter to only keep needed columns
    site_df = site_df[site_features].copy()
    window_df = window_df[window_features].copy()
    
    # Sort both dataframes by SNP_name to ensure proper matching
    site_df = site_df.sort_values('SNP_name')
    window_df = window_df.sort_values('SNP_name')
    
    # Convert SNP_name to integers for proper comparison
    site_df['SNP_name'] = site_df['SNP_name'].astype(int)
    window_df['SNP_name'] = window_df['SNP_name'].astype(int)
    
    # Create window labels based on position rules
    window_df = window_df.copy()
    window_df['window_label'] = 'neutral'  # Default
    
    # Find the sweep window (position 480001)
    sweep_window_idx = window_df[window_df['SNP_name'] == 480001].index
    
    if len(sweep_window_idx) > 0:
        sweep_idx = sweep_window_idx[0]
        
        # Label sweep window
        window_df.loc[sweep_idx, 'window_label'] = 'sweep'
        
        # Label linked windows (before and after sweep window)
        if sweep_idx > 0:  # Window before
            window_df.loc[sweep_idx - 1, 'window_label'] = 'linked'
        if sweep_idx < len(window_df) - 1:  # Window after
            window_df.loc[sweep_idx + 1, 'window_label'] = 'linked'
    
    # Create a new dataframe to store the merged data
    merged_data = []
    
    # Iterate through window rows and find matching site rows
    for i, window_row in window_df.iterrows():
        window_start = window_row['SNP_name']
        window_end = window_start + 40000 - 1  # 40,000 sites per window
        window_label = window_row['window_label']
        
        # Find site rows within this window range
        matching_sites = site_df[(site_df['SNP_name'] >= window_start) & 
                                 (site_df['SNP_name'] <= window_end)].copy()
     
        if not matching_sites.empty:
            # Add window statistics to each matching site row
            for window_col in window_features:
                if window_col != 'SNP_name':  # Skip SNP_name to avoid duplication
                    matching_sites.loc[:, window_col] = window_row[window_col]
            
            # Label sites based on window type and position
            if window_label == 'sweep':
                # In sweep window: all sites are linked except center site (500001)
                matching_sites.loc[:, 'label'] = 'linked'
                center_site_mask = matching_sites['SNP_name'] == 500001
                matching_sites.loc[center_site_mask, 'label'] = 'sweep'
            else:
                # For linked and neutral windows: all sites get the window label
                matching_sites.loc[:, 'label'] = window_label
            
            # Add metadata
            matching_sites.loc[:, 'population'] = population
            matching_sites.loc[:, 'generations'] = generations
            matching_sites.loc[:, 'selcoeff'] = selcoeff
            matching_sites.loc[:, 'simulation'] = sim_num
            
            merged_data.append(matching_sites)
    
    if merged_data:
        return pd.concat(merged_data, ignore_index=True)
    else:
        # Return empty dataframe with correct columns if no matches found
        columns = (site_features + [col for col in window_features if col != 'SNP_name'] + 
                  ['label', 'population', 'generations', 'selcoeff', 'simulation'])
        return pd.DataFrame(columns=columns)

def balance_training_data_downsample(data):
    """
    Balance training data using random downsampling to match the smallest class
    """
    # Get class counts
    class_counts = data['label'].value_counts()
    print(f"Original class distribution:")
    print(class_counts)
    
    # Find the size of the smallest class
    min_class_size = class_counts.min()
    print(f"Smallest class size: {min_class_size}")
    print(f"Will downsample all classes to {min_class_size} samples")
    
    # Extract features for balancing
    feature_columns = ['FST', 'DDAF', 'iHS', 'nSL', 'XP-EHH',
                      'Theta_Pi', 'Theta_W', 'Tajima_D', 'Garud_H', 
                      'Fay_Wu_H', 'Zeng_E', 'PBS', 'NSS']
    
    # Clean data - remove rows with missing values
    data_clean = data.dropna(subset=feature_columns)
    print(f"Data shape after removing NaN: {data_clean.shape}")
    
    # Downsample each class to the minimum class size
    balanced_dfs = []
    
    for label in class_counts.index:
        class_data = data_clean[data_clean['label'] == label]
        
        if len(class_data) > min_class_size:
            # Randomly sample min_class_size samples
            downsampled = class_data.sample(n=min_class_size, random_state=42)
            print(f"Downsampled {label}: {len(class_data)} -> {len(downsampled)}")
        else:
            # Keep all samples if class size is already <= min_class_size
            downsampled = class_data
            print(f"Kept all {label} samples: {len(downsampled)}")
        
        balanced_dfs.append(downsampled)
    
    # Combine all downsampled classes
    balanced_data = pd.concat(balanced_dfs, ignore_index=True)
    
    # Shuffle the final dataset
    balanced_data = balanced_data.sample(frac=1, random_state=42).reset_index(drop=True)
    
    print(f"\nBalanced class distribution:")
    print(balanced_data['label'].value_counts())
    print(f"Total samples: {len(balanced_data)}")
    
    return balanced_data

def train_model(train_data, model_name="dl_model"):
    """Train the deep learning model"""
    # Extract features and labels
    feature_columns = ['FST', 'DDAF', 'iHS', 'nSL', 'XP-EHH',
                      'Theta_Pi', 'Theta_W', 'Tajima_D', 'Garud_H', 
                      'Fay_Wu_H', 'Zeng_E', 'PBS', 'NSS']
    
    # Handle missing values
    train_data_clean = train_data.dropna(subset=feature_columns)
    print(f"Training data shape after removing NaN: {train_data_clean.shape}")
    
    X_train = train_data_clean[feature_columns].values
    
    # Encode labels
    label_encoder = LabelEncoder()
    y_train = label_encoder.fit_transform(train_data_clean['label'])
    
    # Print the encoding mapping for reference
    print(f"Label encoding:")
    for i, label in enumerate(label_encoder.classes_):
        print(f"  {label}: {i}")
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    
    # Build model with regularization
    model = Sequential([
        Dense(64, activation='relu', kernel_regularizer=l2(0.01)),
        BatchNormalization(),
        Dropout(0.5),
        Dense(32, activation='relu', kernel_regularizer=l2(0.01)),
        BatchNormalization(),
        Dropout(0.5),
        Dense(16, activation='relu', kernel_regularizer=l2(0.01)),
        Dropout(0.3),
        Dense(len(label_encoder.classes_), activation='softmax')
    ])

    model.build(input_shape=(None, X_train_scaled.shape[1]))
    
    # Compile model
    model.compile(optimizer='adam',
                  loss='sparse_categorical_crossentropy',
                  metrics=['accuracy'])
    
    print(f"Model architecture:")
    model.summary()
    print(f"Features being used: {feature_columns}")
    
    # Train model
    history = model.fit(
        X_train_scaled, y_train,
        validation_split=0.2,
        epochs=25,
        batch_size=32,
        verbose=1
    )
    
    return model, scaler, label_encoder, history

def evaluate_model(model, scaler, label_encoder, test_data, test_case_name):
    """Evaluate model on test data and generate outputs"""
    # Extract features
    feature_columns = ['FST', 'DDAF', 'iHS', 'nSL', 'XP-EHH',
                      'Theta_Pi', 'Theta_W', 'Tajima_D', 'Garud_H', 
                      'Fay_Wu_H', 'Zeng_E', 'PBS', 'NSS']
    
    # Handle missing values
    test_data_clean = test_data.dropna(subset=feature_columns)
    print(f"Test data shape after removing NaN: {test_data_clean.shape}")
    
    if test_data_clean.empty:
        print(f"No valid test data for {test_case_name}")
        return None, None, None
    
    X_test = test_data_clean[feature_columns].values
    y_true = label_encoder.transform(test_data_clean['label'])
    
    # Scale features
    X_test_scaled = scaler.transform(X_test)
    
    # Make predictions
    y_pred_probs = model.predict(X_test_scaled)
    y_pred = np.argmax(y_pred_probs, axis=1)
    
    # Calculate metrics
    precision = precision_score(y_true, y_pred, average='weighted')
    recall = recall_score(y_true, y_pred, average='weighted')
    f1 = f1_score(y_true, y_pred, average='weighted')
    
    print(f'\nResults for {test_case_name}:')
    print(f'Precision: {precision:.4f}')
    print(f'Recall: {recall:.4f}')
    print(f'F1-Score: {f1:.4f}')
    
    # Confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    labels = label_encoder.classes_
    
    print(f'\nConfusion Matrix for {test_case_name}:')
    print(f'Labels: {labels}')
    print(cm)
    
    # Create detailed classification report
    report = classification_report(y_true, y_pred, target_names=labels, output_dict=True)
    report_df = pd.DataFrame(report).transpose()
    print(f'\nClassification Report for {test_case_name}:')
    print(report_df)
    
    # Create prediction results dataframe for saving
    results_df = test_data_clean.copy()
    results_df['predicted_label'] = label_encoder.inverse_transform(y_pred)
    
    # Add probabilities for each class
    for i, class_name in enumerate(label_encoder.classes_):
        results_df[f'predicted_prob_{class_name}'] = y_pred_probs[:, i]
    
    # Create results directory if it doesn't exist
    os.makedirs('results', exist_ok=True)
    
    # Save raw predictions for plotting
    output_file = f'results/predictions_{test_case_name}.csv'
    results_df.to_csv(output_file, index=False)
    print(f"Saved prediction results to {output_file}")
    
    return test_data_clean, cm, report_df

def main():
    """Main pipeline execution"""
    print("Deep Learning Pipeline for Genomic Selection Detection")
    print("="*70)

    # Load training data
    print("Loading training data...")
    train_data = load_data(data_split="training_data")
    
    if train_data.empty:
        print("No training data loaded. Exiting.")
        return
    
    # Load testing data
    print("\nLoading testing data...")
    test_data = load_data(data_split="testing_data")
    
    if test_data.empty:
        print("No testing data loaded. Exiting.")
        return
    
    # Balance training data using downsampling
    print("\nBalancing training data using downsampling...")
    try:
        balanced_train_data = balance_training_data_downsample(train_data)
    except Exception as e:
        print(f"Failed to balance with downsampling: {e}")
        balanced_train_data = train_data
    
    # Train model
    print("\nTraining model with downsampled balanced data...")
    model, scaler, label_encoder, history = train_model(balanced_train_data, model_name="dl_model_downsampled")
    
    # Test on specific cases
    test_cases = [
        (300, 0.05, "300gens_0.05selcoeff_downsampled"),
        (600, 0.02, "600gens_0.02selcoeff_downsampled")
    ]
    
    for gens, selcoeff, case_name in test_cases:
        print(f"\n{'-'*50}")
        print(f"Testing on {case_name}")
        print(f"{'-'*50}")
        
        # Filter test data for this case
        test_case_data = test_data[
            (test_data['generations'] == gens) & 
            (test_data['selcoeff'] == selcoeff)
        ]
        
        if test_case_data.empty:
            print(f"No test data found for {case_name}")
            continue
            
        print(f"Test data shape for {case_name}: {test_case_data.shape}")
        print(f"Class distribution:")
        print(test_case_data['label'].value_counts())
        
        # Evaluate model on this test case
        pred_results, cm, report = evaluate_model(
            model, scaler, label_encoder, test_case_data, case_name
        )

    print(f"\nTraining History Summary:")
    print(f"Final training loss: {history.history['loss'][-1]:.4f}")
    print(f"Final validation loss: {history.history['val_loss'][-1]:.4f}")
    print(f"Final training accuracy: {history.history['accuracy'][-1]:.4f}")
    print(f"Final validation accuracy: {history.history['val_accuracy'][-1]:.4f}")
    
    print("\nPipeline completed successfully!")
    return model, scaler, label_encoder, history

if __name__ == "__main__":
    model, scaler, label_encoder, history = main()
