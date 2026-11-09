"""
script to run MLP for positive selection scanning classification
references: dhruv raghavan, hannah snell

For pooled exports: python positive_selection_scanning_mlp.py DATA --output RESULTS --sampling both
"""

import pandas as pd
import argparse
from pathlib import Path
import numpy as np
import glob
import os
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from tensorflow.keras.models import Sequential
from tensorflow.keras.regularizers import l2
from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization
from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix, classification_report

CLASSES = ("linked", "neutral", "sweep")  # Alphabetical, as in LabelEncoder.
SITE_FEATURES = ["FST", "DDAF", "iHS", "nSL", "XP-EHH"]
WINDOW_FEATURES = [
    "Theta_Pi", "Theta_W", "Tajima_D", "Garud_H", "Fay_Wu_H", "Zeng_E", "PBS", "NSS",
]
FEATURES = SITE_FEATURES + WINDOW_FEATURES
WINDOW_ALIASES = {
    "Tajimas_D": "Tajima_D",
    "Garuds_H": "Garud_H",
    "Fay_Wus_H": "Fay_Wu_H",
    "Zengs_E": "Zeng_E",
}
WINDOW_SIZE = 40_000


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

    from tqdm import tqdm

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

def read_export(root, level, label, split):
    path = root / level / f"{level}_{label}.{split}"
    table = pd.read_csv(path, sep=r"\s+", na_values=[-998], keep_default_na=True)
    columns = SITE_FEATURES if level == "site" else WINDOW_FEATURES
    table = table.rename(columns={"SNP_name": "position", **WINDOW_ALIASES})
    table = table[["position", *columns]]
    if table.empty or table["position"].isna().any():
        raise ValueError(f"Missing coordinates in {path}")
    positions = table["position"].to_numpy()
    if label == "sweep":
        expected = 500_001 if level == "site" else 480_001
        if not np.all(positions == expected):
            raise ValueError(f"Expected one sweep coordinate per block in {path}")
        table["block"] = np.arange(1, len(table) + 1)
    else:
        # Equal adjacent positions can occur within the same simulation.
        table["block"] = table["position"].diff().lt(0).cumsum() + 1
    table["label"] = label
    return table


def load_export_data(root, split):
    site = {label: read_export(root, "site", label, split) for label in CLASSES}
    window = {label: read_export(root, "window", label, split) for label in CLASSES}
    n_blocks = len(site["sweep"])
    for table in [*site.values(), *window.values()]:
        if table["block"].iat[-1] != n_blocks:
            raise ValueError(f"Inconsistent simulation block counts in {split}")

    sites = pd.concat(site.values(), ignore_index=True)
    windows = pd.concat(window.values(), ignore_index=True).rename(
        columns={"position": "window_start", "label": "window_label"}
    )
    sites["window_start"] = (sites["position"] - 1) // WINDOW_SIZE * WINDOW_SIZE + 1
    paired = sites.merge(
        windows, on=["block", "window_start"], how="left", validate="many_to_one",
        indicator=True, sort=False,
    )
    if not paired["_merge"].eq("both").all():
        raise ValueError(f"Some {split} sites have no corresponding window")
    linked_in_sweep = paired["label"].eq("linked") & paired["window_label"].eq("sweep")
    if not (paired["label"].eq(paired["window_label"]) | linked_in_sweep).all():
        raise ValueError(f"Site and window classes disagree in {split}")

    # Check inferred block order against the number of segregating sites (NSS).
    keys = ["block", "window_start"]
    observed = sites.groupby(keys).size().rename("site_count")
    audit = windows.join(observed, on=keys)
    audit["site_count"] = audit["site_count"].fillna(0)
    concordance = audit.groupby("window_start")[["NSS", "site_count"]].apply(
        lambda group: group["NSS"].corr(group["site_count"])
    )
    if not concordance.ge(0.5).all():
        raise ValueError(f"Could not verify site/window block alignment in {split}")

    usable = np.isfinite(paired[FEATURES]).all(axis=1)
    clean = paired.loc[usable, ["block", "position", "window_start", "label", *FEATURES]]
    counts = pd.DataFrame({
        "raw": paired["label"].value_counts(),
        "complete": clean["label"].value_counts(),
    }).fillna(0).astype(int)
    print(f"{split}: {n_blocks} blocks, minimum NSS / site-count r={concordance.min():.3f}")
    print(counts.to_string(), flush=True)
    return clean, n_blocks


def load_scenarios(path, n_blocks):
    if path is None:
        return None
    mapping = pd.read_csv(path, dtype={"scenario": "string"})
    if (mapping["block"].duplicated().any()
            or set(mapping["block"]) != set(range(1, n_blocks + 1))
            or mapping["scenario"].fillna("").str.strip().eq("").any()):
        raise ValueError("Supply one nonempty scenario per 1-based test block")
    return mapping.set_index("block")["scenario"].str.strip()


def balance_training_data_downsample(data, *, seed=None):
    """
    Balance training data using random downsampling to match the smallest class
    """
    if seed is not None:
        groups = data.groupby('label', sort=True)
        min_class_size = groups.size().min()
        return groups.sample(n=min_class_size, random_state=seed).sample(frac=1, random_state=seed)

    # Get class counts
    class_counts = data['label'].value_counts()
    print(f"Original class distribution:")
    print(class_counts)
    
    # Find the size of the smallest class
    min_class_size = class_counts.min()
    print(f"Smallest class size: {min_class_size}")
    print(f"Will downsample all classes to {min_class_size} samples")
    
    # Extract features for balancing
    feature_columns = FEATURES
    
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

def train_model(train_data, model_name="dl_model", *, validation_data=None):
    """Train the deep learning model"""
    # Extract features and labels
    feature_columns = FEATURES
    
    # Handle missing values
    train_data_clean = train_data.dropna(subset=feature_columns)
    print(f"Training data shape after removing NaN: {train_data_clean.shape}")
    
    X_train = train_data_clean[feature_columns]
    if validation_data is None:
        X_train = X_train.values
    
    # Encode labels
    label_encoder = LabelEncoder()
    y_train = label_encoder.fit_transform(train_data_clean['label'])
    
    # Print the encoding mapping for reference
    print(f"Label encoding:")
    for i, label in enumerate(label_encoder.classes_):
        print(f"  {label}: {i}")
    
    # Scale features
    scaler = StandardScaler()
    if validation_data is None:
        X_train_scaled = scaler.fit_transform(X_train)
        validation_options = {'validation_split': 0.2}
        steps_per_execution = 1
        verbose = 1
    else:
        scaler.fit(X_train)
        X_train_scaled = scaler.transform(X_train).astype(np.float32)
        y_train = y_train.astype(np.int32)
        X_val_scaled = scaler.transform(validation_data[feature_columns]).astype(np.float32)
        y_val = label_encoder.transform(validation_data['label']).astype(np.int32)
        validation_options = {'validation_data': (X_val_scaled, y_val)}
        steps_per_execution = 128
        verbose = 2
    
    # Build model with regularization
    model = Sequential([
        Input(shape=(X_train_scaled.shape[1],)),
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

    # Compile model
    model.compile(optimizer='adam',
                  loss='sparse_categorical_crossentropy',
                  metrics=['accuracy'],
                  steps_per_execution=steps_per_execution)
    
    print(f"Model architecture:")
    model.summary()
    print(f"Features being used: {feature_columns}")
    
    # Train model
    history = model.fit(
        X_train_scaled, y_train,
        epochs=25,
        batch_size=32,
        verbose=verbose,
        **validation_options
    )
    
    return model, scaler, label_encoder, history

def evaluate_model(model, scaler, label_encoder, test_data, test_case_name):
    """Evaluate model on test data and generate outputs"""
    # Extract features
    feature_columns = FEATURES
    
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

def metrics(test, mode):
    cohorts = [("Pooled test", test)]
    if "scenario" in test:
        cohorts.extend(test.groupby("scenario", sort=True))
    for scenario, cohort in cohorts:
        matrix = confusion_matrix(cohort["label"], cohort["predicted_label"], labels=CLASSES)
        for label in ("neutral", "linked", "sweep"):
            idx = CLASSES.index(label)
            tp = matrix[idx, idx]
            fp = matrix[:, idx].sum() - tp
            fn = matrix[idx, :].sum() - tp
            yield {
                "Scenario": scenario, "Level": "Site",
                "Tool": "MLP with downsampling" if mode == "downsample" else "MLP without downsampling",
                "Class": label.title(), "TP": tp, "FP": fp, "FN": fn,
                "Precision": tp / (tp + fp) if tp + fp else 0.0,
                "Recall": tp / (tp + fn) if tp + fn else 0.0,
                "F1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0,
            }


def evaluate_export_model(model, scaler, label_encoder, test_data, mode, scenarios, output):
    """Evaluate the MLP on the pooled export and any supplied test scenarios."""
    X_test_scaled = scaler.transform(test_data[FEATURES]).astype(np.float32)
    probabilities = model.predict(X_test_scaled, batch_size=4096, verbose=0)
    if not np.isfinite(probabilities).all():
        raise ValueError(f"Nonfinite probabilities in {mode} run")

    predictions = test_data[['block', 'position', 'window_start', 'label']].copy()
    if scenarios is not None:
        predictions['scenario'] = predictions['block'].map(scenarios)
    predictions['predicted_label'] = label_encoder.inverse_transform(probabilities.argmax(axis=1))
    for idx, label in enumerate(label_encoder.classes_):
        predictions[f'prob_{label}'] = probabilities[:, idx]

    results = pd.DataFrame(metrics(predictions, mode))
    print(results.to_string(index=False), flush=True)
    results.to_csv(output / f'metrics_{mode}.csv', index=False)
    predictions.to_csv(
        output / f'predictions_{mode}.csv.gz', index=False,
        compression={'method': 'gzip', 'compresslevel': 1},
    )


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

def run_export_pipeline(args):
    """Train and evaluate on the pooled class-split exports."""
    tf.config.threading.set_inter_op_parallelism_threads(2)
    tf.config.threading.set_intra_op_parallelism_threads(2)
    train, n_train = load_export_data(args.data, "train")
    test, n_test = load_export_data(args.data, "test")
    scenarios = load_scenarios(args.test_scenarios, n_test)
    fit_blocks, validation_blocks = train_test_split(
        np.arange(1, n_train + 1), test_size=0.2, random_state=args.seed
    )
    is_fit = train["block"].isin(fit_blocks)
    fit, validation = train.loc[is_fit], train.loc[~is_fit]
    if set(fit["label"]) != set(CLASSES) or set(test["label"]) != set(CLASSES):
        raise ValueError("Expected all three classes in fitting and test data")
    print(f"Fit blocks: {len(fit_blocks)}. Validation blocks: {len(validation_blocks)}.", flush=True)
    print(f"Validation rows: {validation['label'].value_counts().to_dict()}", flush=True)

    args.output.mkdir(parents=True, exist_ok=True)
    modes = ["downsample", "none"] if args.sampling == "both" else [args.sampling]
    for mode in modes:
        tf.keras.backend.clear_session()
        tf.keras.utils.set_random_seed(args.seed)
        train_data = balance_training_data_downsample(fit, seed=args.seed) if mode == 'downsample' else fit
        print(f"{mode} fitting rows: {train_data['label'].value_counts().to_dict()}", flush=True)
        model, scaler, label_encoder, history = train_model(train_data, validation_data=validation)
        history_table = pd.DataFrame(history.history)
        history_table.index = pd.RangeIndex(1, len(history_table) + 1, name='epoch')
        history_table.to_csv(args.output / f'history_{mode}.csv')
        evaluate_export_model(model, scaler, label_encoder, test, mode, scenarios, args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MLP for positive selection scanning classification")
    parser.add_argument('data', nargs='?', type=Path, help='Directory containing the site/ and window/ exports')
    parser.add_argument('--output', type=Path, help='Directory for tables and predictions')
    parser.add_argument('--sampling', choices=['both', 'downsample', 'none'])
    parser.add_argument('--seed', type=int)
    parser.add_argument('--test-scenarios', type=Path, help='Optional CSV with columns block,scenario')
    args = parser.parse_args()

    if args.data is None:
        if any(value is not None for value in (args.output, args.sampling, args.seed, args.test_scenarios)):
            parser.error("Supply the pooled data directory when using export options")
        main()
    else:
        if args.output is None:
            parser.error("--output is required with pooled data")
        args.sampling = args.sampling or 'both'
        args.seed = 42 if args.seed is None else args.seed
        run_export_pipeline(args)
