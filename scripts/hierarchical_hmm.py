"""
script that encodes WINDEX's hierarchical hmm framework
references: dr. lauren sugden, scott mccallum
"""

import numpy as np
from classifiers import Classifiers, Stats
import math
import sys
import pandas as pd
import argparse
import copy
import glob
import re
import os
import copy 
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

_shared_model = None
_worker_model = None 

def _init_worker(): 

	global _worker_model
	_worker_model = _shared_model

def _worker_task(seq_idx, window_trans_mat): 

	global _worker_model
	m = _worker_model 
	seq_data = m.training_sequences[seq_idx]

	m.Svec_sites = seq_data['Svec_sites']
	m.Svec_windows = seq_data['Svec_windows']
	m.positions = seq_data['positions']
	m.windows = seq_data['windows']
	m.window_trans_mat = window_trans_mat # updates each m-step 

	xi, gamma, lp = m.xi_gamma_window()

	# reduce to per-sequence sums so sequences of different lengths combine correctly
	W = len(m.window_state_list)
	T = len(gamma)
	xi_sum_seq = [['n/a'] * W for _ in range(W)]
	gamma_sum_seq = ['n/a'] * W

	for w in range(1, W): 
		g_vals = [gamma[t][w] for t in range(T - 1) if gamma[t][w] != 'n/a']
		if g_vals: 
			gamma_sum_seq[w] = m.logsumexp(g_vals)
		for w2 in range(1, W): 
			x_vals = [xi[t][w][w2] for t in range(len(xi)) if xi[t][w][w2] != 'n/a']
			if x_vals: 
				xi_sum_seq[w][w2] = m.logsumexp(x_vals)

	return xi_sum_seq, gamma_sum_seq, lp

class WINDEX():
    
	def __init__(self, path2trained_sites, path2trained_windows, datafile_sites, datafile_windows, window_size, baum_welch_training, baum_welch_matrix, alpha = 1.0, emissiontype = "ODE"):
		self.name = None
		self.alpha = alpha
		self.path2trained_sites = path2trained_sites
		self.path2trained_windows = path2trained_windows
		self.datafile_sites = datafile_sites
		self.datafile_windows = datafile_windows
		self.window_size = window_size
		self.emissiontype = emissiontype

		if baum_welch_training == False: 
			assert datafile_sites is not None and datafile_windows is not None, \
            	"datafile_sites/datafile_windows are required unless --baum_welch_training is set"
			self.positions, self.Svec_sites, self.classifierobj_sites = self.read_site_statfiles()
			self.windows, self.Svec_windows, self.classifierobj_windows = self.read_window_statfiles()
			self.window_trans_mat, self.window_state_list, self.site_trans_mat_neutral, self.site_states_neutral_list, self.site_trans_mat_linkedleft, self.site_states_linkedleft_list, self.site_trans_mat_sweep, self. site_states_sweep_list, self.site_trans_mat_linkedright, self.site_states_linkedright_list = self.set_transitions(baum_welch_training, window_trans_mat_path = baum_welch_matrix)
		else:
			self.window_trans_mat, self.window_state_list, self.site_trans_mat_neutral, self.site_states_neutral_list, self.site_trans_mat_linkedleft, self.site_states_linkedleft_list, self.site_trans_mat_sweep, self. site_states_sweep_list, self.site_trans_mat_linkedright, self.site_states_linkedright_list = self.set_transitions(baum_welch_training, window_trans_mat_path = baum_welch_matrix)

		self.training_sequences = [] # will get populated by load_training_data if user specifies

	def read_site_statfiles(self, classifierobj_sites = None): 
		if classifierobj_sites is None:
			classifierobj_sites = Classifiers(self.path2trained_sites)

		with open(self.path2trained_sites + 'component_stats.txt') as f_sites: 
			statlist_sites = f_sites.read().strip().splitlines() 

			usecols_sites = ['SNP_name'] + statlist_sites 
			df_sites = pd.read_csv(self.datafile_sites, sep = '\t', usecols = usecols_sites)
			positions_all_sites = df_sites['SNP_name'].astype(int).to_numpy() 
			stat_matrix_sites = df_sites[statlist_sites].to_numpy(dtype = np.float64)

			keep_mask_sites = (stat_matrix_sites != -998).any(axis = 1) # making mask to remove sites with no defined stats
			n_removed_sites = int((~keep_mask_sites).sum())

			if n_removed_sites: 
				removed_sites = positions_all_sites[~keep_mask_sites]
				print(f"removing {n_removed_sites} sites with no defined statistics")

			positions_sites = positions_all_sites[keep_mask_sites].tolist()
			kept_stats_sites = stat_matrix_sites[keep_mask_sites]

			Svec_sites = []
			for row_sites in kept_stats_sites:
				S = Stats(self.path2trained_sites)
				S.stat2score = dict(zip(statlist_sites, row_sites.tolist()))
				Svec_sites.append(S)

			return positions_sites, Svec_sites, classifierobj_sites
	
	def read_window_statfiles(self, classifierobj_windows = None): 
		if classifierobj_windows is None:
			classifierobj_windows = Classifiers(self.path2trained_windows)

		with open(self.path2trained_windows + 'component_stats.txt') as f_windows: 
			statlist_windows = f_windows.read().strip().splitlines() 

			usecols_windows = ['SNP_name'] + statlist_windows 
			df_windows = pd.read_csv(self.datafile_windows, sep = '\t', usecols = usecols_windows)
			positions_all_windows = df_windows['SNP_name'].astype(int).to_numpy() 
			stat_matrix_windows = df_windows[statlist_windows].to_numpy(dtype = np.float64)

			keep_mask_windows = (stat_matrix_windows != -998).any(axis = 1) # making mask to remove windows with no defined stats
			n_removed_windows = int((~keep_mask_windows).sum())

			if n_removed_windows: 
				removed_windows = positions_all_windows[~keep_mask_windows]
				print(f"removing {n_removed_windows} windows with no defined statistics")

			positions_windows = positions_all_windows[keep_mask_windows].tolist()
			kept_stats_windows = stat_matrix_windows[keep_mask_windows]

			Svec_windows = []
			for row_windows in kept_stats_windows:
				S = Stats(self.path2trained_windows)
				S.stat2score = dict(zip(statlist_windows, row_windows.tolist()))
				Svec_windows.append(S)

			return positions_windows, Svec_windows, classifierobj_windows

	def load_training_data(self, training_datafiles_sites, training_datafiles_windows): 
		
		# check if there are the same number of site and window files 
		if len(training_datafiles_sites) != len(training_datafiles_windows): 
			raise ValueError("must have same number of site and window training files!")

		# original_test_sites = self.datafile_sites
		# original_test_windows = self.datafile_windows

		shared_classifier_sites = Classifiers(self.path2trained_sites)
		shared_classifier_windows = Classifiers(self.path2trained_windows)

		self.training_sequences = []

		for df_sites, df_windows in zip(training_datafiles_sites, training_datafiles_windows): 
			self.datafile_sites = df_sites
			self.datafile_windows = df_windows

			positions, Svec_sites, _ = self.read_site_statfiles(shared_classifier_sites)
			windows, Svec_windows, _ = self.read_window_statfiles(shared_classifier_windows)

			self.training_sequences.append({
            'Svec_sites' : Svec_sites,
            'Svec_windows' : Svec_windows,
            'positions' : positions,
            'windows' : windows
        	})

		self.classifierobj_sites = shared_classifier_sites
		self.classifierobj_windows = shared_classifier_windows

		# self.datafile_sites = original_test_sites
		# self.datafile_windows = original_test_windows
		# self.positions, self.Svec_sites, _ = self.read_site_statfiles(shared_classifier_sites)
		# self.windows, self.Svec_windows, _ = self.read_window_statfiles(shared_classifier_windows)

		print(f"Loaded {len(self.training_sequences)} training sequences.\n")

	def load_window_trans_mat(self, path):
		rows = []
		with open(path, 'r') as f_in:
			for line in f_in:
				row = [
					np.nan if val == 'n/a' else float(val)
					for val in line.rstrip('\n').split('\t')
				]
				rows.append(row)
		return np.array(rows)

	def set_transitions(self, baum_welch_training = True, window_trans_mat_path = None):

		if baum_welch_training == True: 
			# set window trans mat initialization
			#window_trans_mat = np.array([[0, 0.9, 0.1, 0, 0], [0.001, 0.998, 0.001, 0, 0], [0, 0, 0.9, 0.1, 0], [0, 0, 0, 0, 1], [0.001, 0.001, 0, 0, 0.998]]) 
			#window_trans_mat = np.array([[0, 0.999, 0.001, 0, 0], [0.001, 0.998, 0.001, 0, 0], [0, 0, 0.998, 0.002, 0], [0, 0, 0, 0, 1], [0.001, 0.001, 0, 0, 0.998]]) 
			#window_trans_mat = np.array([[0, 0.9, 0.1, 0, 0], [0.001, 0.899, 0.1, 0, 0], [0, 0, 0.5, 0.5, 0], [0, 0, 0, 0, 1], [0.001, 0.499, 0, 0, 0.5]]) 
			#window_trans_mat = np.array([[0, 0.98, 0.02, 0, 0], [0.019, 0.962, 0.019, 0, 0], [0, 0, 0.8, 0.2, 0], [0, 0, 0, 0, 1], [0.01, 0.19, 0, 0, 0.8]])
			window_trans_mat = np.array([[0, 0.95, 0.05, 0, 0], [0.02, 0.855, 0.125, 0, 0], [0, 0, 0.875, 0.125, 0], [0, 0, 0, 0, 1], [0.01, 0.125, 0, 0, 0.865]])

		else:
			# window matrix -> read in existing baum-welch trained matrix or use general default 
			if window_trans_mat_path is not None:
				window_trans_mat = self.load_window_trans_mat(window_trans_mat_path)
			else:
				#window_trans_mat = np.array([[0, 0.999, 0.001, 0, 0], [0.0001, 0.9989, 0.001, 0, 0], [0, 0, 0.999, 0.001, 0], [0, 0, 0, 0, 1], [0.001, 0.001, 0, 0, 0.998]]) 
				#window_trans_mat = np.array([[0, 0.0001, 0.9999, 0, 0], [0.0001, 0, 0.9999, 0, 0], [0, 0, 0.96, 0.04, 0], [0, 0, 0, 0, 1], [0.0001, 0.0001, 0, 0, 0.9998]]) 
				#window_trans_mat = np.array([[0, 0.0000001, 0.9999999, 0, 0], [0.0000001, 0.0000001, 0.9999998, 0, 0], [0, 0, 0.96, 0.04, 0], [0, 0, 0, 0, 1], [0.0000001, 0.0000001, 0, 0, 0.9999998]]) 
				#window_trans_mat = np.array([[0, 0.0000000000001, 0.9999999999999, 0, 0], [0.0000000000001, 0.0000000000001, 0.9999999999998, 0, 0], [0, 0, 0.96, 0.04, 0], [0, 0, 0, 0, 1], [0.0000000000001, 0.0000000000001, 0, 0, 0.9999999999998]]) 
				#window_trans_mat = np.array([[0, 0.9, 0.1, 0, 0], [0.001, 0.899, 0.1, 0, 0], [0, 0, 0.5, 0.5, 0], [0, 0, 0, 0, 1], [0.001, 0.499, 0, 0, 0.5]]) 
				#window_trans_mat = np.array([[0, 0.98, 0.02, 0, 0], [0.019, 0.962, 0.019, 0, 0], [0, 0, 0.8, 0.2, 0], [0, 0, 0, 0, 1], [0.01, 0.19, 0, 0, 0.8]])
				#window_trans_mat = np.array([[0, 0.98, 0.02, 0, 0], [0.019, 0.962, 0.019, 0, 0], [0.005, 0.095, 0.8, 0.1, 0], [0, 0, 0, 0, 1], [0.01, 0.19, 0, 0, 0.8]])
				window_trans_mat = np.array([[0, 0.95, 0.05, 0, 0], [0.02, 0.855, 0.125, 0, 0], [0, 0, 0.875, 0.125, 0], [0, 0, 0, 0, 1], [0.01, 0.125, 0, 0, 0.865]])


		window_states = {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
		window_state_list = list(window_states.values())

		# site transitions for window state = neutral
		site_trans_mat_neutral = [[0, 1], [0.001, 0.999]]
		site_states_neutral = {'Start/End':'none', 'Neutral':'neutral'}
		site_states_neutral_list = list(site_states_neutral.values())

		# site transitions for window state = linked left
		site_trans_mat_linkedleft = [[0, 1], [0.001, 0.999]]
		site_states_linkedleft= {'Start/End':'none', 'LinkedLeft': 'linked'}
		site_states_linkedleft_list = list(site_states_linkedleft.values())

		# sweep site transitions are parameterized by specified window size
		### three states: LL, S, LR
		### expectation: 1 sweep state transition per sweep window, so sweep prob should be 1/average num segregating sites per window (using 100 here)

		# site transitions for window state = sweep
		#site_trans_mat_sweep = [[0, 0.999, 0.001, 0], [0, 0.99, 0.01, 0], [0, 0, 0, 1], [0.001, 0, 0, 0.999]]
		site_trans_mat_sweep = [[0, 0.995, 0.005, 0], [0, 0.99, 0.01, 0], [0, 0, 0, 1], [0.01, 0, 0, 0.99]]
		site_states_sweep = {'Start/End':'none', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
		site_states_sweep_list = list(site_states_sweep.values())

		# site matrix for window state = linked right
		site_trans_mat_linkedright = [[0, 1], [0.001, 0.999]]
		site_states_linkedright= {'Start/End':'none', 'LinkedRight':'linked'}
		site_states_linkedright_list = list(site_states_linkedright.values())

		self.window_mask = (window_trans_mat > 0)

		return window_trans_mat, window_state_list, site_trans_mat_neutral, site_states_neutral_list, site_trans_mat_linkedleft, site_states_linkedleft_list, site_trans_mat_sweep, site_states_sweep_list, site_trans_mat_linkedright, site_states_linkedright_list

	def sweep_site_matrix(self, n_sites):
		
		# adjusting the sweep site probability to be adaptive to the number of sites in a given window
		p = 1.0 / max(int(n_sites), 2)
		m = [row[:] for row in self.site_trans_mat_sweep]
		m[0][1] = 1 - p; m[0][2] = p 
		m[1][1] = 1 - p; m[1][2] = p
		
		return m
	
	def get_nonempty_position_indices(self, win_start, win_end):
		
		positions = np.array(self.positions)
		return list(np.where((positions >= win_start) & (positions < win_end))[0])

	def window_baum_welch(self, n_iter = 20, tol = 1e-6, n_workers = None): 
		""" 
		baum-welch training for window level transitions, just for WINDEX transition probabilities 
		code implementation assisted with Claude Sonnet v4.6
		"""

		global _shared_model
		n_workers = n_workers or int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count())) # make this parallelizable with cpus 
		n_seqs = len(self.training_sequences)
		W = len(self.window_state_list)

		# setting alpha parameter to better balance between site and window
		total_pos = sum(len(s['positions']) for s in self.training_sequences)
		total_win = sum(len(s['windows']) for s in self.training_sequences)
		self.alpha = 1.0 / (total_pos / total_win)
		print(f"mean sites/window = {total_pos / total_win:.1f}, alpha = {self.alpha:.6f}", flush = True)

		_shared_model = self

		print(f"starting baum-welch training...\n", flush = True)

		ctx = mp.get_context("fork")
		with ProcessPoolExecutor(
			max_workers = n_workers, 
			mp_context = ctx,
			initializer = _init_worker,
		) as pool: 

			for iteration in range(n_iter): 
				prev_window = copy.deepcopy(self.window_trans_mat)

				futures = [
					pool.submit(_worker_task, i, self.window_trans_mat)
					for i in range(n_seqs)
				]
				results = [f.result() for f in futures] # this is where the _worker_task calls the gamma_xi_window function, which also calls the hierarchical forward/backward algs

				# m-step
				xi_sum = [['n/a'] * W for _ in range(W)]
				gamma_sum = ['n/a'] * W
				total_log_prob = 0.0

				for xi_s, gamma_s, lp in results: 
					total_log_prob += lp 
					for w in range(1, W): 
						if gamma_s[w] != 'n/a': 
							gamma_sum[w] = (self.logsumexp([gamma_sum[w], gamma_s[w]]) if gamma_sum[w] != 'n/a' else gamma_s[w])
						for w2 in range(1, W): 
							if xi_s[w][w2] != 'n/a': 
								xi_sum[w][w2] = (self.logsumexp([xi_sum[w][w2], xi_s[w][w2]]) if xi_sum[w][w2] != 'n/a' else xi_s[w][w2])

				for w in range(1, W): 
					if gamma_sum[w] == 'n/a': 
						continue

					# floor transitions that are allowed
					counts = np.zeros(W)
					for w2 in range(1, W): 
						if xi_sum[w][w2] != 'n/a': 
							counts[w2] = math.exp(xi_sum[w][w2] - gamma_sum[w])
					counts[~self.window_mask[w]] = 0.0
					counts[self.window_mask[w]] += 1e-6

					# rescale to the mass not taken by the end transition in column 0
					target = 1.0 - self.window_trans_mat[w][0]
					self.window_trans_mat[w][1:W] = counts[1:W] * target / counts[1:W].sum()
			
				# convergence check:
				delta_win = max(abs(self.window_trans_mat[w][w2] - prev_window[w][w2]) for w in range(1, W) for w2 in range(W))

				print(f"delta_window = {delta_win:.8f}")
				print(f"Iter {iteration + 1:>3d} | log P(O|λ) = {total_log_prob:.4f} | "f"max |Δ window| = {delta_win:.6f}", flush = True)

				# write window transition matrix every iteration
				output_path = "baum_welch_trained_window_transitions_1kgp_init_revisit2.txt"
				with open(output_path, 'w') as f_out:
					for w in range(W):
						row_vals = [
							repr(float(self.window_trans_mat[w][w2])) if self.window_trans_mat[w][w2] != 'n/a' else 'n/a' for w2 in range(W)]
						f_out.write("\t".join(row_vals) + "\n")

				if delta_win < tol:
					print(f"Converged after {iteration + 1} iterations.", flush = True)
					print(f"Wrote final window trans mat to {output_path}")
					break

		return self.window_trans_mat, self.site_trans_mat_neutral, self.site_trans_mat_linkedleft, self.site_trans_mat_sweep, self.site_trans_mat_linkedright

	def hierarchical_baum_welch(self, n_iter = 20, tol = 1e-6): 
		""" 
		hierarchical extension of baum-welch training, just for WINDEX transition probabilities 
		code implementation assisted with Claude Sonnet v4.6
		"""

		print(f"starting baum-welch training...\n", flush = True)

		W = len(self.window_state_list)

		for iteration in range(n_iter): 
			prev_window = copy.deepcopy(self.window_trans_mat)
			prev_site = {w: copy.deepcopy(self.get_site_config(w)[0]) for w in range(1, W)}

			# accumulate e-step stats across all training sequences 
			xi_acc = gamma_acc = site_xi_acc = site_gamma_acc = None
			total_log_prob = 0.0

			for seq_data in self.training_sequences: 
				self.Svec_sites = seq_data['Svec_sites']
				self.Svec_windows = seq_data['Svec_windows']
				self.positions = seq_data['positions']
				self.windows = seq_data['windows']

				xi, gamma, s_xi, s_gamma, lp = self.xi_gamma_hierarchical_window()
				total_log_prob += lp 
				T = len(gamma)

				if xi_acc is None: 
					xi_acc = xi
					gamma_acc = gamma
					site_xi_acc = s_xi
					site_gamma_acc = s_gamma

				else: 
					# window level
					for t in range(T - 1): 
						for w in range(1, W): 
							for w2 in range(1, W): 
								if xi[t][w][w2] != 'n/a': 
									xi_acc[t][w][w2] = (self.logsumexp([xi_acc[t][w][w2], xi[t][w][w2]]) if xi_acc[t][w][w2] != 'n/a' else xi[t][w][w2])
			
					for t in range(T): 
						for w in range (1, W): 
							if gamma[t][w] != 'n/a': 
								gamma_acc[t][w] = (self.logsumexp([gamma_acc[t][w], gamma[t][w]]) if gamma_acc[t][w] != 'n/a' else gamma[t][w])

					# site level
					for w in range(1, W): 
						_, site_states = self.get_site_config(w)
						N = len(site_states)
						for i in range(1, N): 
							if s_gamma[w][i] != 'n/a': 
								site_gamma_acc[w][i] = (self.logsumexp([site_gamma_acc[w][i], s_gamma[w][i]]) if site_gamma_acc[w][i] != 'n/a' else s_gamma[w][i])
							for j in range(1, N): 
								if s_xi[w][i][j] != 'n/a': 
									site_xi_acc[w][i][j] = (self.logsumexp([site_xi_acc[w][i][j], s_xi[w][i][j]]) if site_xi_acc[w][i][j] != 'n/a' else s_xi[w][i][j])

			# diagnostics
			# for w in range(1, W):
			# 	vals = [gamma_acc[t][w] for t in range(len(gamma_acc)) if gamma_acc[t][w] != 'n/a']
			# 	print(f"window state {w}: {len(vals)} defined gamma entries")

			# for w in range(1, W):
			# 	_, site_states = self.get_site_config(w)
			# 	N = len(site_states)
			# 	n_xi    = sum(1 for i in range(N) for j in range(N)
			# 				if site_xi_acc[w][i][j] != 'n/a')
			# 	n_gamma = sum(1 for i in range(N)
			# 				if site_gamma_acc[w][i] != 'n/a')
			# 	print(f"window state {w}: {n_xi} defined site xi, {n_gamma} defined site gamma")

			# n_xi_win = sum(1 for t in range(len(xi_acc))
			# 			for w in range(W) for w2 in range(W)
			# 			if xi_acc[t][w][w2] != 'n/a')
			# print(f"window xi defined entries: {n_xi_win}")

			# m-step
			self.maximization_hierarchical(xi_acc, gamma_acc, site_xi_acc, site_gamma_acc)

			# print("window_trans_mat sample:", 
			# 	[(w, w2, self.window_trans_mat[w][w2]) 
			# 	for w in range(W) for w2 in range(W) 
			# 	if self.window_trans_mat[w][w2] != 'n/a'][:6])
			
            # convergence check on both levels:
			delta_win = max(abs(self.window_trans_mat[w][w2] - prev_window[w][w2]) for w in range(1, W) for w2 in range(W))
			delta_site = max(abs(self.get_site_config(w)[0][i][j] - prev_site[w][i][j]) for w in range(1, W) for i in range(len(self.get_site_config(w)[1])) for j in range(len(self.get_site_config(w)[1])))
			delta = max(delta_win, delta_site)

			print(f"delta_window = {delta_win:.8f}")
			print(f"delta_site = {delta_site:.8f}")
			print(f"delta = {delta:.8f}")

			print(f"Iter {iteration + 1:>3d} | log P(O|λ) = {total_log_prob:.4f} | "f"max |Δ window| = {delta_win:.6f} | max |Δ site| = {delta_site:.6f}\n", flush = True)

			if delta < tol:
				print(f"Converged after {iteration + 1} iterations.", flush = True)
				break

		return self.window_trans_mat, self.site_trans_mat_neutral, self.site_trans_mat_linkedleft, self.site_trans_mat_sweep, self.site_trans_mat_linkedright

	def get_site_config(self, w, n_sites = None):
			
		# map window state index to its site-level trans matrix and state list
		if w == 1:
			return self.site_trans_mat_neutral, self.site_states_neutral_list
		elif w == 2:
			return self.site_trans_mat_linkedleft, self.site_states_linkedleft_list
		elif w == 3:
			m = self.site_trans_mat_sweep if n_sites is None else self.sweep_site_matrix(n_sites)
			return m, self.site_states_sweep_list
		else:
			return self.site_trans_mat_linkedright, self.site_states_linkedright_list

	def xi_gamma_site(self, Svec, position_indices, trans_mat, state_list, classifier_obj): 

		# get current site positions
		Svec_current = Svec[position_indices[0]:position_indices[-1] + 1]
		total_length = len(Svec_current)
		n_states = len(state_list)

		# get forward and backward probabilities
		logF, log_prob = self.forward(Svec, position_indices, trans_mat, state_list, classifier_obj)
		logB, _ = self.backward(Svec, position_indices, trans_mat, state_list, classifier_obj)

		# initialize xi and gamma data structures
		xi = [[['n/a'] * n_states for _ in range(n_states)] for _ in range(total_length - 1)]		
		gamma = [['n/a'] * n_states for _ in range(total_length)]
		
		# gamma recursion
		for t in range(1, total_length + 1): 
			for i in range(1, n_states): 
				if logF[t][i] != 'n/a' and logB[t][i] != 'n/a': 
					gamma[t - 1][i] = logF[t][i] + logB[t][i] - log_prob

		# xi recursion 
		for t in range(1, total_length): 
			for i in range(1, n_states): 
				if logF[t][i] == 'n/a': 
					continue
				for j in range(1, n_states): 
					if logB[t + 1][j] == 'n/a' or trans_mat[i][j] == 0: 
						continue
					emission = self.emission(state_list[j], Svec_current[t], classifier_obj)
					if emission is None or emission <= 0: 
						continue 
					xi[t - 1][i][j] = logF[t][i] + math.log(trans_mat[i][j]) + math.log(emission) + logB[t + 1][j] - log_prob

		return xi, gamma, log_prob

	def xi_gamma_window(self):
		"""
		window-level-only version: site transitions/emissions are fixed,
		so we only need the site-level forward pass (for log_prob),
		not the site-level backward/xi/gamma.
		"""

		empty_windows = []
		for t in range(len(self.windows)):
			win_start = self.windows[t]
			# win_end = self.windows[t] + self.window_size - 1
			# position_indices = np.where((np.array(self.positions) >= win_start) & (np.array(self.positions) <= win_end))
			win_end = self.windows[t] + self.window_size
			position_indices = self.get_nonempty_position_indices(win_start, win_end)
			if len(position_indices) == 0:
				print("WARNING: window " + str(win_start) + " does not contain any sites. Removing...")
				empty_windows.append(self.windows[t])

		indices_list = [index for index, item in enumerate(self.windows) if item in empty_windows]
		empty_windows = set(empty_windows)
		self.windows = [i for i in self.windows if i not in empty_windows]
		self.Svec_windows = [val for i, val in enumerate(self.Svec_windows) if i not in indices_list]

		n_window_intervals = len(self.Svec_windows)
		n_window_states = len(self.window_state_list)

		window_forward, log_prob = self.hierarchical_forward()
		window_backward, _ = self.hierarchical_backward()

		gamma = [['n/a'] * n_window_states for _ in range(n_window_intervals)]
		for t in range(1, n_window_intervals + 1):
			for w in range(1, n_window_states):
				if window_forward[t][w] != 'n/a' and window_backward[t][w] != 'n/a':
					gamma[t - 1][w] = window_forward[t][w] + window_backward[t][w] - log_prob

		# for t in range(0, n_window_intervals, max(1, n_window_intervals // 20)):
		# 	row = [math.exp(gamma[t][w]) if gamma[t][w] != 'n/a' else 0.0 for w in range(1, n_window_states)]
		# 	print(f"t={t}: {['%.3e' % x for x in row]}", flush = True)

		# only need the site-level marginal log-likelihood per window interval/state
		log_site_prob = [['n/a'] * n_window_states for _ in range(n_window_intervals + 2)]

		for t in range(1, n_window_intervals + 1):
			win_start = self.windows[t - 1]
			win_end = win_start + self.window_size - 1
			pos_win = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) < win_end)]
			position_indices = [self.positions.index(p) for p in pos_win]

			for w in range(1, n_window_states):
				site_trans, site_states = self.get_site_config(w)
				_, site_lp = self.forward(self.Svec_sites, position_indices, site_trans, site_states, self.classifierobj_sites)
				log_site_prob[t][w] = site_lp

		# window-level xi recursion
		xi = [[['n/a'] * n_window_states for _ in range(n_window_states)] for _ in range(n_window_intervals - 1)]
		for t in range(1, n_window_intervals):
			log_emit_next = ['n/a'] * n_window_states
			for w2 in range(1, n_window_states):
				if log_site_prob[t + 1][w2] == 'n/a':
					continue
				win_emission = self.emission(self.window_state_list[w2], self.Svec_windows[t], self.classifierobj_windows)
				if win_emission is not None and win_emission > 0:
					#log_emit_next[w2] = log_site_prob[t + 1][w2] + math.log(win_emission)
					log_emit_next[w2] = self.alpha * log_site_prob[t + 1][w2] + math.log(win_emission)

			for w in range(1, n_window_states):
				if window_forward[t][w] == 'n/a':
					continue
				for w2 in range(1, n_window_states):
					if window_backward[t + 1][w2] == 'n/a':
						continue
					if self.window_trans_mat[w][w2] == 0:
						continue
					if log_emit_next[w2] == 'n/a':
						continue
					xi[t - 1][w][w2] = (window_forward[t][w] + math.log(self.window_trans_mat[w][w2]) + log_emit_next[w2] + window_backward[t + 1][w2] - log_prob)

		return xi, gamma, log_prob

	def xi_gamma_hierarchical_window(self): # hierarchical version of xi_gamma_site

		# first remove windows that don't contain sites 
		empty_windows = []
		for t in range(len(self.windows)): 
			win_start = self.windows[t]
			win_end = self.windows[t] + self.window_size
			position_indices = self.get_nonempty_position_indices(win_start, win_end)
		
			if len(position_indices) == 0: 
				print("WARNING: window " + str(win_start) + " does not contain any sites. Removing...")
				empty_windows.append(self.windows[t])

		indices_list = [index for index, item in enumerate(self.windows) if item in empty_windows]

		empty_windows = set(empty_windows)
		self.windows = [i for i in self.windows if i not in empty_windows]
		self.Svec_windows = [val for i, val in enumerate(self.Svec_windows) if i not in indices_list]
		
		# get window-level forward and backward probabilities
		n_window_intervals = len(self.Svec_windows)
		n_window_states = len(self.window_state_list)

		window_forward, log_prob = self.hierarchical_forward()
		window_backward, _ = self.hierarchical_backward()

		# window-level gamma recursion
		gamma = [['n/a'] * n_window_states for _ in range(n_window_intervals)] # make gamma data structure 

		for t in range(1, n_window_intervals + 1): 
			for w in range(1, n_window_states): # exclude start/end state
				if window_forward[t][w] != 'n/a' and window_backward[t][w] != 'n/a': # check if t and w placement is right here !!!
					gamma[t - 1][w] = window_forward[t][w] + window_backward[t][w] - log_prob

		site_xi_acc = {}
		site_gamma_acc = {}
		for w in range(1, n_window_states): 
			_, site_states = self.get_site_config(w)
			N = len(site_states)
			site_xi_acc[w] = [['n/a'] * N for _ in range(N)]
			site_gamma_acc[w] = ['n/a'] * N

		log_site_prob = [['n/a'] * n_window_states for _ in range(n_window_intervals + 2)]

		for t in range(1, n_window_intervals + 1): 
			win_start = self.windows[t - 1]
			win_end = win_start + self.window_size - 1
			pos_win  = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) <  win_end)]
			position_indices = [self.positions.index(p) for p in pos_win]

			for w in range(1, n_window_states): 
				site_trans, site_states = self.get_site_config(w)

				# site-level forward and backward calculations 
				site_xi, site_gamma_t, site_lp = self.xi_gamma_site(self.Svec_sites, position_indices, site_trans, site_states, self.classifierobj_sites)

				log_site_prob[t][w] = site_lp

				if gamma[t - 1][w] == 'n/a': 
					continue

				gw = min(gamma[t - 1][w], 0.0) # clipped 
				N = len(site_states)
				T_site = len(site_gamma_t)

				# accumulate site xi and site gamma weighted by window xi and window gamma respectively
				for i in range(1, N): 
					for j in range(1, N): 
						vals = [site_xi[ts][i][j] for ts in range(T_site - 1) if site_xi[ts][i][j] != 'n/a']
						if vals: 
							weighted = gw + self.logsumexp(vals)
							if site_xi_acc[w][i][j] != 'n/a': 
								site_xi_acc[w][i][j] = self.logsumexp([site_xi_acc[w][i][j], weighted])
							else: 
								site_xi_acc[w][i][j] = weighted		

				for i in range(1, N): 
					vals = [site_gamma_t[ts][i] for ts in range(T_site - 1) if site_gamma_t[ts][i] != 'n/a']
					if vals: 
						weighted = gw + self.logsumexp(vals)
						if site_gamma_acc[w][i] != 'n/a': 
							site_gamma_acc[w][i] = self.logsumexp([site_gamma_acc[w][i], weighted])
						else: 
							site_gamma_acc[w][i] = weighted

		# check for bad gammas to see if there are overshoots (will print in log)
		for t in range(1, n_window_intervals + 1):
			for w in range(1, n_window_states):
				if gamma[t-1][w] != 'n/a' and gamma[t-1][w] > 1e-6:
					print(f"Warning: Bad gamma[{t - 1}][{w}] = {gamma[t - 1][w]:.4f}"
						f"fwd = {window_forward[t][w]:.4f}"
						f"bwd = {window_backward[t][w]:.4f}"
						f"log_prob = {log_prob:.4f}")

		# window-level xi recursion 
		xi = [[['n/a'] * n_window_states for _ in range(n_window_states)] for _ in range(n_window_intervals - 1)]

		for t in range(1, n_window_intervals): 
			log_emit_next = ['n/a'] * n_window_states
			for w2 in range(1, n_window_states): 
				if log_site_prob[t + 1][w2] == 'n/a': 
					continue
				win_emission = self.emission(self.window_state_list[w2], self.Svec_windows[t], self.classifierobj_windows)
				if win_emission is not None and win_emission > 0: 
					log_emit_next[w2] = log_site_prob[t + 1][w2] + math.log(win_emission)

			for w in range(1, n_window_states): 
				if window_forward[t][w] == 'n/a': 
					continue
				for w2 in range(1, n_window_states): 
					if window_backward[t + 1][w2] == 'n/a': 
						continue 
					if self.window_trans_mat[w][w2] == 0: 
						continue 
					if log_emit_next[w2] == 'n/a': 
						continue
					xi[t - 1][w][w2] = window_forward[t][w] + math.log(self.window_trans_mat[w][w2]) + log_emit_next[w2] + window_backward[t + 1][w2] - log_prob

		return xi, gamma, site_xi_acc, site_gamma_acc, log_prob

	def maximization_hierarchical(self, xi_acc, gamma_acc, site_xi_acc, site_gamma_acc): 
		
		# update window_trans_mat and all four site-level trans mats
		# emissions are pre-trained and stay fixed

		W = len(self.window_state_list)
		T = len(gamma_acc)

		# window-level update 
		xi_sum = [['n/a'] * W for _ in range(W)]
		gamma_sum = ['n/a'] * W

		for w in range(1, W): 
			g_vals = [gamma_acc[t][w] for t in range(T - 1) if gamma_acc[t][w] != 'n/a']
			if g_vals: 
				gamma_sum[w] = self.logsumexp(g_vals)
			for w2 in range(1, W): 
				x_vals = [xi_acc[t][w][w2] for t in range(T - 1) if xi_acc[t][w][w2] != 'n/a']
				if x_vals: 
					xi_sum[w][w2] = self.logsumexp(x_vals)

		for w in range(1, W): 
			if gamma_sum[w] == 'n/a':
				continue
			for w2 in range(1, W): 
				if xi_sum[w][w2] != 'n/a': 
					self.window_trans_mat[w][w2] = math.exp(xi_sum[w][w2] - gamma_sum[w])
				else: 
					self.window_trans_mat[w][w2] = 0.0

			# clipping probabilities so they end up between [0:1]
			row = self.window_trans_mat[w][1:W]
			now = np.clip(row, 0.0, 1.0)
			row_sum = row.sum()
			if row_sum > 0: 
				row = row / row_sum
			self.window_trans_mat[w][1:W] = row
		
		# site-level update
		for w in range(1, W): 
			site_trans, site_states = self.get_site_config(w)
			N = len(site_states)

			for i in range(1, N): 
				if site_gamma_acc[w][i] == 'n/a': 
					continue
				for j in range(1, N): 
					if site_xi_acc[w][i][j] != 'n/a':
						site_trans[i][j] =  math.exp(site_xi_acc[w][i][j] - site_gamma_acc[w][i])
					else: 
						site_trans[i][j] = 0.0

		for w in range(1, W):
				for w2 in range(1, W):
					v = self.window_trans_mat[w][w2]
					if v != 'n/a' and v > 1.0:
						print(f"BAD: window_trans_mat[{w}][{w2}] = {v}  (xi = {xi_sum[w][w2]:.4f}, gamma = {gamma_sum[w]:.4f})")


		return self.window_trans_mat, self.site_trans_mat_neutral, self.site_trans_mat_linkedleft, self.site_trans_mat_sweep, self.site_trans_mat_linkedright

	def emission(self, scenario, S, classifierobj): 

		# call emission probabilities from classifiers.py 
		if self.emissiontype == 'NB':
			Likelihood = classifierobj.nb_likelihood(scenario, S)
		else:
			Likelihood = classifierobj.aode_likelihood(scenario, S)
		try:
			if Likelihood < 1e-200:
				Likelihood = 1e-200
			return Likelihood
		
		except TypeError: 
			print("ERROR: one or more input rows have no defined statistics. Please remove from your data or recalculate.")
			sys.exit(0)

	def logsumexp(self, sumvec): 
		
        # sumvec is log of entries in the sum (x_i's in log-sum-exp)
		# e.g. logV[i-1][j']+log(t[j'][j])
		if len(sumvec) == 0:
			return 'n/a'
		
		a = max(sumvec)
		exponents = [x - a for x in sumvec]
		newsum = sum([math.exp(x) for x in exponents])
		newlogsum = math.log(newsum)
		
		return a + newlogsum

	def viterbi(self, Svec, position_indices, trans_mat, state_list, classifierobj): 

		Svec_current = Svec[position_indices[0]:position_indices[-1]+1]
		pos_output = self.positions[position_indices[0]:position_indices[-1]+1]
		total_length = len(Svec_current)

		logV = [['n/a' for _ in range(len(state_list))] for _ in range(total_length + 2)] 
		P = [['n/a' for _ in range(len(state_list))] for _ in range(total_length + 2)]

		# initialization
		logV[0][0] = 0 # no emission

		# recurrence 
		for i in range(1, total_length + 1):
			for j in range(1, len(state_list)):
				maxpathprobs = -1*np.inf
				for k in range(len(state_list)):
						if logV[i-1][k] != 'n/a' and trans_mat[k][j] != 0:
							log_prob = float(logV[i-1][k]) + float(math.log(trans_mat[k][j]))
							if log_prob > maxpathprobs:
								maxpathprobs = log_prob
								argmax = k 

				if maxpathprobs == -1*np.inf:
					next
				else:
					logV[i][j] = maxpathprobs + math.log(self.emission(state_list[j], Svec_current[i-1], classifierobj))
					P[i][j] = argmax

		# end state
		maxpathprobs = -1*np.inf
		for k in range(len(state_list)):
			if logV[total_length][k] != 'n/a' and trans_mat[k][0] != 0: 
				logpp = float(logV[total_length][k]) + float(math.log(trans_mat[k][0]))
				if logpp > maxpathprobs:
					maxpathprobs = logpp
					argmax = k
		logV[total_length + 1][0] = maxpathprobs
		P[total_length + 1][0] = argmax

		# termination/backtracking
		state_path = [0 for _ in range(total_length + 2)]
		state_path[total_length] = P[total_length + 1][0]
		for i in range(total_length -1,-1,-1):
			state_path[i] = P[i+1][state_path[i+1]]
		
		return logV[total_length + 1][0], state_path, pos_output

	def hierarchical_viterbi(self, windows_out, sites_out):  

		# check if there are sites in every window before running recursion, remove empty windows and throw warning
		empty_windows = []
		for t in range(len(self.windows)): 
			win_start = self.windows[t]
			win_end = self.windows[t] + self.window_size
			position_indices = self.get_nonempty_position_indices(win_start, win_end)
		
			# if len(position_indices) == 0: 
			# 	print("WARNING: window " + str(win_start) + " does not contain any sites. Removing...")
			# 	empty_windows.append(self.windows[t])

			if len(position_indices) < 2: 
				print("WARNING: window " + str(win_start) + " has fewer than 2 sites. Removing...")
				empty_windows.append(self.windows[t])

		indices_list = [index for index, item in enumerate(self.windows) if item in empty_windows]

		empty_windows = set(empty_windows)
		self.windows = [i for i in self.windows if i not in empty_windows]
		self.Svec_windows = [val for i, val in enumerate(self.Svec_windows) if i not in indices_list]
		
		# initialize data structures 
		n_window_states = len(self.window_state_list)
		n_window_intervals = len(self.windows)
		
		window_viterbi = [['n/a' for _ in range(n_window_states)] for _ in range(n_window_intervals + 2)] 
		window_backpointer = [['n/a' for _ in range(n_window_states)] for _ in range(n_window_intervals + 2)]

		window_viterbi[0][0] = 0 # no emission for start/end state
		position_indices = None

		# intialize site-level data structure for current window
		site_viterbi_paths = [[[] for _ in range(n_window_states)] for _ in range(n_window_intervals + 2)]
		site_positions = []
		empty_windows = []
		margin_log = []

		# recursion over all windows
		print("Performing hierarchical Viterbi algorithm...", flush = True)
		for t in range(1, n_window_intervals + 1):
			win_start = self.windows[t-1]
			win_end = self.windows[t-1] + self.window_size - 1

			# check if there is a gap between the previous window and current window, and report gap size/throw warning	
			# if t > 1:
			# 	if (self.windows[t - 1] - self.windows[t - 2]) > self.window_size: 
			# 		print("WARNING: Gap between windows of size " + str(self.windows[t - 1] - self.windows[t - 2]))

			# get position indices 
			position_indices = np.where((np.array(self.positions) >= win_start) & (np.array(self.positions) <= win_end))
			position_indices = list(position_indices[0])
			
			positions = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) <= win_end)]
			site_positions.append(positions)

			for w in range(1, n_window_states):
				if w == 1: 
					max_prob = float('-inf')
					best_prev_state = -1

					site_log_probs, site_path, site_pos = self.viterbi(self.Svec_sites, position_indices, self.site_trans_mat_neutral, self.site_states_neutral_list, self.classifierobj_sites)

					# ### check
					# if win_start == 1480001:
					# 	we = math.log(self.emission(self.window_state_list[w],
					# 								self.Svec_windows[t-1],
					# 								self.classifierobj_windows))
					# 	print(f"[focal] state={self.window_state_list[w]:8s} "
					# 		f"site_lp={site_log_probs:10.1f}  win_emit={we:7.2f}  "
					# 		f"n_sites={len(position_indices)}", flush=True)
					
					## margins check 
					# we = math.log(self.emission(self.window_state_list[w],
                    #         self.Svec_windows[t-1],
                    #         self.classifierobj_windows))
					# margin_log.append((win_start, w, site_log_probs, we))

					for prev_state in range(n_window_states):
						if window_viterbi[t-1][prev_state] != 'n/a' and self.window_trans_mat[prev_state][w] != 0:
							prob = (float(window_viterbi[t-1][prev_state]) + float(np.log(self.window_trans_mat[prev_state, w])) + self.alpha * float(site_log_probs))
							if prob > max_prob:
								max_prob = prob
								best_prev_state = prev_state

					if max_prob == float('-inf'):
						next
					else:
						window_backpointer[t][w] = best_prev_state
						window_viterbi[t][w] = max_prob + math.log(self.emission(self.window_state_list[w], self.Svec_windows[t-1], self.classifierobj_windows)) 
						
						site_viterbi_paths[t][w] = site_path 

				elif w == 2: 
					max_prob = float('-inf')
					best_prev_state = -1
 
					site_log_probs, site_path, site_pos = self.viterbi(self.Svec_sites, position_indices, self.site_trans_mat_linkedleft, self.site_states_linkedleft_list, self.classifierobj_sites)

					# ### check
					# if win_start == 1480001:
					# 	we = math.log(self.emission(self.window_state_list[w],
					# 								self.Svec_windows[t-1],
					# 								self.classifierobj_windows))
					# 	print(f"[focal] state={self.window_state_list[w]:8s} "
					# 		f"site_lp={site_log_probs:10.1f}  win_emit={we:7.2f}  "
					# 		f"n_sites={len(position_indices)}", flush=True)

					## margins check 
					# we = math.log(self.emission(self.window_state_list[w],
                    #         self.Svec_windows[t-1],
                    #         self.classifierobj_windows))
					# margin_log.append((win_start, w, site_log_probs, we))

					for prev_state in range(n_window_states):
						if window_viterbi[t-1][prev_state] != 'n/a' and self.window_trans_mat[prev_state][w] != 0:
							prob = (float(window_viterbi[t-1][prev_state]) + float(np.log(self.window_trans_mat[prev_state, w])) + self.alpha * float(site_log_probs))
							if prob > max_prob:
								max_prob = prob
								best_prev_state = prev_state

					if max_prob == float('-inf'):
						next
					else:
						window_backpointer[t][w] = best_prev_state
						window_viterbi[t][w] = max_prob + math.log(self.emission(self.window_state_list[w], self.Svec_windows[t-1], self.classifierobj_windows)) 
						
						site_viterbi_paths[t][w] = site_path

				elif w == 3: 
					max_prob = float('-inf')
					best_prev_state = -1

					site_log_probs, site_path, site_pos = self.viterbi(self.Svec_sites, position_indices, self.sweep_site_matrix(len(position_indices)), self.site_states_sweep_list, self.classifierobj_sites)

					# ### check
					# if win_start == 1480001:
					# 	we = math.log(self.emission(self.window_state_list[w],
					# 								self.Svec_windows[t-1],
					# 								self.classifierobj_windows))
					# 	print(f"[focal] state={self.window_state_list[w]:8s} "
					# 		f"site_lp={site_log_probs:10.1f}  win_emit={we:7.2f}  "
					# 		f"n_sites={len(position_indices)}", flush=True)

					## margins check 
					# we = math.log(self.emission(self.window_state_list[w],
                    #         self.Svec_windows[t-1],
                    #         self.classifierobj_windows))
					# margin_log.append((win_start, w, site_log_probs, we))

					for prev_state in range(n_window_states):
						if window_viterbi[t-1][prev_state] != 'n/a' and self.window_trans_mat[prev_state][w] != 0:
							prob = (float(window_viterbi[t-1][prev_state]) + float(np.log(self.window_trans_mat[prev_state, w])) + self.alpha * float(site_log_probs))
							if prob > max_prob:
								max_prob = prob
								best_prev_state = prev_state

					if max_prob == float('-inf'):
						next
					else:
						window_backpointer[t][w] = best_prev_state
						window_viterbi[t][w] = max_prob + math.log(self.emission(self.window_state_list[w], self.Svec_windows[t-1], self.classifierobj_windows)) 
						
						site_viterbi_paths[t][w] = site_path

				else: 
					max_prob = float('-inf')
					best_prev_state = -1

					site_log_probs, site_path, site_pos = self.viterbi(self.Svec_sites, position_indices, self.site_trans_mat_linkedright, self.site_states_linkedright_list, self.classifierobj_sites)

					# ### check
					# if win_start == 1480001:
					# 	we = math.log(self.emission(self.window_state_list[w],
					# 								self.Svec_windows[t-1],
					# 								self.classifierobj_windows))
					# 	print(f"[focal] state={self.window_state_list[w]:8s} "
					# 		f"site_lp={site_log_probs:10.1f}  win_emit={we:7.2f}  "
					# 		f"n_sites={len(position_indices)}", flush=True)

					for prev_state in range(n_window_states):
						if window_viterbi[t-1][prev_state] != 'n/a' and self.window_trans_mat[prev_state][w] != 0:
							prob = (float(window_viterbi[t-1][prev_state]) + float(np.log(self.window_trans_mat[prev_state, w])) + self.alpha * float(site_log_probs))
							if prob > max_prob:
								max_prob = prob
								best_prev_state = prev_state


					if max_prob == float('-inf'):
						next
					else:
						window_backpointer[t][w] = best_prev_state
						window_viterbi[t][w] = max_prob + math.log(self.emission(self.window_state_list[w], self.Svec_windows[t-1], self.classifierobj_windows)) 
						
						site_viterbi_paths[t][w] = site_path
		
		#print(window_viterbi)
		#print('\n')

		# window end state 
		maxpathprobs = -1*np.inf
		for w in range(n_window_states):
			if window_viterbi[len(self.windows)][w] != 'n/a' and self.window_trans_mat[w][0] != 0:
				logpp = float(window_viterbi[len(self.windows)][w]) + float(math.log(self.window_trans_mat[w][0]))
				if logpp > maxpathprobs:
					maxpathprobs = logpp
					argmax = w
		window_viterbi[len(self.windows) + 1][0] = maxpathprobs
		window_backpointer[len(self.windows) + 1][0] = argmax


		# backtracking
		window_paths = np.zeros(n_window_intervals + 2, dtype=int)
		window_paths[len(self.windows)] = window_backpointer[len(self.windows) + 1][0]
		for t in range(n_window_intervals -1, -1, -1):
			window_paths[t] = window_backpointer[t + 1][window_paths[t + 1]]

		site_paths = []
		counter = 1 
		for w in window_paths: 
			if w != 0: 
				if counter < n_window_intervals + 1:
					site_paths.append(site_viterbi_paths[counter][w])
					counter += 1

		# print output file
		print("Writing site and window path files...", flush = True)
		site_results = []
		window_results = []
		window_list = []

		for j in range(1, len(window_paths)-1): 
			window_results.append(window_paths[j])
			window_list.append(self.windows[j-1])

			site_path_without_zeros = site_paths[j-1][1:-1]
			site_results.append(site_path_without_zeros)

		site_pos_flat = [item for sublist in site_positions for item in sublist]
		site_results_flat = [item for sublist in site_results for item in sublist]

		site_output = open(sites_out, 'w')
		site_output.write("site\tsite_pred\n")
		for s in range(len(site_results_flat)): 
			site_output.write('\t'.join([str(site_pos_flat[s]), str(site_results_flat[s]) + ("\n")]))
	
		window_output = open(windows_out, 'w')
		window_output.write("win\tpred\n")
		for i in range(len(window_results)):
			window_output.write('\t'.join([str(window_list[i]), str(window_results[i]) + ("\n")]))

		# rep = os.path.basename(windows_out).replace(".window_results", "")
		# outdir = "windex_results/stdpopsim_strong_selection/training_margins"
		# os.makedirs(outdir, exist_ok=True)
		# with open(os.path.join(outdir, f"margins_{rep}.txt"), "w") as fh:
		# 	for win, w, slp, we in margin_log:
		# 		fh.write(f"{rep}\t{win}\t{w}\t{slp}\t{we}\n")
		
		print("Done.\n")

	def forward(self, Svec, position_indices, trans_mat, state_list, classifierobj):

		Svec_current = Svec[position_indices[0]:position_indices[-1]+1]
		logV = [['n/a' for _ in range(len(state_list))] for _ in range(len(Svec_current) + 2)] 
		logV[0][0] = 0

		# recursion	
		for i in range(1, len(Svec_current) + 1): 
			for j in range(1, len(state_list)):
					sumvec = []
					for k in range(len(state_list)):
						if logV[i-1][k] != 'n/a' and trans_mat[k][j] != 0:
							sumvec.append(logV[i-1][k]+math.log(trans_mat[k][j]))
					
					if len(sumvec) > 0:
						Sum = self.logsumexp(sumvec)
						logV[i][j] = Sum + math.log(self.emission(state_list[j], Svec_current[i-1], classifierobj))
						
		Endsumvec = []
		for k in range(1, len(state_list)):
			if logV[len(Svec_current)][k] != 'n/a' and trans_mat[k][0] != 0:
				Endsumvec.append(float(logV[len(Svec_current)][k]) + float(math.log(trans_mat[k][0])))

		EndSum = self.logsumexp(Endsumvec)
		logV[len(Svec_current)+1][0] = EndSum
		
		return logV, EndSum

	def hierarchical_forward(self):

		# check if there are sites in every window before running recursion, remove empty window and throw warning 
		empty_windows = []
		for t in range(len(self.windows)): 
			win_start = self.windows[t]
			win_end = self.windows[t] + self.window_size
			position_indices = self.get_nonempty_position_indices(win_start, win_end)
		
			# if len(position_indices) == 0: 
			# 	print("WARNING: window " + str(win_start) + " does not contain any sites. Removing...")
			# 	empty_windows.append(self.windows[t])

			if len(position_indices) < 2: 
				print("WARNING: window " + str(win_start) + " has fewer than 2 sites. Removing...")
				empty_windows.append(self.windows[t])

		indices_list = [index for index, item in enumerate(self.windows) if item in empty_windows]

		empty_windows = set(empty_windows)
		self.windows = [i for i in self.windows if i not in empty_windows]
		self.Svec_windows = [val for i, val in enumerate(self.Svec_windows) if i not in indices_list]
		
		# initialize data structures
		n_window_states = len(self.window_state_list)
		n_window_intervals = len(self.Svec_windows)

		window_forward = [['n/a' for _ in range(n_window_states)] for _ in range(n_window_intervals + 2)] 

		window_forward[0][0] = 0 # no emission

		# recursion 
		for t in range(1, n_window_intervals + 1):

			win_start = self.windows[t-1]
			win_end = win_start + self.window_size

			pos_window = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) < win_end)]
			position_indices = [self.positions.index(pos) for pos in pos_window]

			for w in range(1, n_window_states):

				site_trans, site_states = self.get_site_config(w, n_sites = len(position_indices))
				_, site_log_probsfinal = self.forward(self.Svec_sites, position_indices, site_trans, site_states, self.classifierobj_sites)

				if site_log_probsfinal == 'n/a':
					continue

				sumvec = []
				for k in range(n_window_states):
					if window_forward[t-1][k] != 'n/a' and self.window_trans_mat[k][w] != 0:
						sumvec.append(float(window_forward[t-1][k]) + float(math.log(self.window_trans_mat[k][w])) + self.alpha * float(site_log_probsfinal))

				if len(sumvec) > 0:
					Sum = self.logsumexp(sumvec)
					window_forward[t][w] = Sum + math.log(self.emission(self.window_state_list[w], self.Svec_windows[t-1], self.classifierobj_windows))

		Endsumvec = []
		for k in range(1, len(self.window_state_list)):
			if window_forward[n_window_intervals][k] != 'n/a' and self.window_trans_mat[k][0] != 0:
				Endsumvec.append(float(window_forward[n_window_intervals][k]) + float(math.log(self.window_trans_mat[k][0])))

		EndSum = self.logsumexp(Endsumvec)
		window_forward[n_window_intervals + 1][0] = EndSum

		return window_forward, EndSum

	def backward(self, Svec, position_indices, trans_mat, state_list, classifierobj): 

		Svec_current = Svec[position_indices[0]:position_indices[-1]+1]
		logV = [['n/a' for _ in range(len(state_list))] for _ in range(len(Svec_current) + 2)] 
		
		# initialization
		for s in range(0, len(state_list)): 
			if s != 0 and trans_mat[s][0] != 0: 
				logV[len(Svec_current)][s] = math.log(trans_mat[s][0])

		# recursion 
		for i in range(len(Svec_current) - 1, 0, -1): 
			for j in range(1, len(state_list)):
				sumvec = [] 
				for k in range(len(state_list)): 
					if logV[i + 1][k] != 'n/a' and trans_mat[j][k] != 0: 
						sumvec.append(logV[i + 1][k] +  math.log(trans_mat[j][k]) + math.log(self.emission(state_list[k], Svec_current[i], classifierobj)))

				if len(sumvec) > 0: 
					Sum = self.logsumexp(sumvec)
					logV[i][j] = Sum

		Beginsumvec = [] 
		for k in range(len(state_list)): 
			if logV[1][k] != 'n/a' and trans_mat[0][k] != 0: 
				Beginsumvec.append(logV[1][k] + math.log(trans_mat[0][k]) + math.log(self.emission(state_list[k], Svec_current[0], classifierobj)))
		
		BeginSum = self.logsumexp(Beginsumvec)
		logV[len(Svec_current) + 1][0] = BeginSum

		return logV, BeginSum

	def hierarchical_backward(self): 

		# check if there are sites in every window before running recursion, remove empty window and throw warning 
		empty_windows = []
		for t in range(len(self.windows)): 
			win_start = self.windows[t]
			win_end = self.windows[t] + self.window_size
			position_indices = self.get_nonempty_position_indices(win_start, win_end)
		
			# if len(position_indices) == 0: 
			# 	print("WARNING: window " + str(win_start) + " does not contain any sites. Removing...")
			# 	empty_windows.append(self.windows[t])

			if len(position_indices) < 2: 
				print("WARNING: window " + str(win_start) + " has fewer than 2 sites. Removing...")
				empty_windows.append(self.windows[t])

		indices_list = [index for index, item in enumerate(self.windows) if item in empty_windows]

		empty_windows = set(empty_windows)
		self.windows = [i for i in self.windows if i not in empty_windows]
		self.Svec_windows = [val for i, val in enumerate(self.Svec_windows) if i not in indices_list]
		
		# initialize data structures 
		n_window_states = len(self.window_state_list)
		n_window_intervals = len(self.windows)

		window_state_path = [0 for _ in range(n_window_intervals + 2)]
		window_backward = [['n/a' for _ in range(n_window_states)] for _ in range(n_window_intervals + 2)]
		
		position_indices = None
		site_paths = [[[] for _ in range(n_window_states)] for _ in range(n_window_intervals + 2)]
		site_positions = []

		# window-level backward algorithm initialization
		for s in range(n_window_states): 
			if s != 0 and self.window_trans_mat[s][0] != 0: 
				window_backward[len(self.Svec_windows)][s] = math.log(self.window_trans_mat[s][0])

		# recursion
		for t in range(n_window_intervals - 1, 0, -1):

			win_start = self.windows[t]
			win_end = win_start + self.window_size

			pos_window = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) < win_end)]
			position_indices = [self.positions.index(pos) for pos in pos_window]

			# compute each site-level marginal once per window, not once per (w, k) pair
			site_lp = ['n/a'] * n_window_states
			for k in range(1, n_window_states):
				site_trans, site_states = self.get_site_config(k, n_sites = len(position_indices))
				_, site_lp[k] = self.backward(self.Svec_sites, position_indices, site_trans, site_states, self.classifierobj_sites)

			for w in range(1, n_window_states):
				sumvec = []
				for k in range(1, n_window_states):
					if window_backward[t + 1][k] == 'n/a' or self.window_trans_mat[w][k] == 0:
						continue
					if site_lp[k] == 'n/a':
						continue

					win_emission = self.emission(self.window_state_list[k], self.Svec_windows[t], self.classifierobj_windows)
					if win_emission is None or win_emission <= 0:
						continue

					sumvec.append(float(window_backward[t + 1][k]) + float(math.log(self.window_trans_mat[w][k])) + self.alpha * float(site_lp[k]) + math.log(win_emission))

				if len(sumvec) > 0:
					window_backward[t][w] = self.logsumexp(sumvec)

		win_start = self.windows[0]
		win_end = win_start + self.window_size

		pos_window = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) < win_end)]
		position_indices = [self.positions.index(pos) for pos in pos_window]

		Beginsumvec = []
		for k in range(n_window_states):
			if window_backward[1][k] != 'n/a' and self.window_trans_mat[0][k] != 0:

				site_trans, site_states = self.get_site_config(k, n_sites = len(position_indices))
				_, site_log_probsfinal = self.backward(self.Svec_sites, position_indices, site_trans, site_states, self.classifierobj_sites)

				if site_log_probsfinal == 'n/a':
					continue

				Beginsumvec.append(window_backward[1][k] + math.log(self.window_trans_mat[0][k]) + self.alpha * float(site_log_probsfinal) + math.log(self.emission(self.window_state_list[k], self.Svec_windows[0], self.classifierobj_windows)))


		BeginSum = self.logsumexp(Beginsumvec)
		window_backward[len(self.Svec_windows) + 1][0] = BeginSum 

		return window_backward, BeginSum
		
	def stochastic_backtrace(self, Svec, logV, positions, position_indices, trans_mat, state_list):
		
		# initialize and create Probs vector
		Svec_current = Svec[position_indices[0]:position_indices[-1]+1]
		pos_output = positions[position_indices[0]:position_indices[-1]+1]

		state_path = [0 for _ in range(len(Svec_current) + 2)]
		unnormalized_logprobs = ['n/a' for _ in range(len(state_list))]
		
		for jprime in range(len(state_list)):
			if logV[len(Svec_current)][jprime] != 'n/a' and trans_mat[jprime][0] != 0:
				unnormalized_logprobs[jprime] = logV[len(Svec_current)][jprime] + math.log(trans_mat[jprime][0]) 
		normalization = self.logsumexp([x for x in unnormalized_logprobs if x != 'n/a'])
		normalized_logprobs = ['n/a' for _ in range(len(state_list))]
		for jprime in range(len(state_list)):
			if unnormalized_logprobs[jprime] != 'n/a' and trans_mat[jprime][0] != 0:
				normalized_logprobs[jprime] = logV[len(Svec_current)][jprime] + math.log(trans_mat[jprime][0]) -normalization 
		Probs = [0 for _ in range(len(state_list))]
		for jprime in range(len(state_list)):
			if normalized_logprobs[jprime] != 'n/a':
				Probs[jprime] = math.exp(normalized_logprobs[jprime])

		# randomly select the state at the last entry
		state = np.random.choice(len(state_list),p=Probs)
		state_path[len(Svec_current)] = state

		# start from last state and iterate backwards
		for i in range(len(Svec_current)-1,-1,-1):
			unnormalized_logprobs = ['n/a' for _ in range(len(state_list))]
			for jprime in range(len(state_list)):
				if logV[i][jprime] != 'n/a' and trans_mat[jprime][state_path[i+1]] != 0:
					unnormalized_logprobs[jprime] = logV[i][jprime] + math.log(trans_mat[jprime][state_path[i+1]]) 
			normalization = self.logsumexp([x for x in unnormalized_logprobs if x != 'n/a'])
			normalized_logprobs = ['n/a' for _ in range(len(state_list))]
			for jprime in range(len(state_list)):
				if unnormalized_logprobs[jprime] != 'n/a' and trans_mat[jprime][state_path[i+1]]:
					normalized_logprobs[jprime] = logV[i][jprime] + math.log(trans_mat[jprime][state_path[i+1]]) -normalization 
			Probs = [0 for _ in range(len(state_list))]
			for jprime in range(len(state_list)):
				if normalized_logprobs[jprime] != 'n/a':
					Probs[jprime] = math.exp(normalized_logprobs[jprime])

			state_path[i] = np.random.choice(len(state_list),p=Probs)

		return state_path, pos_output

	def hierarchical_stochastic_backtrace(self, logV_window): 

		# check if there are sites in every window before running recursion, remove empty window and throw warning 
		empty_windows = []
		for t in range(len(self.windows)): 
			win_start = self.windows[t]
			win_end = self.windows[t] + self.window_size
			position_indices = self.get_nonempty_position_indices(win_start, win_end)
		
			if len(position_indices) == 0: 
				print("WARNING: window " + str(win_start) + " does not contain any sites. Removing...")
				empty_windows.append(self.windows[t])

		indices_list = [index for index, item in enumerate(self.windows) if item in empty_windows]

		empty_windows = set(empty_windows)
		self.windows = [i for i in self.windows if i not in empty_windows]
		self.Svec_windows = [val for i, val in enumerate(self.Svec_windows) if i not in indices_list]
		
		# initialize data structures 
		n_window_states = len(self.window_state_list)
		n_window_intervals = len(self.windows)

		window_state_path = [0 for _ in range(n_window_intervals + 2)]
		unnormalized_logprobs = ['n/a' for _ in range(len(self.window_state_list))]
		position_indices = None
		site_paths = [[[] for _ in range(n_window_states)] for _ in range(n_window_intervals + 2)]
		site_positions = []

		# initialization of unnormalized logprobs for first window state
		for jprime in range(len(self.window_state_list)):
			if logV_window[len(self.Svec_windows)][jprime] != 'n/a' and self.window_trans_mat[jprime][0] != 0:
				unnormalized_logprobs[jprime] = logV_window[len(self.Svec_windows)][jprime] + math.log(self.window_trans_mat[jprime][0]) 
				
		normalization = self.logsumexp([x for x in unnormalized_logprobs if x != 'n/a'])
		normalized_logprobs = ['n/a' for _ in range(len(self.window_state_list))]
			
		for jprime in range(len(self.window_state_list)):
			if unnormalized_logprobs[jprime] != 'n/a' and self.window_trans_mat[jprime][0] != 0:
				normalized_logprobs[jprime] = logV_window[len(self.Svec_windows)][jprime] + math.log(self.window_trans_mat[jprime][0]) - normalization 
		
		Probs = [0 for _ in range(len(self.window_state_list))]
		for jprime in range(len(self.window_state_list)):
			if normalized_logprobs[jprime] != 'n/a':
				Probs[jprime] = math.exp(normalized_logprobs[jprime])
		
		state = np.random.choice(len(self.window_state_list),p = Probs)
		window_state_path[len(self.Svec_windows)] = state

		# get site paths of last initialized window state 
		win_start = self.windows[len(self.Svec_windows)-1]
		win_end = win_start + self.window_size -1

		# keep track of positions for sites 
		pos_window = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) < win_end)]
		position_indices = [self.positions.index(pos) for pos in pos_window]

		if state == 1: 
			# run forward algorithm to get logV for site-based: 
			logV_site, _ = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_neutral, self.site_states_neutral_list, self.classifierobj_sites)
			
			# run stochastic backtrace with logV_site to get path: 
			sb_site_path, pos_for_path = self.stochastic_backtrace(self.Svec_sites, logV_site, self.positions, position_indices, self.site_trans_mat_neutral, self.site_states_neutral_list)
			site_paths.append(sb_site_path) # add site path in to site path data structure
			site_positions.append(pos_for_path)

		elif state == 2: 
			# run forward algorithm to get logV for site-based: 
			logV_site, _ = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_linkedleft, self.site_states_linkedleft_list, self.classifierobj_sites)
			
			# run stochastic backtrace with logV_site to get path: 
			sb_site_path, pos_for_path = self.stochastic_backtrace(self.Svec_sites, logV_site, self.positions, position_indices, self.site_trans_mat_linkedleft, self.site_states_linkedleft_list)
			site_paths.append(sb_site_path) # add site path in to site path data structure
			site_positions.append(pos_for_path)

		elif state == 3: 
			# run forward algorithm to get logV for site-based: 
			logV_site, _ = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_sweep, self.site_states_sweep_list, self.classifierobj_sites)
			
			# run stochastic backtrace with logV_site to get path: 
			sb_site_path, pos_for_path = self.stochastic_backtrace(self.Svec_sites, logV_site, self.positions, position_indices, self.site_trans_mat_sweep, self.site_states_sweep_list)
			site_paths.append(sb_site_path) # add site path in to site path data structure
			site_positions.append(pos_for_path)

		elif state == 4: 
			# run forward algorithm to get logV for site-based: 
			logV_site, _ = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_linkedright, self.site_states_linkedright_list, self.classifierobj_sites)
			
			# run stochastic backtrace with logV_site to get path: 
			sb_site_path, pos_for_path = self.stochastic_backtrace(self.Svec_sites, logV_site, self.positions, position_indices, self.site_trans_mat_linkedright, self.site_states_linkedright_list)
			site_paths.append(sb_site_path) # add site path in to site path data structure
			site_positions.append(pos_for_path)

		# recursion; move backwards through the rest of the states, getting site paths per randomly chosen window
		for i in range(len(self.Svec_windows)-1,-1,-1):
			unnormalized_logprobs = ['n/a' for _ in range(len(self.window_state_list))]
			for jprime in range(len(self.window_state_list)):
				if logV_window[i][jprime] != 'n/a' and self.window_trans_mat[jprime][window_state_path[i+1]] != 0:
					unnormalized_logprobs[jprime] = logV_window[i][jprime] + math.log(self.window_trans_mat[jprime][window_state_path[i+1]]) 

			normalization = self.logsumexp([x for x in unnormalized_logprobs if x != 'n/a'])
			normalized_logprobs = ['n/a' for _ in range(len(self.window_state_list))]
			for jprime in range(len(self.window_state_list)):
				if unnormalized_logprobs[jprime] != 'n/a' and self.window_trans_mat[jprime][window_state_path[i+1]]:
					normalized_logprobs[jprime] = logV_window[i][jprime] + math.log(self.window_trans_mat[jprime][window_state_path[i+1]]) -normalization 
			Probs = [0 for _ in range(len(self.window_state_list))]
			for jprime in range(len(self.window_state_list)):
				if normalized_logprobs[jprime] != 'n/a':
					Probs[jprime] = math.exp(normalized_logprobs[jprime])

			state = np.random.choice(len(self.window_state_list), p = Probs)
			window_state_path[i] = state
		
			# get position indices of last window
			win_start = self.windows[i]
			win_end = win_start + self.window_size -1

			# get position indices 
			position_indices = np.where((np.array(self.positions) >= win_start) & (np.array(self.positions) <= win_end))
			position_indices = list(position_indices[0])

			if state == 1: 
				# run forward algorithm to get logV for site-based: 
				logV_site, _ = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_neutral, self.site_states_neutral_list, self.classifierobj_sites,)
				
				# run stochastic backtrace with logV_site to get path: 
				sb_site_path, pos_for_path = self.stochastic_backtrace(self.Svec_sites, logV_site, self.positions, position_indices, self.site_trans_mat_neutral, self.site_states_neutral_list)
				site_paths.append(sb_site_path) # add site path in as last entry to site path data structure
				site_positions.append(pos_for_path)

			elif state == 2: 
				# run forward algorithm to get logV for site-based: 
				logV_site, _ = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_linkedleft, self.site_states_linkedleft_list, self.classifierobj_sites)
				
				# run stochastic backtrace with logV_site to get path: 
				sb_site_path, pos_for_path = self.stochastic_backtrace(self.Svec_sites, logV_site, self.positions, position_indices, self.site_trans_mat_linkedleft, self.site_states_linkedleft_list)
				site_paths.append(sb_site_path) # add site path in as last entry to site path data structure
				site_positions.append(pos_for_path)

			elif state == 3: 
				# run forward algorithm to get logV for site-based: 
				logV_site, _ = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_sweep, self.site_states_sweep_list, self.classifierobj_sites)
				
				# run stochastic backtrace with logV_site to get path: 
				sb_site_path, pos_for_path = self.stochastic_backtrace(self.Svec_sites, logV_site, self.positions, position_indices, self.site_trans_mat_sweep, self.site_states_sweep_list)
				site_paths.append(sb_site_path) # add site path in as last entry to site path data structure
				site_positions.append(pos_for_path)

			elif state == 4: 
				# run forward algorithm to get logV for site-based: 
				logV_site, _ = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_linkedright, self.site_states_linkedright_list, self.classifierobj_sites)
				
				# run stochastic backtrace with logV_site to get path: 
				sb_site_path, pos_for_path = self.stochastic_backtrace(self.Svec_sites, logV_site, self.positions, position_indices, self.site_trans_mat_linkedright, self.site_states_linkedright_list)
				site_paths.append(sb_site_path) # add site path in as last entry to site path data structure
				site_positions.append(pos_for_path)
		
		# reverse site_paths list to get right order
		site_paths_final = site_paths[::-1]
		site_positions_final = site_positions[::-1]

		return window_state_path, site_paths_final, site_positions_final

	def many_hierarchical_backtraces(self, windows_out, sites_out, n_iter=100): 

		iterations_sites = []
		site_results = []
		site_pos = []

		window_results = []
		window_list = []
		iterations_windows = []

		for i in range(1, n_iter+1): 
			logV_window, _ = self.hierarchical_forward()
			window_path, site_paths, site_positions = self.hierarchical_stochastic_backtrace(logV_window)

			for j in range(1, len(window_path)-1): 
				window_results.append(window_path[j])
				window_list.append(self.windows[j-1])
				iterations_windows.append(i)

				site_path_without_zeros = site_paths[j-1][1:-1]
				site_results.append(site_path_without_zeros)
				site_pos.append(site_positions[j-1])
				iters = [i] * len(site_path_without_zeros)
				iterations_sites.append(iters)

		site_pos_flat = [item for sublist in site_pos for item in sublist]
		site_results_flat = [item for sublist in site_results for item in sublist]
		iterations_sites_flat = [item for sublist in iterations_sites for item in sublist]
		
		site_output = open(sites_out, 'w')
		site_output.write("site\tsite_pred\titer\n")
		for s in range(len(site_results_flat)): 
			site_output.write('\t'.join([str(site_pos_flat[s]), str(site_results_flat[s]), str(iterations_sites_flat[s]) + ("\n")]))
	
		window_output = open(windows_out, 'w')
		window_output.write("win\tpred\titer\n")
		for i in range(len(window_results)):
			window_output.write('\t'.join([str(window_list[i]), str(window_results[i]), str(iterations_windows[i]) + ("\n")]))

	def regular_window_hmm_sb(self, windows_out, n_iter): 

		Svec = self.Svec_windows
		position_indices = range(0, len(Svec))
		trans_mat = self.window_trans_mat
		state_list = self.window_state_list
		classifierobj = self.classifierobj_windows
		positions = self.windows

		window_results = []
		window_list = []
		iterations_windows = []

		for i in range(1, n_iter+1): 
			logV, _ = self.forward(Svec, position_indices, trans_mat, state_list, classifierobj)
			window_path, _ = self.stochastic_backtrace(Svec, logV, positions, position_indices, trans_mat, state_list)

			for j in range(1, len(window_path)-1): 
					window_results.append(window_path[j])
					window_list.append(self.windows[j-1])
					iterations_windows.append(i)
	
		window_output = open(windows_out, 'w')
		window_output.write("win\tpred\titer\n")
		for i in range(len(window_results)):
			window_output.write('\t'.join([str(window_list[i]), str(window_results[i]), str(iterations_windows[i]) + ("\n")]))
		
	def regular_site_hmm_sb(self, sites_out, n_iter): 

		Svec = self.Svec_sites
		position_indices = range(0, len(Svec))
		trans_mat = self.site_trans_mat_sweep
		state_list = self.site_states_sweep_list
		classifierobj = self.classifierobj_sites
		positions = self.positions

		site_results = []
		site_list = []
		iterations_sites = []

		for i in range(1, n_iter+1): 
			logV, _ = self.forward(Svec, position_indices, trans_mat, state_list, classifierobj)
			site_path, _ = self.stochastic_backtrace(Svec, logV, positions, position_indices, trans_mat, state_list)

			for j in range(1, len(site_path)-1): 
					site_results.append(site_path[j])
					site_list.append(self.positions[j-1])
					iterations_sites.append(i)
	
		site_output = open(sites_out, 'w')
		site_output.write("site\tsite_pred\titer\n")
		for i in range(len(site_results)):
			site_output.write('\t'.join([str(site_list[i]), str(site_results[i]), str(iterations_sites[i]) + ("\n")]))
	
	def call_viterbi(self, sites_out): # deprecated

		Svec = self.Svec_sites
		position_indices = np.where((np.array(self.positions) >= 480000) & (np.array(self.positions) <= 520000))
		position_indices = list(position_indices[0])
		trans_mat = self.site_trans_mat_sweep
		state_list = self.site_states_sweep_list
		classifierobj = self.classifierobj_sites
		positions = self.positions

		sites = positions[position_indices[0]:position_indices[-1]+1]
		
		_, state_path, _ = W.viterbi(Svec, position_indices, trans_mat, state_list, classifierobj)

		site_output = open(sites_out, 'w')
		site_output.write("site\tsite_pred\n")
		for s in range(1, len(state_path)): 
			site_output.write('\t'.join([str(sites[s-1]), str(state_path[s]) + ("\n")]))

	def hierarchical_posterior(self, windows_out, sites_out, include_probs = False): 

		window_forward, log_prob_fwd = self.hierarchical_forward()
		window_backward, log_prob_bwd = self.hierarchical_backward()

		n_window_states = len(self.window_state_list)
		n_window_intervals = len(self.windows)
		positions_arr = np.array(self.positions)
		SITE_LABEL = {
			(1, 1): 0, # neutral window: neutral site
			(2, 1): 1, # linkedleft window: linkedleft site
			(3, 1): 1, (3, 2): 2, (3, 3): 3, # sweep window: LL, sweep, LR sites
			(4, 1): 3, # linkedright window: linkedright site
		}

		window_results, window_probs_all = [], [] 
		site_results, site_post_all, site_marg_all, site_pos_all = [], [], [], []

		for t in range(1, n_window_intervals + 1): 

			# window posterior 
			log_post = ['n/a'] * n_window_states 
			for w in range(1, n_window_states): 
				if window_forward[t][w] != 'n/a' and window_backward[t][w] != 'n/a': 
					log_post[w] = window_forward[t][w] + window_backward[t][w]

			# check that there are defined posteriors in window
			defined = [x for x in log_post if x != 'n/a']
			if not defined: 
				raise RuntimeError(f"no window state has nonzero posterior at window {self.windows[t - 1]}")
			norm = self.logsumexp(defined)

			win_probs = [0.0] * n_window_states
			for w in range(1, n_window_states): 
				if log_post[w] != 'n/a': 
					win_probs[w] =  math.exp(log_post[w] - norm)

			best_w = int(np.argmax(win_probs))
			window_results.append(best_w)
			window_probs_all.append(win_probs)

			# site posterior within current window, conditional on best w
			win_start = self.windows[t - 1]
			win_end = win_start + self.window_size
			position_indices = list(np.where((positions_arr >= win_start) & (positions_arr < win_end))[0])
			n_sites = len(position_indices)

			site_probs = self._site_posterior(position_indices, best_w)

			for s in range(n_sites): 
				best_i = int(np.argmax(site_probs[s]))
				site_results.append(best_i)
				site_post_all.append(site_probs[s][best_i])

			if include_probs: 
				# needs every window state's site hmm, not just the chosen one
				site_probs_by_w = {best_w: site_probs}
				for w in range(1, n_window_states):
					if w not in site_probs_by_w:
						site_probs_by_w[w] = self._site_posterior(position_indices, w)
 
				for s in range(n_sites):
					marg = [0.0] * 4 # neutral, linkedleft, sweep, linkedright
					for w in range(1, n_window_states):
						for i in range(1, len(site_probs_by_w[w][s])):
							marg[SITE_LABEL[(w, i)]] += win_probs[w] * site_probs_by_w[w][s][i]
					site_marg_all.append(marg)
 
			site_pos_all.extend(positions_arr[position_indices].tolist())

		# write output
		print("Writing window and site posterior files...", flush = True)

		with open(sites_out, 'w') as site_output:
			header = "site\tsite_pred"
			if include_probs:
				header += "\tsite_post\tpost_neutral\tpost_linkedleft\tpost_sweep\tpost_linkedright"
			site_output.write(header + "\n")
			for s in range(len(site_results)):
				row = [str(site_pos_all[s]), str(site_results[s])]
				if include_probs:
					row += [f"{site_post_all[s]:.6g}"] + [f"{p:.6g}" for p in site_marg_all[s]]
				site_output.write('\t'.join(row) + "\n")
 
		with open(windows_out, 'w') as window_output:
			header = "win\tpred"
			if include_probs:
				header += "\tpost_neutral\tpost_linkedleft\tpost_sweep\tpost_linkedright"
			window_output.write(header + "\n")
			for i in range(len(window_results)):
				row = [str(self.windows[i]), str(window_results[i])]
				if include_probs:
					row += [f"{p:.6g}" for p in window_probs_all[i][1:]]
				window_output.write('\t'.join(row) + "\n")
 
		print("Done.\n")

	def _site_posterior(self, position_indices, w): 

		# site-level posterior for every site in window
		site_trans, site_states = self.get_site_config(w, n_sites = len(position_indices))		
		N = len(site_states)
		n_sites = len(position_indices)

		logF, _ = self.forward(self.Svec_sites, position_indices, site_trans, site_states, self.classifierobj_sites)
		logB, _ = self.backward(self.Svec_sites, position_indices, site_trans, site_states, self.classifierobj_sites)

		probs = []
		for s in range(1, n_sites + 1): 
			lp = ['n/a'] * N
			for i in range(1, N): 
				if logF[s][i] != 'n/a' and logB[s][i] != 'n/a': 
					lp[i] = logF[s][i] + logB[s][i]

			defined = [x for x in lp if x != 'n/a']
			p = [0.0] * N
			if defined: 
				norm = self.logsumexp(defined)
				for i in range(1, N): 
					if lp[i] != 'n/a': 
						p[i] =  math.exp(lp[i] - norm)
			probs.append(p)

		return probs

def main():

	def get_iteration_files(directory, ext_suffix):
		pat = re.compile(r'_(\d+)\.' + re.escape(ext_suffix) + r'$')
		all_files = glob.glob(os.path.join(directory, f"*.{ext_suffix}"), recursive = True)
		stem2file = {}
		for f in all_files:
			base = os.path.basename(f)
			m = pat.search(base) # search base name of file
			if not m:
				continue
			it = int(m.group(1))
			stem = base[:-len('.' + ext_suffix)]
			stem2file[stem] = f  
		
		return stem2file

	parser = argparse.ArgumentParser(description= "Arguments for WINDEX")

	parser.add_argument('--path2trained_sites', help = "Path to directory that contains trained site-level emissions")
	parser.add_argument('--path2trained_windows', help = "Path to directory that contains trained window-level emissions")
	parser.add_argument('--sites_out', help = "Desired output file name for site output")
	parser.add_argument('--windows_out', help = "Desired output file name for window output")
	parser.add_argument('--datafile_sites', help = "Path to site-based statistics file for single genome/simulation")
	parser.add_argument('--datafile_windows', help = "Path to window-based statistics file for single genome/simulation")
	parser.add_argument('--window_size', help = "window size used to calculate window-based statistics")
	parser.add_argument('--stochastic_backtrace', action = 'store_true', help = "run hierarchical stochastic backtrace in addition to hierarchical Viterbi evaluation")
	parser.add_argument('--baum_welch_training', action = 'store_true', help = "run Baum-Welch training on training files. Input both site and window files for this.")
	parser.add_argument('--baum_welch_matrix', help = "path to pretrained window transitions from Baum Welch (or other tuned transitions)")
	parser.add_argument('--posterior', action = 'store_true', help = "run hierarchical posterior decoding instead of hierarchical Viterbi")
	parser.add_argument('--posterior_probs', action = 'store_true', help = "append posterior probability columns to posterior output")
	args = parser.parse_args()

	W = WINDEX(args.path2trained_sites, args.path2trained_windows, args.datafile_sites, args.datafile_windows, int(args.window_size), args.baum_welch_training, args.baum_welch_matrix, alpha = 1.0, emissiontype = 'ODE')
	
	if args.baum_welch_training == True: # this is adapted for 1kgp files, generalize when forming into package

		site_dir = "stats/site/1kgp"
		win_dir = "stats/window/1kgp"

		site_stem2file = get_iteration_files(site_dir, "site_stats_std")
		window_stem2file = get_iteration_files(win_dir, "window_stats_std")

		common_stems = sorted(set(site_stem2file) & set(window_stem2file)) # make sure every site file has a corresponding window file 
		common_stems = [s for s in common_stems if s.endswith("_20")]

		#print(len(common_stems))

		training_datafiles_sites = [site_stem2file[i] for i in common_stems]
		training_datafiles_windows = [window_stem2file[i] for i in common_stems]

		W.load_training_data(training_datafiles_sites, training_datafiles_windows) 

		# run baum-welch training
		W.window_baum_welch(n_iter = 50, tol = 1e-6)
	
	# run hierarchical hmm
	if args.posterior:
		W.hierarchical_posterior(args.windows_out, args.sites_out, include_probs = args.posterior_probs)
	elif args.stochastic_backtrace: 
		W.many_hierarchical_backtraces(args.windows_out, args.sites_out, n_iter = 100)
	else:
		W.hierarchical_viterbi(args.windows_out, args.sites_out)
	
if __name__ == '__main__':
	main()

	
	
