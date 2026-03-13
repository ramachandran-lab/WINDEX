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

class WINDEX():
    
	def __init__(self, path2trained_sites, path2trained_windows, datafile_sites, datafile_windows, window_size, mode = "manual", relaxed = False, emissiontype = "NB"):
		self.name = None
		self.path2trained_sites = path2trained_sites
		self.path2trained_windows = path2trained_windows
		self.datafile_sites = datafile_sites
		self.datafile_windows = datafile_windows
		self.window_size = window_size
		self.emissiontype = emissiontype

		self.positions, self.Svec_sites, self.classifierobj_sites = self.read_site_statfiles()
		self.windows, self.Svec_windows, self.classifierobj_windows = self.read_window_statfiles()

		self.window_trans_mat, self.window_state_list, self.site_trans_mat_neutral, self.site_states_neutral_list, self.site_trans_mat_linkedleft, self.site_states_linkedleft_list, self.site_trans_mat_sweep, self. site_states_sweep_list, self.site_trans_mat_linkedright, self.site_states_linkedright_list = self.set_transitions(mode, relaxed)

	def read_site_statfiles(self):
		# open site info
		print("Reading site file...\n")
		classifierobj_sites = Classifiers(self.path2trained_sites)
		file_sites = open(self.path2trained_sites + 'component_stats.txt')
		stats_sites = file_sites.read()
		file_sites.close()
		statlist_sites = stats_sites.strip().splitlines()
		#statlist_sites = ['DDAF', 'iHS', 'nSL', 'FST'] ##### CHANGE THIS BACK ONCE DONE TESTING #########
		positions = []
		Svec_sites = []
		file_site_data = open(self.datafile_sites,'r')
		f = file_site_data.read()
		file_site_data.close()
		f = f.strip().splitlines()
		header = f[0]
		H = header.strip().split('\t')
		stat2index = {}
		for stat in statlist_sites:
			stat2index[stat] = H.index(stat)
		f = f[1:]
		for line in f:
			L = line.strip().split('\t')
			pos = int(L[H.index('SNP_name')])
			S = Stats(self.path2trained_sites)
			for stat in statlist_sites:
				S.stat2score[stat] = float(L[stat2index[stat]])
			#include = True
			nostats_list = []
			for i in range(len(statlist_sites)):
				nostats_list.append(S.stat2score[statlist_sites[i]])
				#if S.stat2score[statlist_sites[i]] == -998:
			#print(nostats_list)
			if any(num != -998 for num in nostats_list):
					include = True
			if include == True:
				positions.append(pos)
				Svec_sites.append(S)
			else: 
				print("removing site " + pos + "because it has no defined statistics")
		return positions, Svec_sites, classifierobj_sites
	
	def read_window_statfiles(self):
		# open window info
		print("Reading window file...\n")
		classifierobj_windows = Classifiers(self.path2trained_windows)
		file_windows = open(self.path2trained_windows + 'component_stats.txt')
		stats_windows = file_windows.read()
		file_windows.close()
		statlist_windows = stats_windows.strip().splitlines()
		windows = []
		Svec_windows = []
		file_window_data = open(self.datafile_windows,'r')
		f = file_window_data.read()
		file_window_data.close()
		f = f.strip().splitlines()
		header = f[0]
		H = header.strip().split('\t')
		stat2index = {}
		for stat in statlist_windows:
			stat2index[stat] = H.index(stat)
		f = f[1:]
		for line in f:
			L = line.strip().split('\t')
			w = float(L[H.index('SNP_name')])
			S = Stats(self.path2trained_windows)
			for stat in statlist_windows:
				S.stat2score[stat] = float(L[stat2index[stat]])
			include = True
			for i in range(len(statlist_windows)):
				if S.stat2score[statlist_windows[i]] == -998:
					include = True
			if include:
				windows.append(w)
				Svec_windows.append(S)

		return windows, Svec_windows, classifierobj_windows
	
	def set_transitions(self, mode = "manual", relaxed = False):

		if mode == "manual": # choose this if you would like to manually tune probability transitions for WINDEX. s
			# window matrix
			#if relaxed == False:
			# window_trans_mat = np.array([[0, 0.99, 0.01, 0, 0], [0.0001, 0.6999, 0.3, 0, 0], [0, 0, 0.6, 0.3, 0], [0, 0, 0, 0, 1], [0.00000000000001, 0.00000000000001, 0, 0, 0.9999999999998]]) # Neutral to S/E --> 0.00001; Neutral to Sweep --> 0.001
			# window_states = {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
			# window_state_list = list(window_states.values())

			window_trans_mat = np.array([[0, 0.999, 0.001, 0, 0], [0.0001, 0.9989, 0.001, 0, 0], [0, 0, 0.999, 0.001, 0], [0, 0, 0, 0, 1], [0.001, 0.001, 0, 0, 0.998]]) # Neutral to S/E --> 0.00001; Neutral to Sweep --> 0.001
			window_states = {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
			window_state_list = list(window_states.values())

			# window_trans_mat = np.array([[0, 1, 0, 0, 0], [0.00000001, 0.99999998, 0.00000001, 0, 0], [0, 0, 0.999, 0.001, 0], [0, 0, 0, 0, 1], [0.00000001, 0.00000001, 0, 0, 0.99999998]]) # Neutral to S/E --> 0.00001; Neutral to Sweep --> 0.001
			# window_states = {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
			# window_state_list = list(window_states.values())

			# window_trans_mat = np.array([[0, 0.99, 0.01, 0, 0], [0.0001, 0.6, 0.4, 0, 0], [0, 0, 0.99, 0.01, 0], [0, 0, 0, 0, 1], [0.001, 0.001, 0, 0, 0.998]]) # Neutral to S/E --> 0.00001; Neutral to Sweep --> 0.001
			# window_states = {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
			# window_state_list = list(window_states.values())

			# window_trans_mat = np.array([[0, 0.999, 0.001, 0, 0], [0.0000001, 0.9999998, 0.0000001, 0, 0], [0, 0, 0.999999, 0.000001, 0], [0, 0, 0, 0, 1], [0.000001, 0.00001, 0, 0, 0.99998]]) # Neutral to S/E --> 0.00001; Neutral to Sweep --> 0.001
			# window_states = {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
			# window_state_list = list(window_states.values())

			#else: 
				# window_trans_mat = np.array([[0, 0.99, 0.01, 0, 0], [1e-30, 0.6999999999, 0.3, 0, 0], [0, 0, 0.7, 0.3, 0], [0, 0, 0, 0.2, 0.8], [0, 0.00000000000000001, 0, 0, 0.9999999999999999]]) # Neutral to S/E --> 0.00001; Neutral to Sweep --> 0.001
				# window_states = {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
				# window_state_list = list(window_states.values())
					
			# site matrix for window state = neutral
			# site_trans_mat_neutral = [[0, 1], [0.001, 0.999]]
			# site_states_neutral = {'Start/End':'none', 'Neutral':'neutral'}
			# site_states_neutral_list = list(site_states_neutral.values())

			### new site matrix for window-state = neutral
			# site_trans_mat_neutral = [[0, 0.999, 0.001], [0.001, 0.998, 0.001], [0.001, 0, 0.999]]
			# site_states_neutral = {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked'}
			# site_states_neutral_list = list(site_states_neutral.values())

			site_trans_mat_neutral = [[0, 1], [0.001, 0.999]]
			site_states_neutral = {'Start/End':'none', 'Neutral':'neutral'}
			site_states_neutral_list = list(site_states_neutral.values())

			# site matrix for window state = linked left
			# site_trans_mat_linkedleft = [[0, 0.75, 0.25], [0, 0.9, 0.1], [0.01, 0, 0.99]]
			# site_states_linkedleft= {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked'}
			# site_states_linkedleft_list = list(site_states_linkedleft.values())

			### new site matrix for window-state = linked left
			site_trans_mat_linkedleft = [[0, 1], [0.001, 0.999]]
			site_states_linkedleft= {'Start/End':'none', 'LinkedLeft': 'linked'}
			site_states_linkedleft_list = list(site_states_linkedleft.values())

			# site matrix for window state = sweep
			site_trans_mat_sweep = [[0, 0.999, 0.001, 0], [0, 0.999999, 0.000001, 0], [0, 0, 0, 1], [0.01, 0, 0, 0.99]]
			site_states_sweep = {'Start/End':'none', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
			site_states_sweep_list = list(site_states_sweep.values())

			# # site matrix for window state = linked right
			# site_trans_mat_linkedright = [[0, 0, 1], [0.1, 0.9, 0], [0, 0.001, 0.999]]
			# site_states_linkedright= {'Start/End':'none', 'Neutral':'neutral', 'LinkedRight':'linked'}
			# site_states_linkedright_list = list(site_states_linkedright.values())

			# site matrix for window state = linked right
			site_trans_mat_linkedright = [[0, 1], [0.001, 0.999]]
			site_states_linkedright= {'Start/End':'none', 'LinkedRight':'linked'}
			site_states_linkedright_list = list(site_states_linkedright.values())

		elif mode == "parameterize": ### FINISH!!! ###
			param_A = 60000 # expected bp size of linked region (either in simulations or taken from real data)
			param_B = 1/1000000 # expected frequency of sweep site

			param_A_win = param_A / self.window_size # scale expected bp size of linked region by window size to get expected window size of linked region
			param_B_win = param_B * self.window_size # scale expected frequency of sweep site by window size to get expected frequency of linked windows

			# window matrix 
			#window_trans_mat = np.array([[0, (1 - (param_B_win / (1 - (param_A_win * param_B_win)))), (param_B_win / (1 - (param_A_win * param_B_win))), 0, 0], [0.0001, (1 - (param_B_win / (1 - (param_A_win * param_B_win)))), 0.9999*(param_B_win / (1 - (param_A_win * param_B_win))), 0, 0], [0, 0, 1 - (param_B_win / (1 - (param_A_win / 2))), (1 / (param_A_win / 2)), 0], [0, 0, 0, 0, 1], [0.0001, 0.9999*(1 / (param_A_win / 2)), 0, 0, (1 - (1 / (param_A_win / 2)))]]) 
			window_trans_mat = np.array([[0, 1, 0, 0, 0], [0.0001, (1 - (param_B_win / (1 - (param_A_win * param_B_win)))), 0.9999*(param_B_win / (1 - (param_A_win * param_B_win))), 0, 0], [0, 0, 1 - (param_B_win / (1 - (param_A_win / 2))), (1 / (param_A_win / 2)), 0], [0, 0, 0, 0, 1], [0.0001, 0.9999*(1 / (param_A_win / 2)), 0, 0, (1 - (1 / (param_A_win / 2)))]]) 
			window_states = {'Start/End':'none', 'Neutral':'neutral', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
			window_state_list = list(window_states.values())

			# site matrix for window state = neutral
			site_trans_mat_neutral = [[0, 1], [0.0001, 0.9999]]
			site_states_neutral = {'Start/End':'none', 'Neutral':'neutral'}
			site_states_neutral_list = list(site_states_neutral.values())

			# site matrix for window state = linked left
			site_trans_mat_linkedleft = [[0, 1], [0.0001, 0.9999]]
			site_states_linkedleft= {'Start/End':'none', 'LinkedLeft': 'linked'}
			site_states_linkedleft_list = list(site_states_linkedleft.values())

			# site matrix for window state = sweep
			site_trans_mat_sweep = [[0, 1, 0, 0], [0, 1 - (1/(param_A/2)), (1/(param_A/2)), 0], [0, 0, 0, 1], [0.0001, 0, 0, 0.9999]]

			site_states_sweep = {'Start/End':'none', 'LinkedLeft':'linked', 'Sweep':'sweep', 'LinkedRight':'linked'}
			site_states_sweep_list = list(site_states_sweep.values())

			# site matrix for window state = linked right
			site_trans_mat_linkedright = [[0, 1], [0.0001, 0.9999]]
			site_states_linkedright= {'Start/End':'none', 'LinkedRight':'linked'}
			site_states_linkedright_list = list(site_states_linkedright.values())

		return window_trans_mat, window_state_list, site_trans_mat_neutral, site_states_neutral_list, site_trans_mat_linkedleft, site_states_linkedleft_list, site_trans_mat_sweep, site_states_sweep_list, site_trans_mat_linkedright, site_states_linkedright_list

	def emission(self, scenario, S, classifierobj): 

		"""
		function to call emission probabilities from classifiers.py 
		"""
		
		if self.emissiontype == 'NB':
			Likelihood = classifierobj.nb_likelihood(scenario, S)
				
		else:
			Likelihood = classifierobj.ode_likelihood(scenario, S)
		
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
			a = max(sumvec)
			exponents = [x-a for x in sumvec]
			newsum = sum([math.exp(x) for x in exponents])
			newlogsum = math.log(newsum)
			return a+newlogsum

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
			win_end = self.windows[t] + self.window_size - 1

			position_indices = np.where((np.array(self.positions) >= win_start) & (np.array(self.positions) <= win_end))
			position_indices = list(position_indices[0])
		
			if len(position_indices) == 0: 
				print("WARNING: window " + str(win_start) + " does not contain any sites. Removing...")
				empty_windows.append(self.windows[t])

		indices_list = [index for index, item in enumerate(self.windows) if item in empty_windows]

		empty_windows = set(empty_windows)
		self.windows = [i for i in self.windows if i not in empty_windows]
		self.Svec_windows = [val for i, val in enumerate(self.Svec_windows) if i not in indices_list]

		#print(len(self.windows))
		#print(len(self.Svec_windows))
		
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

		# recursion over all windows
		print("Performing hierarchical Viterbi algorithm...\n")
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

					for prev_state in range(n_window_states):
						if window_viterbi[t-1][prev_state] != 'n/a' and self.window_trans_mat[prev_state][w] != 0:
							prob = (float(window_viterbi[t-1][prev_state]) + float(np.log(self.window_trans_mat[prev_state, w])) + float(site_log_probs))
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

					for prev_state in range(n_window_states):
						if window_viterbi[t-1][prev_state] != 'n/a' and self.window_trans_mat[prev_state][w] != 0:
							prob = (float(window_viterbi[t-1][prev_state]) + float(np.log(self.window_trans_mat[prev_state, w])) + float(site_log_probs))
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

					site_log_probs, site_path, site_pos = self.viterbi(self.Svec_sites, position_indices, self.site_trans_mat_sweep, self.site_states_sweep_list, self.classifierobj_sites)

					for prev_state in range(n_window_states):
						if window_viterbi[t-1][prev_state] != 'n/a' and self.window_trans_mat[prev_state][w] != 0:
							prob = (float(window_viterbi[t-1][prev_state]) + float(np.log(self.window_trans_mat[prev_state, w])) + float(site_log_probs))
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

					for prev_state in range(n_window_states):
						if window_viterbi[t-1][prev_state] != 'n/a' and self.window_trans_mat[prev_state][w] != 0:
							prob = (float(window_viterbi[t-1][prev_state]) + float(np.log(self.window_trans_mat[prev_state, w])) + float(site_log_probs))
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
		print("Writing site and window path files...\n")
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

		n_window_states = len(self.window_state_list)
		n_window_intervals = len(self.Svec_windows)

		window_forward = [['n/a' for _ in range(n_window_states)] for _ in range(n_window_intervals + 2)] 

		window_forward[0][0] = 0 # no emission

		# recursion 
		for t in range(1, n_window_intervals + 1):
			
			win_start = self.windows[t-1]
			win_end = win_start + self.window_size -1

			# get position indices 
			pos_window = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) < win_end)]
			position_indices = [self.positions.index(pos) for pos in pos_window]

			for w in range(1, n_window_states):

				if w == 1: 
					_, site_log_probsfinal = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_neutral, self.site_states_neutral_list, self.classifierobj_sites)
					sumvec = []
					for k in range(n_window_states):
						if window_forward[t-1][k] != 'n/a' and self.window_trans_mat[k][w] != 0:
							sumvec.append(float(window_forward[t-1][k]) + float(math.log(self.window_trans_mat[k][w])) + float(site_log_probsfinal))
			
					if len(sumvec) > 0:
						Sum = self.logsumexp(sumvec)
						window_forward[t][w] = Sum + math.log(self.emission(self.window_state_list[w], self.Svec_windows[t-1], self.classifierobj_windows))

				if w == 2: 
					_, site_log_probsfinal = self.forward(self.Svec_sites, position_indices,  self.site_trans_mat_linkedleft, self.site_states_linkedleft_list, self.classifierobj_sites)
					sumvec = []
					for k in range(n_window_states):
						if window_forward[t-1][k] != 'n/a' and self.window_trans_mat[k][w] != 0:
							sumvec.append(float(window_forward[t-1][k]) + float(math.log(self.window_trans_mat[k][w])) + float(site_log_probsfinal))
			
					if len(sumvec) > 0:
						Sum = self.logsumexp(sumvec)
						window_forward[t][w] = Sum + math.log(self.emission(self.window_state_list[w], self.Svec_windows[t-1], self.classifierobj_windows))

				if w == 3: 
					_, site_log_probsfinal = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_sweep, self.site_states_sweep_list, self.classifierobj_sites)
					sumvec = []
					for k in range(n_window_states):
						if window_forward[t-1][k] != 'n/a' and self.window_trans_mat[k][w] != 0:
							sumvec.append(float(window_forward[t-1][k]) + float(math.log(self.window_trans_mat[k][w])) + float(site_log_probsfinal))
			
					if len(sumvec) > 0:
						Sum = self.logsumexp(sumvec)
						window_forward[t][w] = Sum + math.log(self.emission(self.window_state_list[w], self.Svec_windows[t-1], self.classifierobj_windows))

				else: 
					_, site_log_probsfinal = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_linkedright, self.site_states_linkedright_list, self.classifierobj_sites)
					sumvec = []
					for k in range(n_window_states):
						if window_forward[t-1][k] != 'n/a' and self.window_trans_mat[k][w] != 0:
							sumvec.append(float(window_forward[t-1][k]) + float(math.log(self.window_trans_mat[k][w])) + float(site_log_probsfinal))
			
					if len(sumvec) > 0:
						Sum = self.logsumexp(sumvec)
						window_forward[t][w] = Sum + math.log(self.emission(self.window_state_list[w], self.Svec_windows[t-1], self.classifierobj_windows))

		Endsumvec = []
		for k in range(1, len(self.window_state_list)):
			if window_forward[n_window_intervals][k] != 'n/a' and self.window_trans_mat[k][0] != 0:
				Endsumvec.append(float(window_forward[n_window_intervals][k]) + float(math.log(self.window_trans_mat[k][0])))

		EndSum = self.logsumexp(Endsumvec)
		window_forward[n_window_intervals+1][0] = EndSum

		return window_forward, EndSum
		
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

		# setup of data structures
		n_window_intervals = len(self.Svec_windows)

		window_state_path = [0 for _ in range(n_window_intervals + 2)]
		unnormalized_logprobs = ['n/a' for _ in range(len(self.window_state_list))]
		position_indices = None
		site_paths = []
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
		
		state = np.random.choice(len(self.window_state_list),p=Probs)
		window_state_path[len(self.Svec_windows)] = state

		# get site paths of last initialized window state 
		win_start = self.windows[len(self.Svec_windows)-1]
		win_end = win_start + self.window_size -1

		pos_window = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) < win_end)]
		position_indices = [self.positions.index(pos) for pos in pos_window]

		if state == 1: 
			# run forward algorithm to get logV for site-based: 
			logV_site, _ = self.forward(self.Svec_sites, position_indices, self.site_trans_mat_neutral, self.site_states_neutral_list, self.classifierobj_sites)
			
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

		# move backwards through the rest of the states, getting site paths per randomly chosen window
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

			state = np.random.choice(len(self.window_state_list),p=Probs)
			window_state_path[i] = state
		
			# get position indices of last window
			win_start = self.windows[i]
			win_end = win_start + self.window_size -1

			# get position indices 
			position_indices = np.where((np.array(self.positions) >= win_start) & (np.array(self.positions) <= win_end))
			position_indices = list(position_indices[0])

			#pos_window = np.array(self.positions)[(np.array(self.positions) >= win_start) & (np.array(self.positions) < win_end)]
			#position_indices = [self.positions.index(pos) for pos in pos_window]

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
	
	def call_viterbi(self, sites_out):

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

	def main():

		parser = argparse.ArgumentParser(description= "Arguments for WINDEX")

		parser.add_argument('--path2trained_sites', help = "Path to directory that contains trained site-level emissions")
		parser.add_argument('--path2trained_windows', help = "Path to directory that contains trained window-level emissions")
		parser.add_argument('--sites_out', help = "Desired output file name for site output")
		parser.add_argument('--windows_out', help = "Desired output file name for window output")
		parser.add_argument('--datafile_sites', help = "Path to site-based statistics file for single genome/simulation")
		parser.add_argument('--datafile_windows', help = "Path to window-based statistics file for single genome/simulation")
		parser.add_argument('--window_size', help = "window size used to calculate window-based statistics")

		args = parser.parse_args()

		W = WINDEX(args.path2trained_sites, args.path2trained_windows, args.datafile_sites, args.datafile_windows, args.window_size)
		W.hierarchical_viterbi(args.windows_out, args.sites_out)
    	# W.many_hierarchical_backtraces(windows_out, sites_out, n_iter) for stochastic backtrace

	if __name__ == '__main__':
		main()
	
	
	
