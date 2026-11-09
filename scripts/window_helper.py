import read_posfile
import os
import sys
import allel
import pandas as pd
import numpy as np
import math

class SFS_statistics:
	def __init__(self,filename,start,end,numindividuals=100):
		self.filename = filename
		self.start = start
		self.end = end
		self.n = numindividuals
		self.etas = self.get_SFS()

	def get_SFS(self):
		'''
		:return: vector eta, where eta[i] is the number of segregating sites with DAF=i, i=0:self.n
		'''

		DAFs = read_posfile.dafs(self.filename,subset=True,start=self.start,end=self.end)
		DAFs = [int(round(x*self.n)) for x in DAFs]
		etas = [len([x for x in DAFs if x==eta]) for eta in range(self.n+1)]
		return etas

	def thetaW(self):
		'''
		:return: theta_W
		'''

		etas_nonfixed = [self.etas[i] for i in range(1,self.n)]
		S = sum(etas_nonfixed)
		a = float(1)/sum([float(1)/i for i in range(1,self.n)])
		return a*S

	def thetaPi(self):
		'''
		:return: theta_pi
		'''

		a = float(2)/(self.n*(self.n-1))
		S = 0
		for i in range(1,self.n):
			S += i*(self.n-i)*self.etas[i]
		return a*S

	def thetaH(self):
		'''
		:return: theta_H (Fay & Wu)
		'''

		a = float(2)/(self.n*(self.n-1))
		S = 0
		for i in range(1,self.n):
			S += i*i*self.etas[i]
		return a*S

	def thetaL(self):
		'''
		:return: theta_L (Zeng)
		'''
		a = float(1)/(self.n-1)
		S = 0
		for i in range(1,self.n):
			S += i*self.etas[i]
		return a*S

	def TajimaD(self):
		'''
		:return: Tajima's D, unnormalized (numerator only)
		'''
		return self.thetaPi()-self.thetaW()

	def FayWuH(self):
		'''
		:return: Fay & Wu's H, unnormalized (numerator only)
		'''
		return self.thetaPi()-self.thetaH()

	def ZengE(self):

		'''
		:return: Zeng's E, unnormalized (numerator only)
		'''
		return self.thetaL()-self.thetaW()
	
def scikitallel(pos, h1, ac1, ac2, ac3, size, start, end, step): 
	
	if step is None: 
		# windowed nucleotide diversity 
		pi, _, _, _ = allel.windowed_diversity(pos, ac1, size=size, start=start, stop=end)

		# windowed watterson theta 
		w, _, _, _ = allel.windowed_watterson_theta(pos, ac1, size=size, start=start, stop=end)

		# tajima's d
		d, _, _ = allel.windowed_tajima_d(pos, ac1, size=size, start=start, stop=end)

		return pi, w, d
	
	else: 
		# windowed nucleotide diversity 
		pi, _, _, _ = allel.windowed_diversity(pos, ac1, size=size, start=start, stop=end, step=step)

		# windowed watterson theta 
		w, _, _, _ = allel.windowed_watterson_theta(pos, ac1, size=size, start=start, stop=end, step=step)

		# tajima's d
		d, _, _ = allel.windowed_tajima_d(pos, ac1, size=size, start=start, stop=end, step=step)

		return pi, w, d
		

def allstatsline(posfilename, start, end, numindividuals = 200):

	'''
	:param posfilename: .pos file
	:param start: start position for window
	:param end: end position for window
	:param numindividuals: number of sequences
	:param normed: set to true to use values is means and stdevs to normalize
	:param means: dict of stat name to mean, learned using normalize_stats
	:param stdevs: dict of standard deviations, learned using normalize_stats
	:return: tab-delim line of <FILENAME,MIDDLEPOS,TW,TPi,TH,TL,TAJD,FWH,ZE> to put into allstats file
	'''
	S = SFS_statistics(posfilename,start,end,numindividuals=numindividuals)
	tW = S.thetaW()
	tPi = S.thetaPi()
	tH = S.thetaH()
	tL = S.thetaL()
	tajD = S.TajimaD()
	fwH = S.FayWuH()
	zE = S.ZengE()

	#middlepos = (start+end)/2
	#return '\t'.join([posfilename,str(middlepos),str(fwH),str(zE)])
	return fwH, zE

# from dave peede.
def fst_for_pbs(ac1, ac2):

	# compute
	num, den = allel.hudson_fst(ac1, ac2)

	# account for the denominator being zero
	fst = np.nansum(num) / np.nansum(den) if np.nansum(den) != 0 else 0

	# correct for negative Fst values
	return max(fst, 0)

# from dave peede.
def calc_pbs_per_region(g1, g2, g3): ### remove window stuff in here

	# determine if any sites are segregating in each population
	g1_has_seg = g1.count_alleles().is_segregating().any()
	g2_has_seg = g2.count_alleles().is_segregating().any()
	g3_has_seg = g3.count_alleles().is_segregating().any()

	# if any population has no segregating sites, PBS is undefined
	if not (g1_has_seg and g2_has_seg and g3_has_seg):
		return np.nan
	
	# else there are segregating sites
	else:

		# determine allele counts
		a_ac = g1.count_alleles()
		b_ac = g2.count_alleles()
		c_ac = g3.count_alleles()
		
		# calculate Fst.
		a_b_fst = fst_for_pbs(a_ac, b_ac)
		a_c_fst = fst_for_pbs(a_ac, c_ac)
		c_b_fst = fst_for_pbs(c_ac, b_ac)
		
		# correct for Fst values of 1 that will lead to inf
		a_b_fst = min(a_b_fst, 0.99999)
		a_c_fst = min(a_c_fst, 0.99999)
		c_b_fst = min(c_b_fst, 0.99999)
		
		# calculate PBS.
		pbs = (
			((np.log(1.0 - a_b_fst) * -1.0) +\
				(np.log(1.0 - a_c_fst) * -1.0) -\
				(np.log(1.0 - c_b_fst) * -1.0)) / 2.0
		)

	return pbs
	
def garuds_h(g1):
	
	# garud's h 
	h_1, h_12, h_123, h2_h1 = allel.garud_h(g1)

	return h_12

def segregating_sites(g1):

	seg_sites = g1.count_alleles().is_segregating()

	return sum(seg_sites)

if __name__ == '__main__':

	### setup 
	targetpop = allel.read_vcf(sys.argv[1])
	refpop = allel.read_vcf(sys.argv[2])
	altpop = allel.read_vcf(sys.argv[3])
	filename = sys.argv[4]
	output_path = sys.argv[5]

	gt1 = targetpop['calldata/GT']
	h1 = allel.HaplotypeArray(gt1.reshape(gt1.shape[0], -1))
	gt2 = refpop['calldata/GT']
	h2 = allel.HaplotypeArray(gt2.reshape(gt2.shape[0], -1))
	gt3 = altpop['calldata/GT']
	h3 = allel.HaplotypeArray(gt3.reshape(gt3.shape[0], -1))

	pos = targetpop['variants/POS'].astype(np.int64)
	pos_idx = allel.SortedIndex(pos)
	
	ac1 = h1.count_alleles()
	ac2 = h2.count_alleles()
	ac3 = h3.count_alleles()

	windowsize = int(sys.argv[6])
	start = 1
	end = int(sys.argv[7])

	numwindows = int(math.ceil(end / windowsize))

	start_ends = []
	for i in range(numwindows):
		win_start = start + i * windowsize
		win_end = win_start + windowsize - 1 
		start_ends.append([win_start, win_end])

	# # sliding window code
	# sliding = False  # true for sliding windows
	# step = 100 if sliding else windowsize
	# if sliding:
	# 	# only full windows
	# 	numwindows = int(math.floor((end - windowsize) / step)) + 1

	### get fwh, ze, pbs, garud's h, nss per win
	fwh_list = []
	ze_list = []
	garuds_list = []
	pbs_list = []
	nss_list = []

	#print(len(start_ends))
	for i in range(len(start_ends)):
		try: 
			positions = pos_idx.locate_range(start_ends[i][0], start_ends[i][1])
			
		except KeyError:
			# print("Error: window not found in data!")
			fwh_list.append(-998.0)
			ze_list.append(-998.0)
			garuds_list.append(-998.0)
			pbs_list.append(-998.0)
			nss_list.append(-998.0)
			
		else:
			g1 = allel.GenotypeArray(targetpop['calldata/GT'][positions])
			g2 = allel.GenotypeArray(refpop['calldata/GT'][positions])
			g3 = allel.GenotypeArray(altpop['calldata/GT'][positions])
			pbs = calc_pbs_per_region(g1, g2, g3)
			garuds = garuds_h(h1[positions])
			nss = segregating_sites(g1)
			fwh_out, ze_out = allstatsline(filename, start_ends[i][0],start_ends[i][1])
            
			fwh_list.append(fwh_out)
			ze_list.append(ze_out)
			garuds_list.append(garuds)
			pbs_list.append(pbs)
			nss_list.append(nss)
	
	fwh_final = np.array(fwh_list)
	ze_final = np.array(ze_list)
	garuds_final = np.array(garuds_list)
	pbs_final = np.array(pbs_list)
	nss_final = np.array(nss_list)

	### get pi, w, d, pbs, garud's h, seg sites per win
	pi, w, d = scikitallel(pos, h1, ac1, ac2, ac3, windowsize, start, end, step = None)

	def clean(x):
		return np.nan_to_num(np.asarray(x, dtype = float), nan = -998.0, posinf = -998.0, neginf = -998.0)

	pi, w, d = clean(pi), clean(w), clean(d)
	fwh_final = clean(fwh_final)
	ze_final = clean(ze_final)
	garuds_final = clean(garuds_final)
	pbs_final = clean(pbs_final)
	nss_final = clean(nss_final)

	#print("outputting file...")
	output = open(output_path, 'w')
	output.write("win_start\ttheta_pi\ttheta_w\ttajimas_d\tfay_wus_h\tzengs_e\tgaruds_h\tpbs\tnss\n")
	seq = range(start, end - windowsize + 2, windowsize)
	for i in range(len(seq)): 
		output.write('\t'.join([str(seq[i]),str(pi[i]),str(w[i]),str(d[i]),str(fwh_final[i]),str(ze_final[i]),str(garuds_final[i]),str(pbs_final[i]),str(nss_final[i])]) + ("\n"))
	output.close()



