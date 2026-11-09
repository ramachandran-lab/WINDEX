"""
script to generate site and window-based statistics for WINDEX on GPU 
steps to manually add pops into haplotype matrix and parse windows for FST were assisted with Claude Sonnet v5
"""

import pandas as pd
import numpy as np 
from pg_gpu import HaplotypeMatrix, GenotypeMatrix, selection, divergence
from pg_gpu.windowed_analysis import windowed_statistics
import sys
import tskit as ts

def load_data(vcf_file): 

    # create haplotype matrix for haplotype stats
    h = HaplotypeMatrix.from_vcf(vcf_file)

    # load in populations (generic for all simulations) [ not working! ]
    # h.load_pop_file("scripts/simulations/pops")

    # change back to simulation pops when done 
    pop_df = pd.read_csv("pops.txt", sep = "\t", header = 0, names = ["sample", "pop"])
   
    # h.samples is the list of diploid sample names from the VCF
    # sample_sets maps pop_name to list of haplotype indices
    sample_index = {s: i for i, s in enumerate(h.samples)}

    sample_sets = {}
    for pop, group in pop_df.groupby("pop"):
        hap_indices = []
        for sample_id in group["sample"]:
            if sample_id in sample_index:
                i = sample_index[sample_id]
                hap_indices.extend([i, i + len(h.samples)])
        sample_sets[pop] = hap_indices

    h.sample_sets = sample_sets

    # create genotype matrix for diversity stats
    g = GenotypeMatrix.from_haplotype_matrix(h)

    return g, h

def compute_site_statistics(g, h, focal_pop, ref_pop, alt_pop, pol_file = None):
    
    # cross-population EHH (USED SELSCAN INSTEAD!)
    # xpehh_scores = selection.xpehh(h, focal_pop, ref_pop)

    # nSL
    nsl_scores = selection.nsl(h)

    # iHS
    ihs_scores = selection.ihs(h, include_edges = True)

    # FST (need to make it windowed to grab per-SNP FST)
    pos = h.get_positions().get()
    
    edges = np.concatenate([[pos[0] - 1], (pos[:-1] + pos[1:]) / 2, [pos[-1] + 1]])
    fst_scores = windowed_statistics(h, edges, statistics = ['fst'], pop1 = focal_pop, pop2 = ref_pop)['fst']

    if pol_file is not None:
        pol = (pd.read_csv(pol_file, sep = "\t", header = None, names = ["pos", "alt_derived"])
                 .drop_duplicates("pos", keep = "first"))
        m = pol.set_index("pos")["alt_derived"].reindex(np.asarray(pos).astype(np.int64)).to_numpy()
    else:
        m = np.ones(len(pos)) 

    def daf(p):
        f = np.asarray(h.haplotypes[h.sample_sets[p], :].mean(axis = 0).get())
        return np.where(m == 1, f, np.where(m == 0, 1.0 - f, np.nan))

    daf_focal, daf_ref, daf_alt = daf(focal_pop), daf(ref_pop), daf(alt_pop)

    ddaf_scores = pd.DataFrame({
        "DDAF": daf_focal - ((daf_ref + daf_alt) / 2)
    })

    return xpehh_scores, nsl_scores, ihs_scores, fst_scores, ddaf_scores, pos

def main(): 

    vcf_file = sys.argv[1]
    focal_pop = sys.argv[2]
    ref_pop = sys.argv[3]
    alt_pop = sys.argv[4]
    out_file = sys.argv[5]
    pol_file = sys.argv[6] if len(sys.argv) > 6 else None

    g, h = load_data(vcf_file)

    xpehh, nsl, ihs, fst, ddaf, pos = compute_site_statistics(g, h, focal_pop, ref_pop, alt_pop, pol_file)

    # combine to dataframe
    site_stats = pd.DataFrame({'SNP_name': np.array(pos), 'Physical_Distance': pos, 'Map_Distance': pos/1000000000, 'XP-EHH': xpehh, 'DDAF': np.squeeze(np.array(ddaf)), 'nSL': nsl, 'iHS': ihs, 'FST': fst})
    site_stats = site_stats.fillna(-998)
    site_stats.to_csv(out_file, sep = '\t', index = False)


if __name__ == '__main__':
	main()
