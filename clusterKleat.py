from Kleat import *
import argparse
import os
import sys
import time
import numpy as np

from scipy.cluster.hierarchy import fclusterdata

def gen_cluster_index(x, window):
    return fclusterdata(x.reshape(-1, 1), 
                        window,
                        criterion='distance',
                        metric='euclidean',
                        method='centroid')

def cluster(x, window=15):
    cluster_idx = gen_cluster_index(x, window)
    dd = {}
    for i in np.unique(cluster_idx):
        vals = x[cluster_idx == i]
        mean = int(vals.mean())
        dd.update(zip(vals, np.repeat(mean, vals.shape[0])))
    return dd

def ahc(x, window=15):
    return np.unique(cluster(np.array(x), window).values()).tolist()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Cluster multiple KLEAT outputs into a bed file')
    
    parser.add_argument('kleat', nargs='+', help='KLEAT output')
    parser.add_argument('-cw', '--cluster_window', type=int, default=15, help='Set the window size for clustering KLEAT cleavage sites. Default = 20')
    parser.add_argument('-s', '--stdout', action='store_true', help='Output to stdout instead of file')
    parser.add_argument('-n', '--name', default='clusters', help='Name for the file output. Default is "clusters"')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))
    
    args = parser.parse_args()

    if not args.stdout and not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    sites = []

    N = len(args.kleat)
    for i,kleat in enumerate(args.kleat):
        sites += Kleat.parseKleatFast(kleat)

    sites = Kleat.groupKleat(sites)

    strands = {}

    kleats = {}

    handle = open(os.path.join(args.outdir, args.name), 'w') if not args.stdout else sys.stdout
    # Clustering
    start = time.time()
    for chrom in sites:
        kleats[chrom] = {}
        for gene in sites[chrom]:
            if not sites[chrom][gene]:
                continue
            # Get the strandedness of each gene
            strands[gene] = sites[chrom][gene][0].transcript_strand
            sites[chrom][gene] = np.unique([x.cleavage_site for x in sites[chrom][gene]])
            try:
                sites[chrom][gene] = ahc(sites[chrom][gene], args.cluster_window)
            except ValueError:
                sites[chrom][gene] = np.array(sites[chrom][gene])
            for site in sites[chrom][gene]:
                handle.write('{}\t{}\t{}\t{}\t0\t{}\n'.format(chrom, site - 1, site, gene, strands[gene]))
    #print 'Cluster time: {}'.format(time.time()-start)
    handle.close()
