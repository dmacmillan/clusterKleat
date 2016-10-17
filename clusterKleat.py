from Kleat import *
import argparse, os, sys, time

# list must be sorted
# Iterative Agglomerative Heirarchical Clustering
def iterAHC(mylist, window=20):
    while True:
        mymin = float('inf')
        length = len(mylist)
        if length <= 1:
            return mylist
        for i in xrange(length-1):
            dist = mylist[i+1] - mylist[i]
            if dist < mymin:
                r = i
                s = i + 1
                mymin = dist
        if mymin <= window:
            mylist[r] = sum([mylist[r], mylist[s]])/2
            del(mylist[s])
        if mymin > window:
            return mylist

# Agglomerative Heirarchical Clustering
# Use iterative method instead to avoid recursion limit
#def AHC(mylist, window=20):
#    length = len(mylist)
#    if length <= 1:
#        return mylist
#    mymin = float('inf')
#    r = s = None
#    for i in xrange(length-1):
#        dist = mylist[i+1] - mylist[i]
#        if dist < mymin:
#            r = i
#            s = i + 1
#            mymin = dist
#    if mymin <= window:
#        mylist[r] = sum([mylist[r],mylist[s]])/2
#        del(mylist[s])
#        AHC(mylist, window)
#    return mylist

def sprint(text):
    sys.stdout.write(text)
    sys.stdout.flush()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Cluster multiple KLEAT outputs into a bed file')
    
    parser.add_argument('kleat', nargs='+', help='KLEAT output')
    parser.add_argument('-cw', '--cluster_window', type=int, default=15, help='Set the window size for clustering KLEAT cleavage sites. Default = 20')
    #parser.add_argument('-m', '--merged_kleat_output', action='store_true', help='Enable to output merged KLEAT file.')
    parser.add_argument('-s', '--stdout', action='store_true', help='Output to stdout instead of file')
    parser.add_argument('-n', '--name', default='clusters', help='Name for the file output. Default is "clusters"')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))
    
    args = parser.parse_args()

    if not args.stdout and not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    sites = []

    N = len(args.kleat)
    for i,kleat in enumerate(args.kleat):
        #sprint('Loading KLEAT {}/{}\r'.format(i,N))
        sites += Kleat.parseKleatFast(kleat)
    #print 'Loading KLEAT {}/{}\r'.format(N,N)

    #sprint('Grouping kleat data ...')
    sites = Kleat.groupKleat(sites)
    #print 'DONE'

    # Get the strandedness of each gene
    strands = {}

    kleats = {}

    # Clustering
    start = time.time()
    #sprint('Clustering cleavage events ...')
    for chrom in sites:
        kleats[chrom] = {}
        for gene in sites[chrom]:
            if not sites[chrom][gene]:
                continue
            # Get the strandedness of each gene
            strands[gene] = sites[chrom][gene][0].transcript_strand
            sites[chrom][gene] = sorted(set([x.cleavage_site for x in sites[chrom][gene]]))
            sites[chrom][gene] = iterAHC(sites[chrom][gene], args.cluster_window)

    handle = open(os.path.join(args.outdir, args.name), 'w') if not args.stdout else sys.stdout
    for chrom in sites:
        for gene in sites[chrom]:
            for site in sites[chrom][gene]:
                handle.write('{}\t{}\t{}\t{}\t0\t{}\n'.format(chrom, site - 1, site, gene, strands[gene]))
    handle.close()

    #print 'DONE ({}s)'.format(time.time()-start)
    #print 'SUCCESS, EXITING'
