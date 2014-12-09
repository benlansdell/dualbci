#!/usr/bin/env ipython
import sys, os
import argparse
import requests
import json
import networkx as nx
import numpy as np
import scipy.io as scio
from xgmml_networkx import XGMMLWriter
from convnetx import *

def grangerNetwork(matfile, fn_out, threshold=1, pvalname='GCpvalB', weightname='GCdevB', clustername='clustersB', bciunits = 'brainunits'):
	#Load mat file
	a = scio.loadmat(matfile)
	weights = a[weightname]
	pvals = a[pvalname]
	names = [str(name[0]) for name in a['unitnames'][0]]
	bci = [str(name[0]) for name in a[bciunits][0]]
	usedinbci = [int(u in bci) for u in names]
	clusters = a[clustername][0]
	nodeclusters = {}
	for i in range(len(clusters)):
		for j in range(np.size(clusters[i][0])):
			node = clusters[i][0][j]
			nodeclusters[node] = i

	nU = np.size(weights,1)
	#Set pvals of zero to something very small
	weights[weights==0] = 1e-20
	#Set diagonal elements to 0
	for i in range(nU):
		weights[i,i]=0
	maxw = np.max(weights)
	#Create network
	g = nx.MultiDiGraph()
	#Color nodes according to sum of input and output change in deviance
	directionality = np.zeros((nU,1))
	sumcols = weights.sum(0)
	sumrows = weights.sum(1)
	for i in range(nU):
		directionality[i] = sumcols[i]-sumrows[i]
	#Add info about whether node was used in the BCI or not
	##
	round_to_n = lambda x, n: round(x, -int(np.floor(np.log10(x))) + (n - 1))
	[g.add_node(i, cluster = float(nodeclusters[i+1]), unitname = names[i], dirs=float(directionality[i]), bci=usedinbci[i]) for i in range(nU)]
	for i in range(nU):
		for j in range(nU):
			if weights[i,j]>0:
				g.add_edge(i,j, weight=round_to_n(float(weights[i,j]),4), pval=round_to_n(float(pvals[i,j]),4))
	#Write to graphml file
	nx.write_graphml(g, fn_out, prettyprint = False);
	#fn_out = fn_out + '.gml'
	#nx.write_gml(g, fn_out)
	#to_sif(g, fn_out)
	#to_graphml(g, fn_out)
	#to_gexf(g, fn_out)

def runGranger():
	matfile = './GLMGrangerB.mat'; fn_out = './GLMGrangerB.xml';
	grangerNetwork(matfile, fn_out, threshold=1, pvalname='GCpvalB', weightname='GCdevB', clustername='clustersB', bciunits = 'brainunits')
	matfile = './GLMGrangerBP.mat'; fn_out = './GLMGrangerBP.xml';
	grangerNetwork(matfile, fn_out, threshold=1, pvalname='GCpvalBP', weightname='GCdevBP', clustername='clustersBP', bciunits = 'brainunits')
	matfile = './GLMGrangerM.mat'; fn_out = './GLMGrangerM.xml';
	grangerNetwork(matfile, fn_out, threshold=1, pvalname='GCpvalM', weightname='GCdevM', clustername='clustersM', bciunits = 'manualunits')
	matfile = './GLMGrangerMP.mat'; fn_out = './GLMGrangerMP.xml';
	grangerNetwork(matfile, fn_out, threshold=1, pvalname='GCpvalMP', weightname='GCdevMP', clustername='clustersMP', bciunits = 'manualunits')

def main(argv):
    usage = """networkx_granger.py <matfile> <xmlfile> --threshold [val] --matname [name]

	Create networkx object from clustered Granger causality code

	Requires a .mat file containing a square matrix (by default named GCpval) which contains
	the p-values testing for Granger causality. The entry in row i and column j
	indicates unit i's "G-causal effect" on unit j.	    

	Writes to xml file for importing into cytoscape.
   
    Ben Lansdell: ben.lansdell@gmail.com
    2014
    """
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('matfile', type=str,
        help='matfile containing Granger causality data')
    parser.add_argument('fn_out', type=str,
        help='filename for output XML file')
    parser.add_argument('--threshold', type=float, default=1.0,
        help='threshold p-value above which to ignore connection')
    parser.add_argument('--matname', type=str, default='GCpval',
        help='name of matrix within mat file to plot')
    args = parser.parse_args()
    grangerNetwork(args.matfile, args.fn_out, args.threshold, args.matname)

if __name__ == "__main__":
    sys.exit(main(sys.argv))