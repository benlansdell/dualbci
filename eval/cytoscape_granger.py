#!/usr/bin/env ipython
import sys, os
import argparse
import requests
import json
import networkx as nx
from py2cytoscape import util as cy 
from collections import OrderedDict
import numpy as np
import scipy.io as scio

class CytoScape:
	"""Class to interface with Cytoscape via cyREST JS API"""
	def __init__(self):
		self.PORT_NUMBER = 1234
		self.BASE = 'http://localhost:' + str(self.PORT_NUMBER) + '/v1/'
		self.HEADERS = {'Content-Type': 'application/json'}
		self.suid = None
		try:
			self.clearNetworks()
			print 'Successfully connected to Cytoscape.'
		except requests.ConnectionError:
			print 'ConnectionError: Make sure Cytoscape is running and cyREST is installed'

	def clearNetworks(self):
		"""Delete all networks in current session"""
		return requests.delete(self.BASE + 'networks')

	def addNetwork(self, g):
		"""Convert a NetworkX Graph object and add it to cytoscape""" 
		cytoscape_network = cy.from_networkx(g)
		res1 = requests.post(self.BASE + 'networks', data=json.dumps(cytoscape_network), headers=self.HEADERS)
		res1_dict = json.loads(res1.content)
		self.suid = res1_dict['networkSUID']

	def layoutForceDirected(self):
		"""Change layout to force directed"""
		requests.get(self.BASE + 'apply/layouts/force-directed/' + str(self.suid))

	def layoutGroupAttribute(self, attr):
		requests.get(self.BASE + 'apply/layouts/attributes-layout/' + str(self.suid))
		#requests.get(cs.BASE + 'apply/layouts/attributes-layout/' + str(cs.suid))

	def saveView(self, fn_out):
		"""Save to eps file"""
		imgurl = self.BASE+'networks/' + str(self.suid) + '/views/first.png'
		self.__system('wget ' + imgurl + ' > ' + fn_out)

	def __system(self, cmd):
		import subprocess
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out = p.stdout.read()
		err = p.stderr.read()
		return (out, err)

def testNetwork(cs):
	"""Test network"""
	g = nx.MultiDiGraph()
	g.add_node(1)
	g.add_node(2)
	g.add_node(3)
	g.add_edge(1, 2)
	g.add_edge(2, 1)
	g.add_edge(1, 3)
	g.add_edge(3, 2)
	cs.addNetwork(g)

def grangerNetwork(matfile, cs, fn_out, threshold=1, pvalname='GCpvalB', weightname='GCdevB', clustername='clustersB'):
	#Load mat file
	a = scio.loadmat(matfile)
	weights = a[weightname]
	pvals = a[pvalname]
	names = [str(name[0]) for name in a['unitnames'][0]]
	clusters = a[clustername][0]
	nodeclusters = {}
	for i in range(len(clusters)):
		for j in range(np.size(clusters[i][0])):
			node = clusters[i][0][j]
			nodeclusters[node] = i

	nU = np.size(weights,1)
	#print names
	#print nU
	#Set pvals of zero to something very small
	weights[weights==0] = 1e-20
	#Set pvals above threshold to 1
	#weights[weights>threshold] = 1
	#Set diagonal elements to 1
	for i in range(nU):
		weights[i,i]=0
	maxw = np.max(weights)
	#Create network with weights as -log(p-value)
	logpvals = -np.log(weights)
	#Create network
	cs.clearNetworks()
	g = nx.MultiDiGraph()
	[g.add_node(i, cluster = nodeclusters[i+1], unitname = names[i]) for i in range(nU)]
	for i in range(nU):
		for j in range(nU):
			if weights[i,j]>0:
				g.add_edge(i,j, weight=weights[i,j], pval=pvals[i,j])
	cs.addNetwork(g)
	#Change view to group by cluster
	#cs.layoutForceDirected()
	#cs.layoutGroupAttribute('cluster')

	#Change display to universe, directed

	#Save to file
	#cs.saveView(fn_out)

def runGranger():
	matfile = './GLMGrangerB.mat'; fn_out = './GLMGrangerB.png';
	grangerNetwork(matfile, cs, fn_out, threshold=1, pvalname='GCpvalB', weightname='GCdevB', clustername='clustersB')
	matfile = './GLMGrangerBP.mat'; fn_out = './GLMGrangerBP.png';
	grangerNetwork(matfile, cs, fn_out, threshold=1, pvalname='GCpvalBP', weightname='GCdevBP', clustername='clustersBP')
	matfile = './GLMGrangerM.mat'; fn_out = './GLMGrangerM.png';
	grangerNetwork(matfile, cs, fn_out, threshold=1, pvalname='GCpvalM', weightname='GCdevM', clustername='clustersM')
	matfile = './GLMGrangerMP.mat'; fn_out = './GLMGrangerMP.png';
	grangerNetwork(matfile, cs, fn_out, threshold=1, pvalname='GCpvalMP', weightname='GCdevMP', clustername='clustersMP')

def main(argv):
    usage = """cytoscape_granger.py <matfile> <epsfile> --threshold [val] --matname [name]

	Draw Granger Causality connectivity graph using force-directed visualization:
	Connections with stronger connectivity are clustered together.

	Requires a .mat file containing a square matrix (by default named GCpval) which contains
	the p-values testing for Granger causality. The entry in row i and column j
	indicates unit i's "G-causal effect" on unit j.	    

	Writes image to .eps file.
   
    Ben Lansdell: ben.lansdell@gmail.com
    2014
    """
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('matfile', type=str,
        help='matfile containing Granger causality data')
    parser.add_argument('fn_out', type=str,
        help='filename for output EPS file')
    parser.add_argument('--threshold', type=float, default=1.0,
        help='threshold p-value above which to ignore connection')
    parser.add_argument('--matname', type=str, default='GCpval',
        help='name of matrix within mat file to plot')
    args = parser.parse_args()
    cs = CytoScape()
    grangerNetwork(args.matfile, cs, args.fn_out, args.threshold, args.matname)

if __name__ == "__main__":
    sys.exit(main(sys.argv))