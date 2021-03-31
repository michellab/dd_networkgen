#!/bin/python

import networkx as nx
import numpy as np
from cvxopt import matrix

import sys
sys.path.append("DiffNet")
import diffnet as dn
import graph as gph
from netbfe import *

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
matplotlib.rcParams['font.size'] = 20
matplotlib.rcParams['figure.dpi']= 100
import matplotlib.image as mpimg
from numpy.lib.stride_tricks import as_strided


import argparse
import csv


# argparse.
parser = argparse.ArgumentParser()
parser.add_argument("-matrix", type=str, help="NumPy array file to run DiffNet on.")
parser.add_argument("-max_sampling", type=int, help="Maximum degree of sampling (100-1000).", default=None)
parser.add_argument("-lig_img_path", type=str, help="Path to ligand images.")
parser.add_argument("-lig_order_path", type=str, help="Path to CSV file containing the order of ligands in s_\{ij\}.")
parser.add_argument("-ntwkpath", type=str, help="Path to write DiffNet graph to.")

args = parser.parse_args()

sij_path = args.matrix
max_sampling = args.max_sampling
imgs_path = args.lig_img_path
lig_order_path = args.lig_order_path
output_path = args.ntwkpath




def draw_diffnet_graph( g, lig_imgs, ligand_order, pos=None, ax=None, fig=None,
                        widthscale=None, nodescale=2.5, node_color=None, molscale=None, max_sampling=None,
                        origins=['O']):
    '''
    Draw a graph representing the difference network.
    
    Args:
    g: nx.Graph - the graph representing the difference network.
    pos: Kx2 numpy array or dict - the coordinates to place the nodes 
    in the graph. If numpy array, pos[i] is the coordinates for node i,
    excluding origin.  If dict, pos[i] is the coordinate of node i, including
    origin. If None, use a spring layout.
    Returns:
    pos: dict - pos[i] gives the positions of the node i.
    '''
    K = g.number_of_nodes() - len(origins)

    if isinstance( pos, np.ndarray):
        mypos = dict( [(i, pos[i]) for i in xrange(K)])
        if (len(pos) == K):
            for d, o in enumerate(origins):
                mypos.update( {o : (-1.0*d, -1.0*d)})
        else:
            for d, o in enumerate(origins):
                mypos.update( {o : pos[K+d]})
    elif type( pos) == dict:
        mypos = pos
    elif pos == "custom":
        #mypos = nx.nx_agraph.graphviz_layout(g, prog="neato")
        mypos = nx.circular_layout( g)
    else:
        mypos = nx.spring_layout( g)
    
    node_size = nodescale*K
    if node_color is None:
        node_color = 'red'
    nx.draw_networkx_nodes( g, mypos, nodelist=range(K), 
                            node_size=node_size,
                            node_color=node_color,
                            ax=ax, alpha=1,
                          )
#     nodeO = nx.draw_networkx_nodes( g, mypos, nodelist=origins,
#                                     node_size=node_size*2,
#                                     node_color='#FFFFFF',
#                                     width=2.,
#                                     ax=ax)
#     if node_color is None or len(node_color)<=K:
#         nodeO.set_edgecolor( 'red')
#     else:
#         nodeO.set_edgecolor( node_color[K:])

    if widthscale is None:
        widthscale = 5.*K

    weights = np.array( [ w for u, v, w in list(g.edges( data='weight')) ])
    weights[weights<0] = 0 # Set negative numbers to 0.
    

    if max_sampling:
        
        # figure out highest weights:
        sorted_weights = np.sort(weights)[::-1]
        sampling_checker = num_top_edges_to_keep = 0

        # keep adding highest weights until max_sampling threshold is reached
        for weight in sorted_weights:
            sampling_checker += weight
            
            # when threshold is reached, output the weight threshold
            if sampling_checker >= max_sampling:
                weight_threshold = weight
                break
                
        # replace every weight below threshold with nan to hide edges:
        weights[weights<weight_threshold] = np.nan 
    
    
    width = weights*widthscale
    nx.draw_networkx_edges( g, mypos, 
                            width=width,
                            ax=ax)
    
    
    # overlay molecular representations onto nodes:
    trans = ax.transData.transform
    trans2 = fig.transFigure.inverted().transform
    imsize=molscale
    for n, lig_name in zip(g.nodes(), ligand_order):

        img = lig_imgs[lig_name]
        (x,y) = mypos[n]
        xx,yy = trans((x,y)) # figure coordinates
        xa,ya = trans2((xx,yy)) # axes coordinates
        

        # add mol image:
        a = plt.axes([xa-imsize/2.0,ya-imsize/2.0, imsize, imsize ], zorder=100)
        a.imshow(img)
        

        a.axis('off')
    
    
    return mypos

if __name__ == "__main__":


	# load ligand images + the order they appear in sij:
	path_to_ligand_imgs = imgs_path
	path_to_ligand_idx = lig_order_path

	ligand_order = []
	with open(path_to_ligand_idx, "r") as f:
	    reader = csv.reader(f)
	    for row in reader:
	        ligand_order.append(row[0])
	        
	lig_images_dict = {}
	# for each ligand, read in the mol image and build dict:
	for ligand in ligand_order:

	    lig_img = mpimg.imread(path_to_ligand_imgs+ligand+".png")
	    lig_images_dict[ligand] = lig_img

	sij = np.load(sij_path)    

	sij = matrix(sij)

        
	N=1000
	integer_n=True

	n = networkBFEalloc( sij, N)
	if integer_n:
	    nint = round_to_integers( n)
	    n = matrix( nint[:], n.size, tc='d')
	    
	G = gph.diffnet_to_graph( n)
	fig, ax = plt.subplots( figsize=(10, 10))
	draw_diffnet_graph( G, lig_images_dict, ligand_order, 
	                   pos="custom", ax=ax, fig=fig, widthscale=200./N, 
	                   nodescale=700, molscale=0.1, max_sampling=max_sampling, node_color="white")
	ax.set_aspect(1)
	ax.axis('off')

	plt.savefig(output_path, dpi=300)

