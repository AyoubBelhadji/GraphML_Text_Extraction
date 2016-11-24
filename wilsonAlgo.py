# -*- coding: utf-8 -*-
"""

@author : Xi SHEN, Ayoub Belhadji
@compiler: python 2.7, netwokx

Generate Uniform Spanning Tree, Algo de Wilson

Ref : Generating random spanning trees more quickly than the cover time
"""

import numpy as np
import matplotlib.pyplot as pl
import networkx as nx
from random import random, choice, sample ## generate random value

''' certain examples of using nx
---------------------### Example 1 ###-------------------------------

G = nx.Graph()  #Creat a graph                                         
G.add_node(1)   #add a node 1                                   
G.add_edge(2,3) #add nodes : 2 and 3, at the same time add an DIRECTED edge from 2 to 3                                   
G.add_edge(3,2) #add nodes : 3 and 2, at the same time add an DIRECTED edge from 2 to 3
		#than the edge between 2 and 3 become UNDIRECTED                 
print G.nodes()                                      
print G.edges()                                      
print G.number_of_edges()  

nx.draw(G)  # plot
pl.show()

---------------------### Example 2 ###-------------------------------
G=nx.cubical_graph()
pos=nx.spring_layout(G) # positions for all nodes

# nodes
nx.draw_networkx_nodes(G,pos,
                       nodelist=[0,1,2,3],
                       node_color='r',
                       node_size=500,
                   alpha=0.8)
nx.draw_networkx_nodes(G,pos,
                       nodelist=[4,5,6,7],
                       node_color='b',
                       node_size=500,
                   alpha=0.8)

# edges
nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
nx.draw_networkx_edges(G,pos,
                       edgelist=[(0,1),(1,2),(2,3),(3,0)],
                       width=8,alpha=0.5,edge_color='r')
nx.draw_networkx_edges(G,pos,
                       edgelist=[(4,5),(5,6),(6,7),(7,4)],
                       width=8,alpha=0.5,edge_color='b')


# some math labels
labels={}
labels[0]=r'$a$'
labels[1]=r'$b$'
labels[2]=r'$c$'
labels[3]=r'$d$'
labels[4]=r'$\alpha$'
labels[5]=r'$\beta$'
labels[6]=r'$\gamma$'
labels[7]=r'$\delta$'
nx.draw_networkx_labels(G,pos,labels,font_size=16)

plt.axis('off')
plt.savefig("labels_and_colors.png") # save as png
plt.show() # display
'''

## Let's write the algorithm
RG = nx.random_graphs.random_regular_graph(3,20) 
G = RG
G=nx.cubical_graph()
def Uniform_Spanning_Tree ( G ) :
	if not nx.is_connected(G) :
		raise RuntimeError('Current graph is not connected...')
	if G.is_directed() : 
		raise RuntimeError('Current graph is directed...')
	visited_nodes = set(sample(G.nodes(), 1)) ## for the initialization, let's randomly peak two nodes
	unvisited_nodes = set(G.nodes()) - visited_nodes
	edges_path = [] ## store the edges 

	while len(unvisited_nodes) != 0  :
		## randomly peak a unvisited nodes as start
		node_start = sample(unvisited_nodes, 1)[0] 
		nodes_pass = np.array( [ node_start ] )
		next_node = node_start
		while next_node not in visited_nodes : 
			next_node = sample(G[next_node], 1)[0]
			## remove the loop
			if next_node in nodes_pass :
				index_node = np.where(nodes_pass == next_node)[0][0]
				nodes_pass = nodes_pass[: index_node + 1]
			else :
				nodes_pass = np.append(nodes_pass, next_node)
		## Update visited nodes and unvisited nodes
		visited_nodes |= set(nodes_pass)
		unvisited_nodes =unvisited_nodes - set(nodes_pass)
		for i in range(len(nodes_pass) - 1):
			edges_path.append((nodes_pass[i], nodes_pass[i+1]))
	spanning_tree = nx.Graph()
	spanning_tree.add_edges_from(edges_path)
	return edges_path, spanning_tree

edges_path, spanning_tree = Uniform_Spanning_Tree ( G )




## plot
pos=nx.spring_layout(G) # positions for all nodes
pl.figure(1)
ax1 = pl.subplot(121)
nx.draw_networkx_nodes(G, pos, node_color='r', node_size=500, alpha=0.8)

nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)

nx.draw_networkx_edges(G, pos, edgelist=edges_path, width=8, alpha=0.5, edge_color='r')

labels={}
for i in range(len(G.nodes())) :
	labels[i]=i
nx.draw_networkx_labels(G, pos, labels, font_size=16 )
ax1.set_title('Original graphs and spanning tree generated (red edges)')

ax2 = pl.subplot(122)
nx.draw(spanning_tree, pos, node_color='r', node_size=500) 
ax2.set_title('The spanning tree generated from the left figure')

pl.show()


##########-------------- Let's verify certain properties -----------------###########


## First property 
## Counting the spanning : Nb of spanning tree = det(Laplacian remove ith row and ith column)
## Take a simple example to avoid the computational complexity 
G = nx.complete_graph(4)
print 'Nb of spanning tree is : %d'%round(np.linalg.det(nx.laplacian_matrix(G)[1 :, 1: ]))

nbSpanningtree = 16
hist_spanning_tree = [] # dictionary to store the spanning tree and to count the statistics 
hist_spanning_tree_Laplacian = []
couter = np.zeros(nbSpanningtree)
nbIteration = 10000
for i in range(nbIteration) : 
	_, spanning_tree = Uniform_Spanning_Tree ( G )
	laplacian = nx.laplacian_matrix(spanning_tree)
	find_tree = False
	if hist_spanning_tree : 
		for i in range(len(hist_spanning_tree)) : 
			if np.linalg.norm(laplacian - hist_spanning_tree_Laplacian[i]) == 0 :
				couter[i] += 1
				find_tree = True
	if not find_tree : 
		hist_spanning_tree.append(spanning_tree)
		hist_spanning_tree_Laplacian.append(laplacian)
		try :
			couter[len(hist_spanning_tree) - 1] += 1
		except : 
			raise RuntimeError('Number of spanning tree is not correct.')

if len(hist_spanning_tree) == nbSpanningtree : 
	print 'Graphs get %d spanning tree, the propoerty of MTT is verified' %nbSpanningtree
else: 
	raise RuntimeError('Graphs get %d spanning tree, which should be %d' %(len(hist_spanning_tree), nbSpanningtree))

pl.figure (2) 
for i in range(4):
	for j in range(4):
		tree = hist_spanning_tree[i * 4 + j]
		ax = pl.subplot(4, 6, i * 6 + j +1)
		pos=nx.spring_layout(tree)
		nx.draw(tree, pos, node_color='r', node_size=500)
		ax.set_title('%d'%(i * 4 + j + 1))

ax = pl.subplot(3, 3, 6)
pos=nx.spring_layout(G)
nx.draw(G, pos, node_color='r', node_size=500)
ax.set_title('Full graph')
pl.show()

pl.bar(range(nbSpanningtree), couter, 0.35, color='r')
pl.title('The frequency of %d generating spanning tree'%nbIteration)
pl.xlabel('Index of spanning tree')
pl.ylabel('Nb of appearance')
pl.show()







	
