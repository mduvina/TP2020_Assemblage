#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, "r") as fichier:
        for i, ligne in enumerate(fichier):
            if i % 4 == 1:
                ligne=ligne.strip()
                yield(ligne)

def cut_kmer(read, kmer_size):
    n = len(read)
    for i in range(0,n-kmer_size) :
        kmers=read[i:i+kmer_size]
        yield(kmers)


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict={}
    for ligne in read_fastq(fastq_file):
        for kmer in cut_kmer(ligne, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer]+=1
            else:
                kmer_dict[kmer]=1
    return kmer_dict


def build_graph(kmer_dict):
    graph=nx.DiGraph()
    for kmer in kmer_dict:
        prefix=kmer[:-1]
        suffix=kmer[1:]
        graph.add_edge(prefix, suffix, weight=kmer_dict[kmer])
    return graph

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node:
            graph.remove_nodes_from(path[:-1])
        if delete_sink_node:
            graph.remove_nodes_from(path[1:])
    return graph

def std(data):
    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    best1=[]
    best2=[]
    for i, path in enumerate(path_list):
        if path_average_weight(graph, path)==max(weight_avg_list):
            pass
        

def path_average_weight(graph, path):
    poids = []
    subgraph=graph.subgraph(path)
    for edge in subgraph.edges:
        poids.append(edge["weight"])
    average_weight = statistics.mean(poids)
    return average_weight

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    starting_nodes = []
    for node in graph.nodes():
        if any(graph.predecessors(node))==False:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph): 
    ending_nodes = []
    for node in graph.nodes():
        if any(graph.successors(node))==False:
            ending_nodes.append(node)
    return ending_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs_list = []
    for start in starting_nodes:
        for end in ending_nodes:
            for path in nx.all_simple_paths(graph, start, end):
                contig=[]
                for element in path:
                    contig.append(element)
                contigs_list.append(("".join(contig), len(contig)))
    return contigs_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as output:
        for i in range(len(contigs_list)):
            output.write(">contig_Numéro "+i+" len="+contigs_list[i][1]+"\n")
            output.write(fill(contigs_list[i][0])+"\n")
    output.close()

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    kmer_dict=build_kmer_dict(args.fastq_file, 5)
    graph=build_graph(kmer_dict)
    for node in get_starting_nodes(graph):
        print(node)
    for node in get_sink_nodes(graph):
        print(node)

    
if __name__ == '__main__':
    main()
