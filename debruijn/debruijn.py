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
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")

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
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    isfile(fastq_file)
    with open(fastq_file, "r") as file_in:
        for line in file_in:
            yield next(file_in).strip()
            next(file_in)
            next(file_in)


def cut_kmer(read, kmer_size):
    for i in range(len(read)-(kmer_size-1)):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    read = "".join(list(read_fastq(fastq_file)))
    list_kmer = cut_kmer(read, kmer_size)
    kmer_dict = {}

    for kmer in list_kmer:
        if kmer in kmer_dict:
            kmer_dict[kmer]+=1
        else:
            kmer_dict[kmer] = 1

    return kmer_dict


def build_graph(kmer_dict):
    digraph = nx.DiGraph()
    for key, value in kmer_dict.items():
        digraph.add_edge(key[:-1], key[1:], weight=value)
    return digraph


def get_starting_nodes(graph):
    starting_nodes = []
    for node in list(graph.nodes()):
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return starting_nodes


def get_sink_nodes(graph):
    sinking_nodes = []
    for node in list(graph.nodes()):
        if len(list(graph.successors(node))) == 0:
            sinking_nodes.append(node)
    return sinking_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for start in starting_nodes:
        for end in ending_nodes:
            if nx.has_path(graph, start, end):
                for path in nx.all_simple_paths(graph, start, end):
                    seq = path[0]
                    for node in path[1:]:
                        seq+=node[-1]
                    contigs.append(tuple((seq, len(seq))))
    return contigs


def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as file_out:
        for i in range(len(contigs_list)):
            file_out.write(f">contig_{i} len={contigs_list[i][1]}\n")
            file_out.write(f"{fill(contigs_list[i][0])}\n")


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node and delete_sink_node == False:
            graph.remove_nodes_from(path[:-1])
        elif delete_entry_node == False and delete_sink_node:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph
    
        

def std(data):
    return statistics.stdev(data)
  

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    subgraph = graph.subgraph(path).edges(data=True)
    return statistics.mean([d["weight"] for (u, v, d) in subgraph])

def solve_bubble(graph, ancestor_node, descendant_node):
    if nx.has_path(graph, ancestor_node, descendant_node):
        length_list = []
        weight_avg_list = []
        path = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
        for subpath in path:
            length_list.append(len(subpath))
            weight_avg_list.append(path_average_weight(graph, subpath))
        if len(weight_avg_list) > 1:
            return select_best_path(graph, path, length_list, weight_avg_list)
    return graph


def simplify_bubbles(graph):
    bubble = False
    for node in graph.nodes:
        node_predecessor = list(graph.predecessors(node))
        nb_pred = len(node_predecessor)
        if nb_pred > 1:
            for i in range(nb_pred - 1):
                for j in range(i + 1, nb_pred):
                    node_ancestor = nx.lowest_common_ancestor(
                        graph, node_predecessor[i], node_predecessor[j])
                    if node_ancestor:
                        bubble = True
                        break
            if bubble:
                graph = simplify_bubbles(
                    solve_bubble(graph, node_ancestor, node))
                break
    return graph

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


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


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    start = get_starting_nodes(graph)
    end = get_sink_nodes(graph)
    contig_list = get_contigs(graph, start, end)
    save_contigs(contigs_list, args.output_file)
    
    #print(start)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
