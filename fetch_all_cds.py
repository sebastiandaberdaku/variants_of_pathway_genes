#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This script retrieves the CDSs (Coding Region Sequences) for each gene.'''
from Bio import Entrez
from time import sleep
import re
import os

Entrez.email = "sebastian.daberdaku@dei.unipd.it"     # Always tell NCBI who you are

for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")
    output_folder = "%s/%s_CDS" % (pathway_folder, pathway_name.replace(" ", "_"))
    if not os.path.exists(output_folder): os.makedirs(output_folder)

    # load gene ids and names from file
    genes_by_chromosome = {}
    with open("%s/official_gene_symbols_and_names.txt" % pathway_folder, "r") as in_genes :
        header = in_genes.readline()
        for gene in in_genes.readlines() :
            gene_fields = gene.split("\t")
            gene_id = int(gene_fields[2])
            gene_symbol = gene_fields[5]
            chromosome = gene_fields[11]
            if chromosome in genes_by_chromosome :
                genes_by_chromosome[chromosome].add((gene_id, gene_symbol))
            else :
                genes_by_chromosome[chromosome] = {(gene_id, gene_symbol)}
       
    # For each chromosome fetch all of it's CDS regions 
    for chromosome in genes_by_chromosome :
        print chromosome
        for f_format in ["na", "aa"] : 
            handle = Entrez.efetch(db="nucleotide", id=chromosome, rettype="fasta_cds_%s" % f_format, retmode="text")
            chromosome_cds = handle.read()
            handle.close()
            # for each gene of interest in the current chromosome, extract the corresponding CDS
            all_cds = [">" + cds for cds in chromosome_cds.split(">")[1:]]
            for current_cds in all_cds : 
                m = re.search(r'\[gene=(.+)\] \[db_xref=(.*)GeneID:(\d+)\]', current_cds)
                if m : 
                    gene_symbol = m.group(1)
                    gene_id = int(m.group(3)) 
                    if (gene_id, gene_symbol) in genes_by_chromosome[chromosome] : 
                        print "%d\t%s" %(gene_id, gene_symbol)
                        with open("%s/%d_%s_cds_%s.fasta" % (output_folder, gene_id, gene_symbol, f_format), "a") as out_cds :
                            out_cds.write(current_cds)
            sleep(1)
