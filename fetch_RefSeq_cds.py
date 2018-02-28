#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This script retrieves the CDSs (Coding Region Sequences) for each RefSeqGene'''
from Bio import Entrez
from time import sleep
import re
import os

Entrez.email = "sebastian.daberdaku@dei.unipd.it"     # Always tell NCBI who you are

for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")
    output_folder = "%s/%s_CDS_RefSeq" % (pathway_folder, pathway_name.replace(" ", "_"))
    if not os.path.exists(output_folder): os.makedirs(output_folder)

    # load gene ids and names from file
    genes_by_chromosome = {}
    with open("%s/official_gene_symbols_and_names.txt" % pathway_folder, "r") as in_genes :
        header = in_genes.readline()
        for gene in in_genes.readlines() :
            gene_fields = gene.split("\t")
            gene_id = int(gene_fields[2])
            gene_symbol = gene_fields[5]
            chromosome = gene_fields[10]
            if chromosome in genes_by_chromosome :
                genes_by_chromosome[chromosome].add((gene_id, gene_symbol))
            else :
                genes_by_chromosome[chromosome] = {(gene_id, gene_symbol)}
    
    for chromosome in genes_by_chromosome :
        print chromosome
        for gene_id, gene_symbol in genes_by_chromosome[chromosome] :
            handle = Entrez.esearch(db="nucleotide", term="RefSeqGene[keyword] AND chromosome %s[Title]" % chromosome, retmax=100000)
            entries = Entrez.read(handle)
            handle.close()
            sleep(1)
            handle = Entrez.esummary(db="nucleotide", id=",".join(entries["IdList"]))
            summary = Entrez.read(handle)
            handle.close()
            sleep(1)
            for s in summary : 
                if gene_symbol in s["Title"] : 
                    accession = s["AccessionVersion"]
                    for f_format in ["na", "aa"] : 
                        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta_cds_%s" % f_format, retmode="text")
                        current_cds = handle.read()
                        handle.close()
                        cds_list = [">" + cds for cds in current_cds.split(">")[1:]]
                        for gene_cds in cds_list : 
                            m = re.search(r'\[gene=(.+)\] \[db_xref=(.*)GeneID:(\d+)\]', gene_cds)
                            if m : 
                                cds_gene_symbol = m.group(1)
                                cds_gene_id = int(m.group(3))
                                if gene_id == cds_gene_id and gene_symbol == cds_gene_symbol :
                                    print "%d\t%s" %(gene_id, gene_symbol)
                                    with open("%s/%d_%s_RefSeq_cds_%s.fasta" % (output_folder, gene_id, gene_symbol, f_format), "a") as out_cds :
                                        out_cds.write(gene_cds)
                        sleep(1)        