#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This script loads the list of genes from an input file and retrieves, for each gene, the official 
gene symbol according to the HUGO Gene Nomenclature Committee: HGNC database of human gene names.'''
import re
from Bio import Entrez
from time import sleep
Entrez.email = "sebastian.daberdaku@dei.unipd.it"     # Always tell NCBI who you are

for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")
    genes_file = "%s/%s_genes.txt" % (pathway_folder, pathway_name.replace(" ", "_"))
    # load gene ids and names from file
    genes = set()
    with open(genes_file, "r") as in_genes :
        entries = in_genes.read().split("\n///\n")
        for e in entries :
            m = re.search(r'DBLINKS     NCBI-GeneID: (.+)\n', e)
            if m : genes.add(int(m.group(1)))
    
    print "loaded %d genes" % len(genes)
    RetMax = 500000
    all_genes_text = ""
    for gene_id in sorted(genes) : 
        query = '%d[uid]' % gene_id
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="tabular", retmode="text")
        entry = handle.read()
        handle.close()
        if not all_genes_text : 
            all_genes_text += entry
            print entry
        else :
            all_genes_text += "\n" + entry.split("\n")[1]
            print entry.split("\n")[1]
        sleep(1)
    
    with open("%s/official_gene_symbols_and_names.txt" % pathway_folder , "w") as out_gene, open("%s/manually_added_genes.txt" % pathway_folder, "r") as manually_added_genes :
        header = manually_added_genes.readline()
        out_gene.write(all_genes_text + "\n" + manually_added_genes.read())
