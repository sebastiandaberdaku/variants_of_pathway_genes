#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This script retrieves all known genes/proteins involved in the human 
'insulin signaling' and 'tryptophan metabolism' pathways.'''

from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.REST import kegg_get, kegg_list
from time import sleep
import os

for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    outdir_pathway = "./%s" % pathway_name.replace(" ", "_")
    if not os.path.exists(outdir_pathway): os.makedirs(outdir_pathway)
    # find human pathways in the KEGG database
    # /list/pathway/hsa    returns the list of human pathways
    entries = kegg_list("pathway", "hsa").read().split("\n")
    # get the 'pathway' entry identifier
    for e in entries :
        if pathway_name in e.lower() :
            pathway_entry = e.split()[0]
            print e
            break
    # save the pathway to file for future use
    with open('%s/%s.xml' % (outdir_pathway, pathway_entry),'wb') as out_kgml:
        out_kgml.write(kegg_get(pathway_entry, "kgml").read())
    # load the pathway from file
    with open('%s/%s.xml' % (outdir_pathway, pathway_entry), 'r') as in_kgml :
        pathway = KGML_parser.read(in_kgml)
      
    print pathway
     
    # extract all genes involved in this pathway, we use a set because the same gene could participate in multiple reactions
    all_genes = set() 
    for g in pathway.genes :
        all_genes |= set(g.name.split())
    print "total number of retrieved genes: %d" % len(all_genes)
     
    with open("%s/%s_genes.txt" %(outdir_pathway, pathway_name.replace(" ", "_")), "w") as out_genes :
        for g in all_genes : 
            print g
            out_genes.write(kegg_get(g).read())
            sleep(1) # wait between KEGG get requests 
