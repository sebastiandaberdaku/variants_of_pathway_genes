#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This script queries the OMIM database for known mutations affecting a given set of proteins/genes.'''

import os
import requests
from time import sleep

# The guys at Johns Hopkins University are not fond of this kind of crawling on their site and might ban the IP
# Better use Tor/Onion or a SSH tunnel ;)
PROXY_LIST = {
    "http": "socks5://127.0.0.1:9150",
    "https": "socks5://127.0.0.1:9150"
}

# load OMIM db
genemap2 = []
with open("./OMIM/genemap2.txt", "r") as in_omim :
    for entry in in_omim.readlines() :
        if entry[0] == "#" : continue
        genemap2.append(entry)

for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")
    output_folder = "%s/%s_OMIM" % (pathway_folder, pathway_name.replace(" ", "_"))
    if not os.path.exists(output_folder): os.makedirs(output_folder)
        
    # load gene id and names from file
    gene_ids = set()
    gene_symbol = {}
    with open("%s/official_gene_symbols_and_names.txt" % pathway_folder, "r") as in_genes :
        header = in_genes.readline()
        for gene in in_genes.readlines() :
            gene_fields = gene.split("\t")
            gene_id = gene_fields[2]
            gene_symbol[gene_id] = gene_fields[5]
            gene_ids.add(gene_id)
        
    with open("%s/OMIM_gene_entries.txt" % pathway_folder , "w") as out_omim :
        for entry in genemap2 :
            omim_id = entry.split("\t")[5]
            gene_id = entry.split("\t")[9]
            if gene_id in gene_ids :
                out_omim.write(entry)
                url = "https://www.omim.org/allelicVariant/%s?format=tsv" % omim_id
                print "Currently parsing MIM: %s, Gene ID: %s" % (omim_id, gene_id)
                print url
                try : 
                    req = requests.get(url, data=None, proxies=PROXY_LIST, headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/35.0.1916.47 Safari/537.36'})
                    response = req.text.encode('utf-8')
                    if not "OMIM Error" in response :
                        with open("%s/%s_%s_%s_Allelic_Variants.txt" %(output_folder, gene_id, gene_symbol[gene_id], omim_id), "w") as out_file :
                            out_file.write(response)
                except Exception,e:
                    print e
                    continue
                sleep(3) # not too aggressive...