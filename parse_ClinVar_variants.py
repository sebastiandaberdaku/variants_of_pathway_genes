#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This script retrieves the latest version of the ClinVar variant summary database. Variations are extracted for each gene of interest.'''

import gzip
import os
import re
from glob import glob


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)

variant_summary_filename = "variant_summary.txt.gz"

# summary_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
# Download the latest ClinVar db in tab delimited format if "variant_summary.txt.gz" is not present
if not os.path.exists(variant_summary_filename): 
    from ftplib import FTP
    ncbi_server = "ftp.ncbi.nlm.nih.gov"
    clinvar_path = "pub/clinvar/tab_delimited/"
    ftp = FTP(ncbi_server)
    ftp.login() 
    ftp.cwd(clinvar_path)
    with open(variant_summary_filename, 'wb') as out_file :
        print "Downloading %s " % variant_summary_filename
        ftp.retrbinary("RETR " + variant_summary_filename, out_file.write)
        ftp.quit()


# import all variants
print "Loading %s..." % variant_summary_filename
all_variants = {}
with gzip.open(variant_summary_filename, 'r') as f:
    header = f.readline()
    for entry in f.readlines() : 
        variation_fields = entry.replace(",", ";").split("\t")
        allele_id = variation_fields[0]
        gene_id = variation_fields[3]
        gene_symbol = variation_fields[4]
        if variation_fields[9].isdigit() : 
            dbsnp_id = "rs" + variation_fields[9]
        else : 
            dbsnp_id = ""
        if re.match(r"^(nsv|esv)\d+$", variation_fields[10]) : 
            dbvar_id = variation_fields[10]
        else : 
            dbvar_id = ""
        rcv_accession = variation_fields[11]
        assembly = variation_fields[16]
        chromosome_accession = variation_fields[17]
        start = variation_fields[19]
        stop = variation_fields[20]
        m = re.search(r"OMIM Allelic Variant:(\d+\.?\d+)", variation_fields[28])
        if m : omim_id = m.group(1)
        else : omim_id = ""
        identifiers = "%s|%s|%s|%s|%s" % (allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id)
        if assembly == "GRCh38" and start.isdigit() and stop.isdigit():
            if chromosome_accession in all_variants :
                all_variants[chromosome_accession].append((identifiers, (int(start), int(stop)), gene_id, gene_symbol))   
            else :
                all_variants[chromosome_accession] = [(identifiers, (int(start), int(stop)), gene_id, gene_symbol)]   
    
for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")
    CDS_folder = "%s/%s_CDS" % (pathway_folder, pathway_name.replace(" ", "_"))
    output_folder = "%s/%s_ClinVar" % (pathway_folder, pathway_name.replace(" ", "_"))
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    all_ClinVar_filename = "%s/ClinVar_gene_variants.txt" % pathway_folder

    # Load all CDSs of interest
    print "Loading all CDS regions of %s genes..." % pathway_name
    cds_gene_map = {}
    for cds_file in glob("%s/*_na.fasta" % CDS_folder):
        cd_filename = os.path.basename(cds_file)
        gene_id = cd_filename.split("_")[0]
        gene_symbol = cd_filename.split("_")[1]
        with open(cds_file, "r") as in_cds_file :
            for l in in_cds_file.readlines() :
                if l[0] != ">" : continue
                chromosome_accession = re.search(r">lcl\|(.+)\_cds", l).group(1)
                protein_id = re.search(r"\[protein\_id=(.+?)\]", l).group(1)
                cds_list = re.search(r"\[location=(?:complement\()?(?:join\()?(.+?)\){0,2}\]", l).group(1)
                cds_intervals = [(int(x.split("..")[0]), int(x.split("..")[1])) if ".." in x else (int(x), int(x)) for x in cds_list.split(",")]
                if chromosome_accession in cds_gene_map :
                    if (gene_id, gene_symbol) in cds_gene_map[chromosome_accession] : 
                        cds_gene_map[chromosome_accession][(gene_id, gene_symbol)][protein_id] = cds_intervals
                    else :
                        cds_gene_map[chromosome_accession][(gene_id, gene_symbol)] = {protein_id : cds_intervals}
                else : 
                    cds_gene_map[chromosome_accession] = {(gene_id, gene_symbol) : {protein_id : cds_intervals}}

    # Check all known variations in the ClinVar database if they affect the CDS regions of genes of interest.
    with open(all_ClinVar_filename, "w") as out_clinvar:
        for chromosome_accession in cds_gene_map :
            print "Evaluating variations in chromosome %s" % (chromosome_accession)
            for gene_id, gene_symbol in cds_gene_map[chromosome_accession] : 
                for protein_id in cds_gene_map[chromosome_accession][(gene_id, gene_symbol)] : 
                    # match with all variations from ClinVar in the current chromosome
                    ids_list = []
                    print "Evaluating %s:%s:%s" % (gene_id, gene_symbol, protein_id)
                    for variant in all_variants[chromosome_accession] : 
                        identifiers = variant[0]
                        current_interval = variant[1]
                        for cds_interval in cds_gene_map[chromosome_accession][(gene_id, gene_symbol)][protein_id] : 
                            if getOverlap(current_interval, cds_interval) :           
                                with open("%s/%s_%s_%s_alleleIDs.txt" % (output_folder, gene_id, gene_symbol, protein_id), "a") as out_file : 
                                    out_file.write("%s\n" % identifiers)
                                    ids_list.append(identifiers)
                                break
                    if ids_list :
                        out_clinvar.write("%s:%s:%s:%s\n" % (gene_id, gene_symbol, protein_id, ",".join(ids_list)))
