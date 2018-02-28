#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from glob import glob
import re
import gzip


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

# import all variants into a dictionary with the ClinVar accession number (RCVXXXXXXXXX) as key
print "Loading %s..." % variant_summary_filename
all_clinvar_variants = {}
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
        rcv_accessions = variation_fields[11].split(";")
        assembly = variation_fields[16]
        chromosome_accession = variation_fields[17]
        chromosome = variation_fields[18]
        start = variation_fields[19]
        stop = variation_fields[20]
        if assembly == "GRCh38" and start.isdigit() and stop.isdigit() :
            for rcv in rcv_accessions:
                if rcv in all_clinvar_variants :
                    all_clinvar_variants[rcv].add((chromosome_accession, (int(start), int(stop)), gene_id, gene_symbol, allele_id, dbsnp_id, dbvar_id))   
                else :
                    all_clinvar_variants[rcv] = set([(chromosome_accession, (int(start), int(stop)), gene_id, gene_symbol, allele_id, dbsnp_id, dbvar_id)])   


for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")
    OMIM_folder = "%s/%s_OMIM" % (pathway_folder, pathway_name.replace(" ", "_"))
    CDS_folder = "%s/%s_CDS" % (pathway_folder, pathway_name.replace(" ", "_"))

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
                if (gene_id, gene_symbol) in cds_gene_map : 
                    cds_gene_map[(gene_id, gene_symbol)][protein_id] = (chromosome_accession, cds_intervals)
                else :
                    cds_gene_map[(gene_id, gene_symbol)] = {protein_id : (chromosome_accession, cds_intervals)}

    omim_variants = {}
    for allelic_variants in sorted(glob("%s/*_Allelic_Variants.txt" % OMIM_folder)):
        fields = os.path.basename(allelic_variants)[:-21].split("_")
        gene_id = fields[0]
        gene_symbol = fields[1]
        omim_id = fields[2]
        with open(allelic_variants, "r") as in_omim :
            for l in in_omim.readlines() :
                if l[0] != "." : continue
                f = l.replace(",", ";").rstrip("\n").split("\t")
                if len(f) < 5 : continue
                variant_number = f[0]
                dbsnp_id = f[3]
                match_found = False
                clinvar_ids = f[5].split(";;;")
                for rcv in clinvar_ids :
                    if rcv not in all_clinvar_variants :
                        print "%s not present in %s" % (rcv, variant_summary_filename)
                        continue
                    clinvar_variants = all_clinvar_variants[rcv]
                    for v in clinvar_variants : 
                        chromosome_accession = v[0]
                        interval1 = v[1]
                        allele_id = v[4]
                        if not dbsnp_id : dbsnp_id = v[5]
                        dbVar_id = v[6]
                        if (gene_id, gene_symbol) in cds_gene_map :
                            for protein_id in cds_gene_map[(gene_id, gene_symbol)] : 
                                if chromosome_accession == cds_gene_map[(gene_id, gene_symbol)][protein_id][0] :
                                    for interval2 in cds_gene_map[(gene_id, gene_symbol)][protein_id][1] :
                                        if getOverlap(interval1, interval2) :
                                            identifiers = "%s|%s|%s|%s|%s%s" % (allele_id, rcv, dbsnp_id, dbvar_id, omim_id, variant_number)
                                            if (gene_id, gene_symbol, protein_id) in omim_variants :
                                                omim_variants[(gene_id, gene_symbol, protein_id)].add(identifiers)
                                            else : 
                                                omim_variants[(gene_id, gene_symbol, protein_id)] = set([identifiers])
                                            with open("%s/%s_%s_%s_alleleIDs.txt" % (OMIM_folder, gene_id, gene_symbol, protein_id), "a") as out_ids : 
                                                out_ids.write("%s\n" % identifiers)
                                            break
                        else :
                            print "Unknown CDS for gene: %s:%s" %(gene_id, gene_symbol)
                            
    with open("%s/OMIM_gene_variants.txt" % pathway_folder, "w") as out_omim :
        for gene_id, gene_symbol, protein_id in omim_variants : 
            out_omim.write("%s:%s:%s:%s\n" % (gene_id, gene_symbol, protein_id, ",".join(sorted(omim_variants[(gene_id, gene_symbol, protein_id)]))))



