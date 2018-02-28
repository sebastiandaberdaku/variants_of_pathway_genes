#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This script retrieves the latest version of the ClinVar database in VCF format. Variations are extracted for each gene of interest.'''

import vcf
import os
from glob import glob
import gzip

vcf_filename = "clinvar.vcf.gz"

# clinvar_vcf_link = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
# Download the latest ClinVar db in VCF format if "clinvar.vcf.gz" is not present
if not os.path.exists(vcf_filename): 
    from ftplib import FTP
    ncbi_server = "ftp.ncbi.nlm.nih.gov"
    clinvar_path = "pub/clinvar/vcf_GRCh38/"
    ftp = FTP(ncbi_server)
    ftp.login() 
    ftp.cwd(clinvar_path)
    with open(vcf_filename, 'wb') as out_vcf_file :
        ftp.retrbinary("RETR " + vcf_filename, out_vcf_file.write)
        ftp.quit()
'''
INFO keys:
AF_ESP:        Type=Float,Description="allele frequencies from GO-ESP"
AF_EXAC:       Type=Float,Description="allele frequencies from ExAC"
AF_TGP:        Type=Float,Description="allele frequencies from TGP"
ALLELEID:      Type=Integer,Description="the ClinVar Allele ID"
CLNDN:         Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB"
CLNDNINCL:     Type=String,Description="For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB"
CLNDISDB:      Type=String,Description="Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN"
CLNDISDBINCL:  Type=String,Description="For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN"
CLNHGVS:       Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression."
CLNREVSTAT:    Type=String,Description="ClinVar review status for the Variation ID"
CLNSIG:        Type=String,Description="Clinical significance for this single variant"
CLNSIGINCL:    Type=String,Description="Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance."
CLNVC:         Type=String,Description="Variant type"
CLNVCSO:       Type=String,Description="Sequence Ontology id for variant type"
CLNVI:         Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier"
DBVARID:       Type=String,Description="nsv accessions from dbVar for the variant"
GENEINFO:      Type=String,Description="Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)"
MC:            Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence"
ORIGIN:        Type=String,Description="Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other"
RS:            Type=String,Description="dbSNP ID (i.e. rs number)"
SSR:           Type=Integer,Description="Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other"
'''

'''
CLNVCSO - CLNVC

SO:0000667 - Insertion
SO:1000032 - Indel
SO:1000036 - Inversion
SO:0000159 - Deletion
SO:1000035 - Duplication
SO:0000289 - Microsatellite
SO:0001059 - Variation
SO:0001483 - single_nucleotide_variant
SO:0001742 - copy_number_gain
'''

coding_non_syn = set(['SO:0001587|nonsense', 'SO:0001589|frameshift_variant', 'SO:0001583|missense_variant'])

print "Loading all VCF entries..."
# loading clinvar.vcf
with open(vcf_filename, 'r') as in_vcf_file :
    vcf_reader = vcf.Reader(in_vcf_file)
    vcf_entries = {vcf_record.INFO["ALLELEID"] : vcf_record for vcf_record in vcf_reader}

for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")
    OMIM_folder = "%s/%s_ClinVar" % (pathway_folder, pathway_name.replace(" ", "_"))
    all_ClinVar_filename = "%s/ClinVar_gene_variants_in_VCF.txt" % pathway_folder
    if os.path.exists(all_ClinVar_filename): os.remove(all_ClinVar_filename)
    
    # load the AlleleIDs to extract from the VCF
    for alleleIDs_file in sorted(glob("%s/*_alleleIDs.txt" % OMIM_folder)):
        print "Currently parsing : " + os.path.basename(alleleIDs_file)
        with open(alleleIDs_file, "r") as in_allele_ids, open(all_ClinVar_filename, "a") as out_clinvar :
            f = os.path.basename(alleleIDs_file)[:-14]
            fields = f.split("_")
            gene_id = fields[0]
            gene_symbol = fields[1]
            protein_id = fields[2] + "_" + fields[3]
            output_records = []
            ids_list = []
            for l in in_allele_ids.readlines() :
                identifiers = l.rstrip("\n").split("|")
                allele_id = int(identifiers[0]) 
                rcv_accession = identifiers[1]
                dbsnp_id = identifiers[2]
                dbvar_id = identifiers[3]
                omim_id = identifiers[4]
                if allele_id in vcf_entries :
                    vcf_record = vcf_entries[allele_id]
                    if "MC" in vcf_record.INFO : 
                        if coding_non_syn.isdisjoint(vcf_record.INFO["MC"]) : 
                            continue
                        else :
                            print vcf_record.INFO["MC"]
                    output_records.append(vcf_record)
                    ids_list.append("%d|%s|%s|%s|%s" %(allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id))
                    
            if output_records : 
                with gzip.open(alleleIDs_file[:-14] + ".vcf.gz", "w") as vcf_out_file :
                    vcf_writer = vcf.Writer(vcf_out_file, vcf_reader)
                    for vcf_record in output_records : 
                        vcf_writer.write_record(vcf_record)
            if ids_list :
                out_clinvar.write("%s:%s:%s:%s\n" % (gene_id, gene_symbol, protein_id, ",".join(ids_list)))
