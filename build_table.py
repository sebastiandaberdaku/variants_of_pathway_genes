#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Load genes in Rule-based insulin signaling model
insulin_signaling_RBA_genes = set()
with open("./insulin_signaling/genes_in_rule_based_model.txt", "r") as in_model : 
    for l in in_model.readlines() : 
        if l[0] == "#" : continue
        genes = l.rstrip().split("\t")[1].split(",")
        for g in genes : 
            gene_id = g.split(":")[0]
            gene_symbol = g.split(":")[1]
            insulin_signaling_RBA_genes.add((gene_id, gene_symbol))
# Load genes in full FBA tryptophan metabolism model
tryptophan_metabolism_FBA_genes = set()
with open("./tryptophan_metabolism/genes_in_FBA_model.txt", "r") as in_model : 
    for l in in_model.readlines() : 
        if l[0] == "#" : continue
        genes = l.rstrip().split("\t")[1].split(",")
        for g in genes : 
            gene_id = g.split(":")[0]
            gene_symbol = g.split(":")[1]
            tryptophan_metabolism_FBA_genes.add((gene_id, gene_symbol))
# Load genes in full reduced FBA tryptophan metabolism model
tryptophan_metabolism_reduced_FBA_genes = set()
with open("./tryptophan_metabolism/genes_in_reduced_FBA_model.txt", "r") as in_model : 
    for l in in_model.readlines() : 
        if l[0] == "#" : continue
        genes = l.rstrip().split("\t")[1].split(",")
        for g in genes : 
            gene_id = g.split(":")[0]
            gene_symbol = g.split(":")[1]
            tryptophan_metabolism_reduced_FBA_genes.add((gene_id, gene_symbol))
            
                    
# Load variants from ClinVar
ClinVar_variants = {}
ClinVar_variants_by_dbsnp_id = {}

for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")

    with open("%s/ClinVar_gene_variants_in_VCF.txt" % pathway_folder, "r") as in_ClinVar : 
        for l in in_ClinVar.readlines() :
            fields = l.rstrip().split(":")
            gene_id = fields[0]
            gene_symbol = fields[1]
            protein_id = fields[2]
            variation_ids = [tuple(identifiers.split("|")) for identifiers in fields[3].split(",")]
            model = ""
            if pathway_name == "insulin signaling" : 
                if (gene_id, gene_symbol) in insulin_signaling_RBA_genes : 
                    model = "Rule-based insulin signaling pathway model"
            else :
                if (gene_id, gene_symbol) in tryptophan_metabolism_FBA_genes : 
                    model = "FBA tryptophan metabolism pathway model"
                if (gene_id, gene_symbol) in tryptophan_metabolism_reduced_FBA_genes : 
                    if model : 
                        model += "|reduced FBA tryptophan metabolism pathway model"
                    else : 
                        model = "reduced FBA tryptophan metabolism pathway model"

            for allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id in variation_ids :
                if allele_id in ClinVar_variants :
                    ClinVar_variants[allele_id].append((gene_id, gene_symbol, protein_id, pathway_name + " (from KEGG)", model, (allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id)))
                else :
                    ClinVar_variants[allele_id] = [(gene_id, gene_symbol, protein_id, pathway_name + " (from KEGG)", model, (allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id))]
                if dbsnp_id : 
                    if dbsnp_id in ClinVar_variants_by_dbsnp_id : 
                        ClinVar_variants_by_dbsnp_id[dbsnp_id].append((gene_id, gene_symbol, protein_id, pathway_name + " (from KEGG)", model, (allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id)))
                    else : 
                        ClinVar_variants_by_dbsnp_id[dbsnp_id] = [(gene_id, gene_symbol, protein_id, pathway_name + " (from KEGG)", model, (allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id))]

                        
            # Load variants from OMIM
OMIM_variants = {}
for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")

    with open("%s/OMIM_gene_variants_in_VCF.txt" % pathway_folder, "r") as in_OMIM : 
        for l in in_OMIM.readlines() :
            fields = l.rstrip().split(":")
            gene_id = fields[0]
            gene_symbol = fields[1]
            protein_id = fields[2]
            variation_ids = [tuple(identifiers.split("|")) for identifiers in fields[3].split(",")]

            model = ""
            if pathway_name == "insulin signaling" : 
                if (gene_id, gene_symbol) in insulin_signaling_RBA_genes : 
                    model = "Rule-based insulin signaling pathway model"
            else :
                if (gene_id, gene_symbol) in tryptophan_metabolism_FBA_genes : 
                    model = "FBA tryptophan metabolism pathway model"
                if (gene_id, gene_symbol) in tryptophan_metabolism_reduced_FBA_genes : 
                    if model : 
                        model += "|reduced FBA tryptophan metabolism pathway model"
                    else : 
                        model = "reduced FBA tryptophan metabolism pathway model"

            for allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id in variation_ids :
                if omim_id in OMIM_variants :
                    OMIM_variants[omim_id].append((gene_id, gene_symbol, protein_id, pathway_name + " (from KEGG)", model, (allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id)))
                else :
                    OMIM_variants[omim_id] = [(gene_id, gene_symbol, protein_id, pathway_name + " (from KEGG)", model, (allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id))]
     
# Load variants from dbSNP
dbSNP_variants = {}
for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")

    with open("%s/dbSNP_gene_variants.txt" % pathway_folder, "r") as in_dbSNP : 
        for l in in_dbSNP.readlines() :
            fields = l.rstrip().split(":")
            gene_id = fields[0]
            gene_symbol = fields[1]
            protein_id = fields[2]
            snp_list = fields[3].split(",")
            
            model = ""
            if pathway_name == "insulin signaling" : 
                if (gene_id, gene_symbol) in insulin_signaling_RBA_genes : 
                    model = "Rule-based insulin signaling pathway model"
            else :
                if (gene_id, gene_symbol) in tryptophan_metabolism_FBA_genes : 
                    model = "FBA tryptophan metabolism pathway model"
                if (gene_id, gene_symbol) in tryptophan_metabolism_reduced_FBA_genes : 
                    if model : 
                        model += "|reduced FBA tryptophan metabolism pathway model"
                    else : 
                        model = "reduced FBA tryptophan metabolism pathway model"

            for snp_id in snp_list :
                allele_ids = rcv_accessions = dbsnp_ids = dbvar_ids = omim_ids = ""
                if snp_id in ClinVar_variants_by_dbsnp_id : 
                    #allele_id, rcv_accession, dbsnp_id, dbvar_id, omim_id
                    allele_ids = ";".join(sorted(set([e[5][0] for e in ClinVar_variants_by_dbsnp_id[snp_id]])))
                    rcv_accessions = ";".join(sorted(set([e[5][1] for e in ClinVar_variants_by_dbsnp_id[snp_id]])))
                    dbsnp_ids = ";".join(sorted(set([e[5][2] for e in ClinVar_variants_by_dbsnp_id[snp_id]])))
                    dbvar_ids = ";".join(sorted(set([e[5][3] for e in ClinVar_variants_by_dbsnp_id[snp_id]])))
                    omim_ids = ";".join(sorted(set([e[5][4] for e in ClinVar_variants_by_dbsnp_id[snp_id]])))

                if snp_id in dbSNP_variants :
                    dbSNP_variants[snp_id].add((gene_id, gene_symbol, protein_id, pathway_name + " (from KEGG)", model, (allele_ids, rcv_accessions, dbsnp_ids, dbvar_ids, omim_ids)))
                else :
                    dbSNP_variants[snp_id] = set([(gene_id, gene_symbol, protein_id, pathway_name + " (from KEGG)", model, (allele_ids, rcv_accessions, dbsnp_ids, dbvar_ids, omim_ids))])
                    
                    

    
with open("variants_table.txt", "w") as output_table : 
    header = "#Variation ID\tDatabase\tGene ID\tGene symbol\tProteinID\tpathway\t[AlleleID]|[ClinVar accession]|[Reference SNP ID]|[Variant Region Accession]|[Allelic Variant]\tmodel\n"
    output_table.write(header)
    for dbsnp_id in sorted(dbSNP_variants) : 
        for entry in dbSNP_variants[dbsnp_id] : 
            gene_id, gene_symbol, protein_id, pathway_name, model, (allele_ids, rcv_accessions, dbsnp_ids, dbvar_ids, omim_ids) = entry
            other_ids = "%s|%s|%s|%s|%s" % (allele_ids, rcv_accessions, "", dbvar_ids, omim_ids)
            output_table.write("%s\tdbSNP\t%s\t%s\t%s\t%s\t%s\t%s\n" % (dbsnp_id, gene_id, gene_symbol, protein_id, pathway_name, other_ids, model))
    for allele_id in sorted(ClinVar_variants) :
        for entry in ClinVar_variants[allele_id] : 
            gene_id, gene_symbol, protein_id, pathway_name, model, (allele_ids, rcv_accessions, dbsnp_ids, dbvar_ids, omim_ids) = entry
            other_ids = "%s|%s|%s|%s|%s" % ("", rcv_accessions, dbsnp_ids, dbvar_ids, omim_ids)
            output_table.write("%s\tClinVar\t%s\t%s\t%s\t%s\t%s\t%s\n" % (allele_id, gene_id, gene_symbol, protein_id, pathway_name, other_ids, model))
    for omim_id in sorted(OMIM_variants) :
        for entry in OMIM_variants[omim_id] :
            gene_id, gene_symbol, protein_id, pathway_name, model, (allele_ids, rcv_accessions, dbsnp_ids, dbvar_ids, omim_ids) = entry
            other_ids = "%s|%s|%s|%s|%s" % (allele_id, rcv_accessions, dbsnp_ids, dbvar_ids, "")
            output_table.write("%s\tOMIM\t%s\t%s\t%s\t%s\t%s\t%s\n" % (omim_id, gene_id, gene_symbol, protein_id, pathway_name, other_ids, model))
