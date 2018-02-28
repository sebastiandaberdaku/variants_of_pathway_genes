#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This script queries the NCBI dbSNP database for coding non-synomous SNPs for a given set of proteins/genes.
The database can be searched programmatically using the Entrez module from the BioPython package.

The "Function Class" search field accepts the following indexed search terms list, 
an asterix (*) preceeds the coding non-synomous sequence variants:
*cds indel
    indel snp with length of multiple of 3bp, not causing frameshift
downstream variant 500b : replaced dbSNP term: near-gene-3
near gene 3 : replaced by downstream variant 500b
    A sequence variant located within a half KB of the end of a gene.
upstream variant 2kb : replaced dbSNP term: near-gene-3
near gene 5 : replaced by upstream variant 2kb
    A sequence variant located within 2KB 5' of a gene.
*frameshift : replaced by frameshift variant
    An attribute describing a sequence that contains a mutation involving the deletion or insertion of one or more bases, where this number is not divisible by 3.
*frameshift variant : replaced dbSNP term frameshift
    A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three.
gene segment : replaced by nc transcript variant
nc transcript variant : replaced dbSNP term gene-segment
    A transcript variant of a non-coding RNA gene.
intron : replaced by intron variant
    A region of a primary transcript that is transcribed, but removed from within the transcript by splicing together the sequences (exons) on either side of it.
intron variant : eplaced dbSNP term: intron
    A transcript variant occurring within an intron.
*missense
    A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved.
*nonsense : replaced by "stop gained" and "stop lost"
    A point mutation in a sequence of DNA that results in a premature stop codon, or a point-nonsense codon in the transcribed mRNA, and in a truncated, incomplete, and usually nonfunctional protein product
reference
    The allele is the same as the contig (contig reference) and hence causes no change to the translated sequence.
splice 3
splice acceptor variant
    A splice variant that changes the 2 base region at the 3' end of an intron.
splice 5
splice donor variant
    A splice variant that changes the 2 base pair region at the 5' end of an intron.
*stop gain
*stop gained
    A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened polypeptide.
*stop loss
*stop lost
    A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript.
synonymous codon
    A sequence variant where there is no resulting change to the encoded amino acid.
utr 3
utr 5
    Messenger RNA sequences that are untranslated and lie five prime or three prime to sequences which are translated.
utr variant 3 prime
    A UTR variant of the 3' UTR.
utr variant 5 prime
    A UTR variant of the 5' UTR.
'''

import os
from glob import glob
from Bio import Entrez
from time import sleep
import re

Entrez.email = "sebastian.daberdaku@dei.unipd.it"     # Always tell NCBI who you are

all_cls = set(['synonymous-codon', 'splice-acceptor-variant', 'reference', 'stop-gained', 'stop-lost', 'nc-transcript-variant', 'utr-variant-3-prime', 'utr-variant-5-prime', 'missense', 'downstream-variant-500B', 'intron-variant', 'splice-donor-variant', 'upstream-variant-2KB', 'cds-indel', 'frameshift-variant'])
coding_non_syn = set(["non-synonymous-codon", "coding-sequence-variant", "nonsense", "missense", "frameshift-variant", "cds-indel", "stop-gained", "stop-lost", "complex-change-in-transcript", "incomplete-terminal-codon-variant"])

def snps_in_gene(gene_id, gene_symbol, protein_id, flat_file_report):
    '''Returns the list of coding non-synomous SNPs for the current gene.
    
    Keyword arguments:
    gene_id: NCBI gene id
    gene_symbol: official name for the current gene
    protein_id: NCBI protein id
    flat_file_report: string containing the query results for the current gene
    '''
    # list of coding non-synomous SNPs for the current gene: initially empty
    rsIDs = set()
    # each snp entry ends with double newline "\n\n"
    for snp in flat_file_report.split("\n\n") : 
        rsID = re.search(r"\A(rs\d+) | ", snp).group(1)
        molecular_consequences = re.findall(r"LOC \| %s \| locus_id=%s \| fxn-class=(.+?) \| .*?prot_acc=%s\n" %(gene_symbol, gene_id, protein_id.replace(".", "\.")), snp)
        if not coding_non_syn.isdisjoint(molecular_consequences) : 
            rsIDs.add(rsID)
    return sorted(rsIDs)
    
'''Coding - Non-Synonymous: cds indel, frameshift variant, missense, stop gained, stop lost.'''
coding_non_syn_keywords = ["cds indel", "frameshift variant", "missense", "stop gained", "stop lost"]
fnx_class = " OR ".join('"%s"[Function Class]' %x for x in coding_non_syn_keywords)

for pathway_name in ["insulin signaling", "tryptophan metabolism"] :
    pathway_folder = "./%s" % pathway_name.replace(" ", "_")
    output_folder = "%s/%s_dbSNP" % (pathway_folder, pathway_name.replace(" ", "_"))
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    all_SNPs_filename = "%s/dbSNP_gene_variants.txt" % pathway_folder
    if os.path.exists(all_SNPs_filename): os.remove(all_SNPs_filename)
    
    # load genes and proteins
    proteins_per_gene = {}
    for cds_file in sorted(glob("%s/%s_CDS/*_cds_na.fasta" % (pathway_folder, pathway_name.replace(" ", "_")))):
        cd_filename = os.path.basename(cds_file)
        gene_id = cd_filename.split("_")[0]
        gene_symbol = cd_filename.split("_")[1]
        with open(cds_file, "r") as in_cds_file :
            proteins = []
            for l in in_cds_file.readlines() :
                if l[0] != ">" : continue
                protein_id = re.search(r"\[protein\_id=(.+?)\]", l).group(1)
                proteins.append(protein_id)
            proteins_per_gene[(gene_id, gene_symbol)] = proteins
    
    print "loaded %d genes" % len(proteins_per_gene)

    RetMax = 500000
    for gene_id, gene_symbol in sorted(proteins_per_gene) : 
        print "Protein ID: %s" % gene_id 
        # Here we build the query to retrieve coding non-synomous SNPs only for the current gene
        query = '"homo sapiens"[Organism] AND %s[Gene Name] AND (%s)' % (gene_symbol, fnx_class)
        print query
        handle = Entrez.esearch(db="snp", term=query, retmax=RetMax)
        record = Entrez.read(handle)
        handle.close()
        assert (int(record["Count"]) <= RetMax)
        sleep(1)
        print "Retrieved: %s records." % record["Count"]
        assert (len(record["IdList"]) == int(record["Count"]))
        with open(all_SNPs_filename, "a") as out_snp :
            handle = Entrez.efetch(db="snp", id=",".join(record["IdList"]), mode = "text", report='flt')
            # A SNP can belong to more than one gene. The Entrez query might accept SNPs which  
            # are coding non-synomous for other genes but not for the one(s) we are interested in.
            content = handle.read().rstrip()
            handle.close()
            for protein_id in proteins_per_gene[(gene_id, gene_symbol)] :
                print "Evaluating: " + gene_id + ":" + gene_symbol + ":" + protein_id
                snp_ids = snps_in_gene(gene_id, gene_symbol, protein_id, content)
                if snp_ids : 
                    print "Coding non-synomous: %d records." % len(snp_ids)
                    with open("%s/%s_%s_%s_rsIDs.txt" % (output_folder, gene_id, gene_symbol, protein_id), "a") as out_ids:
                        for id_snp in snp_ids :
                            out_ids.write("%s\n" % id_snp)
                    out_snp.write("%s:%s:%s:%s\n" % (gene_id, gene_symbol, protein_id, ",".join(snp_ids)))
#                     with open("%s/%s_%s_%s.xml" % (output_folder, gene_id, gene_symbol, protein_id), "a") as out_xml :
#                         sleep(1)
#                         handle = Entrez.efetch(db="snp", id=",".join(snp_ids), mode="text", report='xml')
#                         out_xml.write(handle.read())
#                         handle.close()
        # Wait a little bit between queries so we don't upset NCBI.
        sleep(1)
