# Variants of pathway genes
## This is a collection of scripts used to identify and download all known mutations in the dbSNP, ClinVar and OMIM databases affecting genes/proteins in the human insulin signaling and tryptophan metabolism pathways. All genomic locations refer to the GRCh38.p7 assembly.

	Prerequisites:
	- python 2.7
	- BioPython (pip install biopython)
	- PyVCF 	 (pip install pyvcf)
	- vcftools		https://vcftools.github.io/index.html
	- zlib module https://docs.python.org/3/library/zlib.html#module-zlib
	- a SSH tunnel or the Tor browser (https://www.torproject.org/) to 
	  crawl the OMIM web page without getting your IP banned.

	Usage:

	The scripts are meant to be executed in a certain order. Folders and subfolders 
	will be created automatically. Here I will describe each script in the intended 
	order of execution:

		1- retrieve_and_parse_KGML_pathways.py
		This script retrieves all known genes/proteins involved in the human 
		'insulin signaling' and 'tryptophan metabolism' pathways from the KEGG
		database. The results are saved in the 
		'./insulin_signaling/insulin_signaling_genes.txt' 
		and in the
		'./tryptophan_metabolism/tryptophan_metabolism_genes.txt' 
		files respectively.
		
		2- get_gene_names_and_symbols.py
		This script loads the list of genes from the previously computed files
		and retrieves, for each gene, the official gene symbol according to the 
		HUGO Gene 	Nomenclature Committee: HGNC database of human gene names.
		The output is saved in the 'official_gene_symbols_and_names.txt' file.
		The entries in the 'manually_added_genes.txt' files are added to the end
		of the output file. These genes were extracted from the known models of
		the tryptophan metabolism and insulin signaling pathways 
		(FBA, RBA, dynamic models).
		
		3- fetch_all_cds.py
		This script retrieves the CDSs (Coding Region Sequences) for each gene.
		
		4- fetch_RefSeq_cds.py
		This script retrieves the CDSs (Coding Region Sequences) for each RefSeqGene.
		
		5- query_dbSNP_NCBI.py
		This script queries the NCBI dbSNP database for coding non-synomous SNPs for 
		a given set of proteins/genes. The database can be searched programmatically 
		using the Entrez module from the BioPython package.
		
		6- extract_dbSNP_VCF.sh
		This script extracts the snp records in VCF format from the complete dbSNP Human 
		database. The script uses multiprocessing to speed-up the process. The number of 
		parallel processes can be changed by setting the variable N at line 29.
		For the first usage: uncomment line 27 to download the current dbSNP Human database 
		in VCF format. 
		(i.e. #wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz)
		NOTE: The complete human dbSNP database in VCF format '00-All.vcf.gz' is a very
		large file: +7GB compressed and +50GB uncompressed. Parsing it takes quite a lot 
		of time, so it is a good idea to use parallelism and multithreading. However, 
		keep in mind that too many thread could end-up using too much memory and you
		probably want to avoid ending up using the swap/paging memory areas. 
		I advise to check the available RAM and set the N variable accordingly. 12 parallel
		threads should be quite safe for a computer with 4GB of RAM.
		
		7- parse_ClinVar_variants.py
		This script retrieves the latest version of the ClinVar variant summary database,
		available at:
		'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz'.
		Variations are extracted for each gene of interest. Each candidate variation is
		matched agains the known CDS coordinates for all genes of interest.
		
		8- extract_ClinVar_VCF.py
		This script retrieves the latest version of the ClinVar database in VCF format,
		available at:
		'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz' 
		Variations are extracted for each gene of interest. However, only the allele ids 
		extracted during step 7 are considered. 
		
		9- fetch_OMIM_allelic_variants.py
		This script queries the OMIM database for known mutations affecting a given set 
		of proteins/genes. Entris for each NCBI gene are present in the OMIM database,
		however, only few have annotated allelic variants. The allelic variants, when
		present, are always located in urls like: 
		https://www.omim.org/allelicVariant/<gene-ID>?format=tsv.
		The current script is a simple web crawler that tries to download all the allelic
		variants of interest from the OMIM database.
		You can also check the files in the './OMIM/' folder, especially 'genemap2.txt',
		which contains the list of all genes in the OMIM database. These files were obtained
		by written request to the OMIM administrators, and should not be made publicly 
		available, otherwise we would be breaking their license agreement. 
		
		10- parse_OMIM_variants.py
		Very similar to number 8. This script extracts the variants of interest for each
		gene/protein from the OMIM allelic variants computed in the previous step.
		
		11-extract_OMIM_VCF.py
		This script retrieves the latest version of the ClinVar database in VCF format,
		available at:
		'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz' 
		Variations are extracted for each gene of interest. However, only the allele ids 
		extracted during step 10 are considered. 
		
		12- build_table.py
		Ugly script that builds a table with all the obtained genetic variations.	
