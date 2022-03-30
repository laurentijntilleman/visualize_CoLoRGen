Extra scripts for visualizing CoLoRGen output
=============================================

CoLoRGen: comprehensive long read genotyping pipeline.
This tool generates consensus sequence of the targeted regio per allele
and returns vcf files for each gene that is found. If star-alleles are
present, the correct star-allles are assigned to the genes.

The scripts in this repository expect that the [CoLoRGen](https://github.com/laurentijntilleman/CoLoRGen) is already ran.


Installation
------------

Most of the script use python3, therefore python3 need to be installed on the
system. Also the virtualenv package need to be present as CoLoRGen runs
everything in a virtual environment.

Usage
-----

To visulize the data run

    visualize.sh -p <parameters.sh>

The following paramters need to be adjust in the `parameters.sh` file:

    number_treats=40
    nanopore_sequencing_files='fastq_pass'  # raw fastq files
    analyse_folder='./' # working directory
    ref_folder="ref/" # reference folder
    reference_genome='GRCh38.fa' # fasta file of the reference genome
    minimap2_index='GRCh38.fa.mmi' # minimap2 index file
    reference_gtf='GRCh38.gtf' # gtf file with the exon regions
    # edit paths to external software tools
    export PATH=samtools-1.11:$PATH
    export PATH=bcftools-1.11:$PATH
    export PATH=htslib-1.11:$PATH
    hap_py='hap.py' #path to the hap.py script
    gene_json='cyp2d6_hybrids.json' # json file with gene info and hybrid info
    star_alleles="${ref_folder}CYP2D6.NC_000022.11.haplotypes.tsv" # star allele nomenclature downloaded from PharmVar
    # position of the targeted region
    target_start=42121165
    target_stop=42149371
    target_chr='chr22'
    # medaka model
    model_m='r941_min_hac_g507' # model for medaka variant calling and consensus
    log_level='debug' # log level

In the git repository, there is a example of the json file with the gene info

The references set composed by Krushe *et al.* \[1\] of the GM12878 cell line is added to the ref folder.


Output
------

The figures will be outputed in the haplotype subfolder of the output folder of the CoLoRGen pipeline


## References

\[1\] Krusche P, Trigg L, Boutros PC, Mason CE, De La Vega FM, Moore BL, *et al.*, *Best practices for benchmarking germline small-variant calls in human genomes.*, Nat Biotechnol. 2019;37(5):555–60. 
