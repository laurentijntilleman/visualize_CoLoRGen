Help()
{
    echo 'usage: base_script.sh [-h] [-c]'
    echo
    echo 'Visualize the variants generated with the CoLoRGen pipeline of the GM12878 cell line'
    echo
    echo 'arguments:'
    echo '  -p  path to the parameters file'
    echo
    echo 'optinal arguments:'
    echo '  -h  show this help message and exit'
    exit 1
}

while getopts p:hd flag
do
    case "${flag}" in
        p) parameters_path=${OPTARG};;
        h) Help;;
        d) debug=true
    esac
done

if [ -z ${parameters_path} ]; then
  echo 'Error no parameter -p'
  echo
  Help
fi

source ${parameters_path}

FILTER="(((TYPE=\"snp\" | TYPE=\"mnp\") & (QUAL > 9)) | ((TYPE!=\"snp\" & TYPE!=\"mnp\") & (QUAL > 9)))"

haplotype_folder=${analyse_folder}"haplotype/" #haplotype data
log_folder=${analyse_folder}"log/" #log folder
gene_folder=${analyse_folder}"gene/" #log folder

## start logging
echo "$(date) start script"
default_log_file=${log_folder}stout.log
error_log_file=${log_folder}error.log

# start virtualenv
source venv/bin/activate

sample=${haplotype_folder}consensus_3
bowtie2 -a -x ${sample} -p ${number_treats} -f ${ref_folder}/gRNAs.fasta -S ${sample}_gRNA.sam
samtools view -Sb ${sample}_gRNA.sam > ${sample}_gRNA.bam
samtools sort ${sample}_gRNA.bam > ${sample}_gRNA_sorted.bam
samtools index ${sample}_gRNA_sorted.bam
bedtools bamtobed -i ${sample}_gRNA_sorted.bam > ${sample}_gRNA_sorted.bed

# define haplotype
if [[ "`ls -dq ${haplotype_folder}haplotypes_*_CYP2D6.vcf | wc -l`" == 2 ]]; then
  medaka tools haploid2diploid ${haplotype_folder}haplotypes_*_CYP2D6.vcf ${reference_genome} ${haplotype_folder}haplotypes_CYP2D6.vcf --split_mnp
  bcftools filter --threads ${number_treats}  -s lowqual -i "${FILTER}" ${haplotype_folder}haplotypes_CYP2D6.vcf > ${haplotype_folder}haplotypes_CYP2D6_filtered.vcf
  bgzip ${haplotype_folder}haplotypes_CYP2D6_filtered.vcf
  tabix ${haplotype_folder}haplotypes_CYP2D6_filtered.vcf.gz
  python2 ${hap_py} ${ref_folder}/hg38.hybrid.vcf.gz ${haplotype_folder}haplotypes_CYP2D6_filtered.vcf.gz --threads 2 -r ${reference_genome} -R ${gene_folder}CYP2D6.bed  -o ${haplotype_folder}haplotypes_hap_CYP2D6
  python visualise.py --log_level ${log_level}  --log_file ${default_log_file} --vcf ${haplotype_folder}haplotypes_CYP2D6_filtered.vcf.gz -r ${ref_folder}/hg38.hybrid.vcf.gz --start ${CYP2D6_start} --stop ${CYP2D6_stop} -c ${CYP2D6_chr} -f ${haplotype_folder} --name 'CYP2D6'
fi

if [[ "`ls -dq ${haplotype_folder}haplotypes_*_CYP2D7.vcf | wc -l`" == 2  ]]; then
  medaka tools haploid2diploid ${haplotype_folder}haplotypes_*_CYP2D7.vcf ${reference_genome} ${haplotype_folder}haplotypes_CYP2D7.vcf --split_mnp
  bcftools filter --threads ${number_treats}  -s lowqual -i "${FILTER}" ${haplotype_folder}haplotypes_CYP2D7.vcf > ${haplotype_folder}haplotypes_CYP2D7_filtered.vcf
  bgzip ${haplotype_folder}haplotypes_CYP2D7_filtered.vcf
  tabix ${haplotype_folder}haplotypes_CYP2D7_filtered.vcf.gz
  python2 ${hap_py} ${ref_folder}/hg38.hybrid.vcf.gz ${haplotype_folder}haplotypes_CYP2D7_filtered.vcf.gz --threads 2 -r ${reference_genome} -R ${gene_folder}CYP2D7.bed  -o ${haplotype_folder}haplotypes_hap_CYP2D7
  python visualise.py --log_level ${log_level}  --log_file ${default_log_file} --vcf ${haplotype_folder}haplotypes_CYP2D7_filtered.vcf.gz -r ${ref_folder}/hg38.hybrid.vcf.gz --start ${CYP2D7_start} --stop ${CYP2D7_stop} -c ${CYP2D7_chr} -f ${haplotype_folder} --name 'CYP2D7'
fi

bgzip ${haplotype_folder}medaka_round_0/round_1_phased.vcf
tabix ${haplotype_folder}medaka_round_0/round_1_phased.vcf.gz

python2 ${hap_py} ${ref_folder}/hg38.hybrid.vcf.gz ${haplotype_folder}medaka_round_0/round_1_phased.vcf.gz --threads 2 -r ${reference_genome} -R ${gene_folder}target.bed  -o ${haplotype_folder}medaka_round_0_hap

python2 ${hap_py} ${ref_folder}/hg38.hybrid.vcf.gz ${haplotype_folder}medaka_round_0/round_1_phased.vcf.gz  --threads 2 -r ${reference_genome} -R ${gene_folder}CYP2D6.bed  -o ${haplotype_folder}medaka_round_0_hap_CYP2D6
python2 ${hap_py} ${ref_folder}/hg38.hybrid.vcf.gz ${haplotype_folder}medaka_round_0/round_1_phased.vcf.gz  --threads 2 -r ${reference_genome} -R ${gene_folder}CYP2D7.bed  -o ${haplotype_folder}medaka_round_0_hap_CYP2D7
python visualise.py --log_level ${log_level}  --log_file ${default_log_file} --vcf ${haplotype_folder}medaka_round_0/round_1_phased.vcf.gz -r ${ref_folder}/hg38.hybrid.vcf.gz --start ${CYP2D6_start} --stop ${CYP2D6_stop} -c ${CYP2D6_chr} -f ${haplotype_folder} --name 'medaka_CYP2D6'
python visualise.py --log_level ${log_level}  --log_file ${default_log_file} --vcf ${haplotype_folder}medaka_round_0/round_1_phased.vcf.gz -r ${ref_folder}/hg38.hybrid.vcf.gz --start ${CYP2D7_start} --stop ${CYP2D7_stop} -c ${CYP2D7_chr} -f ${haplotype_folder} --name 'medaka_CYP2D7'

# end virtualenv
deactivate
