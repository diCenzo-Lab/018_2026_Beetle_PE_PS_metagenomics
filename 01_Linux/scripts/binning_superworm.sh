#!/usr/bin/env bash

# sh Binning.sh threads sample_name

# Make output directory
if [ ! -d binning ]; then
  mkdir binning
fi

# Rename scaffolds
if [ ! -d binning/renamed_scaffolds ]; then
  mkdir binning/renamed_scaffolds
fi
sed "s/>/>${2}_/g" metagenome_assemblies_100/reduced_assemblies/${2}_2.5K.fasta > binning/renamed_scaffolds/${2}_2.5K.renamed.fasta

# Build bowtie2 index
if [ ! -d binning/bowtie2_indexes ]; then
  mkdir binning/bowtie2_indexes
fi
bowtie2-build --threads ${1} binning/renamed_scaffolds/${2}_2.5K.renamed.fasta binning/bowtie2_indexes/${2}

# Map all reads
if [ ! -d binning/read_mapping ]; then
  mkdir binning/read_mapping
fi
mkdir binning/read_mapping/${2}
# Repeat the two lines below for every set of paired-end reads
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Standard_1a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Standard_1a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Standard_1a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Standard_1a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Standard_1a_to_${2}.bam -o binning/read_mapping/${2}/superworm_Standard_1a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Control_1a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Control_1a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Control_1a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Control_1a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Control_1a_to_${2}.bam -o binning/read_mapping/${2}/superworm_Control_1a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_PB_1a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_PB_1a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_PB_1a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_PB_1a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_PB_1a_to_${2}.bam -o binning/read_mapping/${2}/superworm_PB_1a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Styrofoam_1a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Styrofoam_1a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Styrofoam_1a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Styrofoam_1a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Styrofoam_1a_to_${2}.bam -o binning/read_mapping/${2}/superworm_Styrofoam_1a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Standard_1b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Standard_1b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Standard_1b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Standard_1b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Standard_1b_to_${2}.bam -o binning/read_mapping/${2}/superworm_Standard_1b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Control_1b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Control_1b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Control_1b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Control_1b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Control_1b_to_${2}.bam -o binning/read_mapping/${2}/superworm_Control_1b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_PB_1b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_PB_1b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_PB_1b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_PB_1b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_PB_1b_to_${2}.bam -o binning/read_mapping/${2}/superworm_PB_1b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Styrofoam_1b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Styrofoam_1b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Styrofoam_1b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Styrofoam_1b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Styrofoam_1b_to_${2}.bam -o binning/read_mapping/${2}/superworm_Styrofoam_1b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Standard_2a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Standard_2a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Standard_2a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Standard_2a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Standard_2a_to_${2}.bam -o binning/read_mapping/${2}/superworm_Standard_2a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Control_2a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Control_2a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Control_2a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Control_2a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Control_2a_to_${2}.bam -o binning/read_mapping/${2}/superworm_Control_2a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_PB_2a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_PB_2a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_PB_2a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_PB_2a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_PB_2a_to_${2}.bam -o binning/read_mapping/${2}/superworm_PB_2a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Styrofoam_2a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Styrofoam_2a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Styrofoam_2a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Styrofoam_2a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Styrofoam_2a_to_${2}.bam -o binning/read_mapping/${2}/superworm_Styrofoam_2a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Standard_2b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Standard_2b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Standard_2b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Standard_2b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Standard_2b_to_${2}.bam -o binning/read_mapping/${2}/superworm_Standard_2b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Control_2b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Control_2b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Control_2b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Control_2b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Control_2b_to_${2}.bam -o binning/read_mapping/${2}/superworm_Control_2b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_PB_2b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_PB_2b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_PB_2b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_PB_2b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_PB_2b_to_${2}.bam -o binning/read_mapping/${2}/superworm_PB_2b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Styrofoam_2b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Styrofoam_2b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Styrofoam_2b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Styrofoam_2b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Styrofoam_2b_to_${2}.bam -o binning/read_mapping/${2}/superworm_Styrofoam_2b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Standard_3a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Standard_3a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Standard_3a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Standard_3a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Standard_3a_to_${2}.bam -o binning/read_mapping/${2}/superworm_Standard_3a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Control_3a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Control_3a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Control_3a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Control_3a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Control_3a_to_${2}.bam -o binning/read_mapping/${2}/superworm_Control_3a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_PB_3a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_PB_3a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_PB_3a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_PB_3a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_PB_3a_to_${2}.bam -o binning/read_mapping/${2}/superworm_PB_3a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Styrofoam_3a_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Styrofoam_3a_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Styrofoam_3a_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Styrofoam_3a_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Styrofoam_3a_to_${2}.bam -o binning/read_mapping/${2}/superworm_Styrofoam_3a_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Standard_3b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Standard_3b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Standard_3b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Standard_3b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Standard_3b_to_${2}.bam -o binning/read_mapping/${2}/superworm_Standard_3b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Control_3b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Control_3b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Control_3b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Control_3b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Control_3b_to_${2}.bam -o binning/read_mapping/${2}/superworm_Control_3b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_PB_3b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_PB_3b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_PB_3b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_PB_3b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_PB_3b_to_${2}.bam -o binning/read_mapping/${2}/superworm_PB_3b_to_${2}.sorted.bam
bowtie2 -p ${1} -x binning/bowtie2_indexes/${2} -1 read_preprocessing/cleaned_reads/superworm_Styrofoam_3b_R1.cleaned.fastq.gz -2 read_preprocessing/cleaned_reads/superworm_Styrofoam_3b_R2.cleaned.fastq.gz --no-unal 2> binning/read_mapping/${2}/superworm_Styrofoam_3b_to_${2}.txt | samtools view -Sb -@ 2 -o binning/read_mapping/${2}/superworm_Styrofoam_3b_to_${2}.bam
samtools sort -@ ${1} -m 3G binning/read_mapping/${2}/superworm_Styrofoam_3b_to_${2}.bam -o binning/read_mapping/${2}/superworm_Styrofoam_3b_to_${2}.sorted.bam
# End of part to be repeated
ls binning/read_mapping/${2}/*.sorted.bam | parallel samtools index '{}'
rm binning/read_mapping/${2}/*a.bam
rm binning/read_mapping/${2}/*b.bam

# Generate MetaBat2 bins
if [ ! -d binning/depth ]; then
  mkdir binning/depth
fi
if [ ! -d binning/metabat2_bins ]; then
  mkdir binning/metabat2_bins
fi
jgi_summarize_bam_contig_depths --outputDepth binning/depth/${2}_depth.txt binning/read_mapping/${2}/*.sorted.bam
metabat2 -i binning/renamed_scaffolds/${2}_2.5K.renamed.fasta -a binning/depth/${2}_depth.txt -t ${1} -o binning/metabat2_bins/${2}/
rename '.' '' binning/metabat2_bins/${2}/.*

# Generate CONCOCT bins
if [ ! -d binning/cut_scaffolds ]; then
  mkdir binning/cut_scaffolds
fi
if [ ! -d binning/coverage ]; then
  mkdir binning/coverage
fi
if [ ! -d binning/concoct_binning ]; then
  mkdir binning/concoct_binning
fi
mkdir binning/concoct_binning/${2}
if [ ! -d binning/concoct_bins ]; then
  mkdir binning/concoct_bins
fi
mkdir binning/concoct_bins/${2}
cut_up_fasta.py binning/renamed_scaffolds/${2}_2.5K.renamed.fasta -c 10000 -o 0 --merge_last -b binning/cut_scaffolds/${2}_2.5K.renamed_10K.bed > binning/cut_scaffolds/${2}_2.5K.renamed_10K.fa
concoct_coverage_table.py binning/cut_scaffolds/${2}_2.5K.renamed_10K.bed binning/read_mapping/${2}/*.sorted.bam > binning/coverage/${2}_concoct_coverage_table.tsv
concoct --composition_file binning/cut_scaffolds/${2}_2.5K.renamed_10K.fa --coverage_file binning/coverage/${2}_concoct_coverage_table.tsv -b binning/concoct_binning/${2}/${2} -t ${1}
merge_cutup_clustering.py binning/concoct_binning/${2}/${2}_clustering_gt1000.csv > binning/concoct_binning/${2}/${2}_clustering_merged.csv
extract_fasta_bins.py binning/renamed_scaffolds/${2}_2.5K.renamed.fasta binning/concoct_binning/${2}/${2}_clustering_merged.csv --output_path binning/concoct_bins/${2}/

# Generate MaxBin2 bins
if [ ! -d binning/abundance ]; then
  mkdir binning/abundance
fi
mkdir binning/abundance/${2}
if [ ! -d binning/maxbins2_bins ]; then
  mkdir binning/maxbins2_bins
fi
mkdir binning/maxbins2_bins/${2}
for N in 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50
do
  NAME=$(cut -f${N},${N} binning/depth/${2}_depth.txt | sed 's/\.sorted\.bam//' | head -1)
  awk -v i="$N" '{print $1"\t"$i}' binning/depth/${2}_depth.txt | sed "s/$NAME.sorted.bam/abundance/g" | sed 's/contigName/contig header/g' > binning/abundance/${2}/${NAME}_abundance.txt
done
find binning/abundance/${2}/*abundance.txt > binning/abundance/${2}/abundance_list.txt
run_MaxBin.pl -contig binning/renamed_scaffolds/${2}_2.5K.renamed.fasta -out binning/maxbins2_bins/${2}/${2} -abund_list binning/abundance/${2}/abundance_list.txt -thread ${1}

# Get final bins with DAS_Tool
if [ ! -d binning/DAS_Tool_input ]; then
  mkdir binning/DAS_Tool_input
fi
mkdir binning/DAS_Tool_input/${2}
if [ ! -d binning/DAS_Tool_output ]; then
  mkdir binning/DAS_Tool_output
fi
mkdir binning/DAS_Tool_output/${2}
Fasta_to_Contig2Bin.sh -i binning/metabat2_bins/${2}/ -e fa > binning/DAS_Tool_input/${2}/${2}_Metabat_bins.tsv
Fasta_to_Contig2Bin.sh -i binning/concoct_bins/${2}/ -e fa > binning/DAS_Tool_input/${2}/${2}_concoct_bins.tsv
Fasta_to_Contig2Bin.sh -i binning/maxbins2_bins/${2}/ -e fasta > binning/DAS_Tool_input/${2}/${2}_Maxbin_bins.tsv
DAS_Tool -i binning/DAS_Tool_input/${2}/${2}_Metabat_bins.tsv,binning/DAS_Tool_input/${2}/${2}_concoct_bins.tsv,binning/DAS_Tool_input/${2}/${2}_Maxbin_bins.tsv -l concoct,maxbin2,metabat2 -c binning/renamed_scaffolds/${2}_2.5K.renamed.fasta --threads ${1} --write_bins -o binning/DAS_Tool_output/${2}/${2}

# Save space by removing BAM files
#rm binning/read_mapping/${2}/*.bam

