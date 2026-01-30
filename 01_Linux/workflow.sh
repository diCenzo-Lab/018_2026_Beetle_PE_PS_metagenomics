# Preprocess illumina reads
mkdir read_preprocessing/ # create directory to hold all read preprocessing data
mkdir read_preprocessing/bbduk_output/ # create directory to hold bbduk cleaned reads
sh scripts/bbduk_run.sh # perform initial QC filtering of illumina reads
mkdir read_preprocessing/trimmomatic_output/ # create directory to hold trimmed reads
sh scripts/trimmomatic.sh # trim the reads with trimmomatic
mkdir read_preprocessing/insect_genomes/ # create directory to hold insect genome bowtie2 databases
cp ../insect_genomes/final_genomes/mealworm_masked.fasta read_preprocessing/insect_genomes/mealworm.fasta # get mealworm genome
cp ../insect_genomes/final_genomes/superworm_masked.fasta read_preprocessing/insect_genomes/superworm.fasta # get superworm genome
sed -i 's/a/A/g' read_preprocessing/insect_genomes/mealworm.fasta # unmask, just in case
sed -i 's/t/T/g' read_preprocessing/insect_genomes/mealworm.fasta # unmask, just in case
sed -i 's/g/G/g' read_preprocessing/insect_genomes/mealworm.fasta # unmask, just in case
sed -i 's/c/C/g' read_preprocessing/insect_genomes/mealworm.fasta # unmask, just in case
sed -i 's/a/A/g' read_preprocessing/insect_genomes/superworm.fasta # unmask, just in case
sed -i 's/t/T/g' read_preprocessing/insect_genomes/superworm.fasta # unmask, just in case
sed -i 's/g/G/g' read_preprocessing/insect_genomes/superworm.fasta # unmask, just in case
sed -i 's/c/C/g' read_preprocessing/insect_genomes/superworm.fasta # unmask, just in case
bowtie2-build --threads 16 read_preprocessing/insect_genomes/mealworm.fasta read_preprocessing/insect_genomes/mealworm # build bowtie2 database
bowtie2-build --threads 16 read_preprocessing/insect_genomes/superworm.fasta read_preprocessing/insect_genomes/superworm # build bowtie2 database
mkdir read_preprocessing/insect_alignment/ # create directory to hold insect mapping files
sh scripts/bowtie2.sh # map all reads to the corresponding insect genome
sh scripts/extract_unmapped_reads.sh # extract cleaned reads lacking reads that map to the insect genome assemblies
mkdir read_preprocessing/stats/ # create directory to hold information on the sequencing stats
sh scripts/generate_read_stats.sh # get some stats for the illumina data

# Prepare metagenome assemblies
mkdir metagenome_assemblies_100/ # create directory to hold all assembly data
sh scripts/spades.sh # run metaSPAdes to create metagenome assemblies
mkdir metagenome_assemblies_100/stats/ # create directory for assembly stats
sh scripts/stats.sh # get assembly statistics
mkdir metagenome_assemblies_100/reduced_assemblies/ # create directory for metagenome assemblies (contigs ≥ 2.5 kb)
sh scripts/reduced_assemblies.sh # create reduced assemblies consisting only of contigs ≥ 2.5 kb
mkdir metagenome_assemblies_100/stats_reduced/ # create directory for reduced assembly stats
sh scripts/stats_reduced.sh # get reduced assembly statistics

# Bin contigs into prokaryotic MAGs
mkdir binning/ # create directory to hold all binning data
sh scripts/run_binning.sh # bin the samples with MetaBat2, MaxBin2, CONCOCT, and DAS_Tool
sh scripts/dereplicate.sh # dereplicate the MAGs with dRep, and rename
# Final MAGs are found at: binning/final_MAGs/

# Get eukaryotic MAGs
sh scripts/euk_bins.sh # get eukaryotic bins from the already binned data
# Final MAGs are found at: binning/final_MAGs/
mkdir binning/busco_euk/ # create directory to hold final busco data for the one eukaryotic MAG
busco -i binning/final_MAGs/mealworm/Mealworm_eMAG_0001.fasta -o binning/busco_euk/eukaryota/ -l eukaryota_odb10 -c 16 -m genome -f
busco -i binning/final_MAGs/mealworm/Mealworm_eMAG_0001.fasta -o binning/busco_euk/auto/ --auto-lineage-euk -c 16 -m genome -f

# Get my own CheckM scores of the MAGs
mkdir CheckM/
mkdir CheckM/MAGs/
cp binning/final_MAGs/*/*.fasta CheckM/MAGs/
checkm lineage_wf -t 16 -x fasta -f checkm_output.txt CheckM/MAGs/ CheckM/
mv checkm_output.txt CheckM/

# Taxonomically classify MAGs
mkdir GTDB_classification/ # create directory to hold the GTDB-tk output
mkdir GTDB_classification/mealworm/ # create directory to hold the mealworm output
mkdir GTDB_classification/superworm/ # create directory to hold the superworm output
mkdir GTDB_classification/tmp_dir/ # Make a termpoary directory
gtdbtk classify_wf --genome_dir binning/final_MAGs/mealworm/ --out_dir GTDB_classification/mealworm/ --cpus 16 --extension fasta --prefix Mealworm_MAGs --pplacer_cpus 16 --tmpdir GTDB_classification/tmp_dir/ # Run GTDB-tk for the mealworm MAGs
gtdbtk classify_wf --genome_dir binning/final_MAGs/superworm/ --out_dir GTDB_classification/superworm/ --cpus 16 --extension fasta --prefix Superworm_MAGs --pplacer_cpus 12 --tmpdir GTDB_classification/tmp_dir/ # Run GTDB-tk for the superworm MAGs

# Create a maximum likelihood phylogeny of the bacterial MAGs, but first rename the files
mkdir bacterial_phylogeny/
mkdir bacterial_phylogeny/MAGs/
cp binning/final_MAGs/*/*.fasta bacterial_phylogeny/MAGs/
rm bacterial_phylogeny/MAGs/Mealworm_eMAG_0001.fasta
perl scripts/add_species_to_MAGs.pl GTDB_classification/mealworm/Mealworm_MAGs.bac120.summary.tsv
perl scripts/add_species_to_MAGs.pl GTDB_classification/superworm/Superworm_MAGs.bac120.summary.tsv
mv bacterial_phylogeny/MAGs/Bacillaceae_G_sp_Mealworm_MAG_0052.fasta bacterial_phylogeny/MAGs/Bacillaceae_bacterium_Mealworm_MAG_0052.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Mealworm_MAG_0005.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Mealworm_MAG_0005.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Mealworm_MAG_0024.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Mealworm_MAG_0024.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Mealworm_MAG_0038.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Mealworm_MAG_0038.fasta
mv bacterial_phylogeny/MAGs/JALAGN01_sp_Mealworm_MAG_0031.fasta bacterial_phylogeny/MAGs/Mycoplasmataceae_bacterium_Mealworm_MAG_0031.fasta
mv bacterial_phylogeny/MAGs/Burkholderiales_sp_Superworm_MAG_0027.fasta bacterial_phylogeny/MAGs/Burkholderiales_bacterium_Superworm_MAG_0027.fasta
mv bacterial_phylogeny/MAGs/Caloramatoraceae_sp_Superworm_MAG_0052.fasta bacterial_phylogeny/MAGs/Caloramatoraceae_bacterium_Superworm_MAG_0052.fasta
mv bacterial_phylogeny/MAGs/Enterobacteriaceae_sp_Superworm_MAG_0037.fasta bacterial_phylogeny/MAGs/Enterobacteriaceae_bacterium_Superworm_MAG_0037.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Superworm_MAG_0009.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Superworm_MAG_0009.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Superworm_MAG_0020.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Superworm_MAG_0020.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Superworm_MAG_0040.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Superworm_MAG_0040.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Superworm_MAG_0041.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Superworm_MAG_0041.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Superworm_MAG_0060.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Superworm_MAG_0060.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Superworm_MAG_0063.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Superworm_MAG_0063.fasta
mv bacterial_phylogeny/MAGs/Enterococcaceae_sp_Superworm_MAG_0095.fasta bacterial_phylogeny/MAGs/Enterococcus_sp_Superworm_MAG_0095.fasta
mv bacterial_phylogeny/MAGs/Erysipelotrichaceae_sp_Superworm_MAG_0017.fasta bacterial_phylogeny/MAGs/Erysipelotrichaceae_bacterium_Superworm_MAG_0017.fasta
mv bacterial_phylogeny/MAGs/JAJDGJ01_sp_Superworm_MAG_0016.fasta bacterial_phylogeny/MAGs/Bacilli_bacterium_Superworm_MAG_0016.fasta
mv bacterial_phylogeny/MAGs/Mycoplasmataceae_sp_Superworm_MAG_0007.fasta bacterial_phylogeny/MAGs/Mycoplasmataceae_bacterium_Superworm_MAG_0007.fasta
mv bacterial_phylogeny/MAGs/Mycoplasmataceae_sp_Superworm_MAG_0093.fasta bacterial_phylogeny/MAGs/Mycoplasmataceae_bacterium_Superworm_MAG_0093.fasta
mv bacterial_phylogeny/MAGs/Peptostreptococcaceae_sp_Superworm_MAG_0018.fasta bacterial_phylogeny/MAGs/Romboutsia_sp_Superworm_MAG_0018.fasta
mv bacterial_phylogeny/MAGs/Peptostreptococcaceae_sp_Superworm_MAG_0019.fasta bacterial_phylogeny/MAGs/Romboutsia_sp_Superworm_MAG_0019.fasta
mv bacterial_phylogeny/MAGs/Streptococcaceae_sp_Superworm_MAG_0062.fasta bacterial_phylogeny/MAGs/Lactococcus_sp_Superworm_MAG_0062.fasta
rename 's/Superworm_MAG/sMAG/' bacterial_phylogeny/MAGs/*.fasta
rename 's/Mealworm_MAG/mMAG/' bacterial_phylogeny/MAGs/*.fasta

# Now make the phylogeny
gtdbtk identify --genome_dir bacterial_phylogeny/MAGs/ --out_dir bacterial_phylogeny/markers/ --extension fasta --cpus 16 --write_single_copy_genes
gtdbtk align --identify_dir bacterial_phylogeny/markers/ --out_dir bacterial_phylogeny/aligned_markers/ --skip_gtdb_refs --rnd_seed 97 --cpus 16
cd bacterial_phylogeny/
cp aligned_markers/align/gtdbtk.bac120.user_msa.fasta.gz .
gunzip gtdbtk.bac120.user_msa.fasta.gz
iqtree2 -s gtdbtk.bac120.user_msa.fasta -m MF -T 16 --prefix model
best_model=$(grep 'Best-fit' model.log | cut -f3 -d' ') # LG+F+R8
iqtree2 -s gtdbtk.bac120.user_msa.fasta -m $best_model --alrt 1000 -J 1000 --jack-prop 0.4 -T 16 --prefix Phylogeny

# Cluster MAGs into species with ANI
mkdir MAG_clustering/ # create directory to hold the MAG clustering data
sh scripts/cluster_MAGs.sh # run FastANI and then use the data to form matrices

# Get MAG abundance
mkdir MAG_abundance/ # create directory to hold the MAG abundance data
sh scripts/get_MAG_abundances.sh # map the illumina reads to the MAGs

# Annotated metagenomes and MAGs
# Annotation with DRAM done using KBase
mkdir antiSMASH_output/ # make directory ot hold the antiSMASH output
cp -r bacterial_phylogeny/MAGs/ antiSMASH_output/
sed -i 's/Superworm_/s/' antiSMASH_output/MAGs/*sMAG*
sed -i 's/Mealworm_/m/' antiSMASH_output/MAGs/*mMAG*
conda activate antismash # activate the antiSMASH python environment
sh scripts/run_antismash_bacterial.sh # run antiSMASH on the bacterial MAGs
conda deactivate
grep -m 1 'category' antiSMASH_output/*/*region*.gbk | sed 's/  / /g' | sed 's/  / /g' | sed 's/  / /g' | sed 's/  / /g' | sed 's/  / /g' > temp1.txt # Get the classification of all BGCs in all MAGs
grep -m 1 '/product' antiSMASH_output/*/*region*.gbk | sed 's/  / /g' | sed 's/  / /g' | sed 's/  / /g' | sed 's/  / /g' | sed 's/  / /g' > temp2.txt # Get the product of all BGCs in all MAGs
join temp1.txt temp2.txt > temp3.txt # Create one file that has both the categories and the product information
cut -f2 -d' ' temp3.txt | sed 's/\/category=//' | sed 's/\"//g' > temp4.txt # Get a clean column of just the categorization
cut -f1 -d' ' temp3.txt | sed "s/\//\t/g" | sed 's/.gbk://' | cut -f2,3 > temp5.txt # Get columns of the MAG name and BGC locus
cut -f3 -d' ' temp3.txt | sed 's/\/product=//' | sed 's/\"//g' > temp6.txt # Get a clean column of just the product
paste temp5.txt temp4.txt temp6.txt > antiSMASH_output/summary.csv # Create a summary file
cut -f1,3,4 antiSMASH_output/summary.csv | sort -u | cut -f2,3 | sort | uniq -c > antiSMASH_output/summary.statistics.csv # Summarize the number of clusters
rm temp*.txt # Remove temporary files

# Get secretome and find protein families of interest
mkdir secretome/ # make directory to hold the secretome files
mkdir secretome/secreted_proteins/ # make directory for holding the total secretome
mkdir secretome/secreted_proteins/metagenome_proteomes/ # Make folder to hold the DRAM annotations as protein fasta files
# Manually upload the DRAM outputs to the above folder
cat secretome/secreted_proteins/metagenome_proteomes/*.faa > secretome/secreted_proteins/total_metagenome_proteins.fasta # Combine all the proteome files
mkdir secretome/secreted_proteins/cdhit_output/ # Folder to hold CD-HIT output
mkdir secretome/secreted_proteins/signalp_output/ # Folder to hold signalP output
cd-hit -T 28 -aL 0.99 -c 0.99 -M 50000 -d 0 -i secretome/secreted_proteins/total_metagenome_proteins.fasta -o secretome/secreted_proteins/cdhit_output/total_nonredundant_proteins_full # Get a non-redundant protein set
mv secretome/secreted_proteins/cdhit_output/total_nonredundant_proteins_full secretome/secreted_proteins/cdhit_output/total_nonredundant_proteins_full.fasta # Rename file
source /datadisk1/Bioinformatics_programs/signalp6_fast/sp6env/bin/activate # Activate virtual environment for signalP
signalp6 --fastafile secretome/secreted_proteins/cdhit_output/total_nonredundant_proteins_full.fasta --output_dir secretome/secreted_proteins/signalp_output --format none --organism other --mode fast --write_procs 8 --bsize 40 # Predict the secreted proteins
deactivate # Deactivate virtual environment for signalP
cp secretome/secreted_proteins/signalp_output/processed_entries.fasta secretome/full_secretome.fasta # Make a copy of the predicted secreted protein fasta file
mkdir secretome/oxidative_proteins/ # make directory for holding the predicted oxidative proteins
mkdir secretome/oxidative_proteins/input_data/ # Make directory to hold the input data
# Manually add the seed alignments for the following protein families to the above folder: PF00394 (laccases), PF07731 (laccases), PF07732 (laccases), PF00255 (glutathione peroxidases), PF04261 (Dyp-type peroxidases), PF21105 (Dyp-type peroxidases), PF20628 (Dyp-type peroxidases), PF13532 (alkane-1 monooxygenases), PF00487 (fatty acid desaturases)
# Manually add a protein multifasta file to the above folder, that contains the amino acid sequences of eight proteins previously implicated in plastic degradation
sh scripts/get_oxidative_proteins.sh # get putatitve orthologs of the protein families of interest from the secreted proteins (involves a manual step)
# Created SSNs using https://efi.igb.illinois.edu/efi-est/

# Blast secreted proteins of interest against MAGs to get top hits
mkdir secretome/blast # make directory for holding the blast output
cat bacterial_phylogeny/MAGs/*.fasta > secretome/blast/combined_MAGs.fasta # Combine the MAGs as one multifasta file
cat secretome/oxidative_proteins/hmmscan_output/*.fasta > secretome/blast/secreted_proteins.fasta # Combine the secreted proteins of interest as a single multifasta file
cd secretome/blast/ # Change directory
makeblastdb -in combined_MAGs.fasta -out combined_MAGs -title combined_MAGs -dbtype 'nucl' # Make the blast db
tblastn -query secreted_proteins.fasta -db combined_MAGs -out blast.output.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand' -culling_limit 1 -max_target_seqs 1 -max_hsps 1 -num_threads 16 # Run tblastn