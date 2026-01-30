# Get lists of MAGs
find binning/final_MAGs/mealworm/Mealworm_MAG*.fasta > MAG_clustering/genomePath_mealworm.txt # Get the genome paths for mealworms
find binning/final_MAGs/superworm/Superworm_MAG*.fasta > MAG_clustering/genomePath_superworm.txt # Get the genome paths for superworms
cat MAG_clustering/genomePath_mealworm.txt MAG_clustering/genomePath_superworm.txt > MAG_clustering/genomePath_all.txt # Combine the genome paths for both species

# Run FastANI
fastANI --ql MAG_clustering/genomePath_mealworm.txt --rl MAG_clustering/genomePath_mealworm.txt -o MAG_clustering/mealworm_fastANI_output.txt -t 16 # Run fastANI on mealworm MAGs
fastANI --ql MAG_clustering/genomePath_superworm.txt --rl MAG_clustering/genomePath_superworm.txt -o MAG_clustering/superworm_fastANI_output.txt -t 16 # Run fastANI on superworm MAGs
fastANI --ql MAG_clustering/genomePath_all.txt --rl MAG_clustering/genomePath_all.txt -o MAG_clustering/all_fastANI_output.txt -t 16 # Run fastANI on all MAGs

# Parse FastANI output for mealworms
sort -k1,1 -k2,2 MAG_clustering/mealworm_fastANI_output.txt > MAG_clustering/mealworm_fastANI_output_sorted_1.txt # Sort the file by first column then by second column
sort -k2,2 -k1,1 MAG_clustering/mealworm_fastANI_output.txt > MAG_clustering/mealworm_fastANI_output_sorted_2.txt # Sort the file by second column then by first column
cut -f 1,2,3 MAG_clustering/mealworm_fastANI_output_sorted_1.txt > temp.txt # Get the relevant columns of the first sorted file
cut -f 3 MAG_clustering/mealworm_fastANI_output_sorted_2.txt > temp2.txt # Get the relevant columns of the second sorted file
paste temp.txt temp2.txt > MAG_clustering/mealworm_fastANI_output_twoWay.txt # Combine the relevant columns
rm temp* # remove the temporary files
perl scripts/prepareANImatrix.pl MAG_clustering/mealworm_fastANI_output_twoWay.txt MAG_clustering/genomePath_mealworm.txt | sed 's/binning\/final_MAGs\/mealworm\///g' | sed 's/binning\/final_MAGs\/superworm\///g' | sed 's/\.fasta//g' > MAG_clustering/mealworm_ANI_matrix.txt # make a two-way ANI matrix from the fastANI output

# Parse FastANI output for superworms
sort -k1,1 -k2,2 MAG_clustering/superworm_fastANI_output.txt > MAG_clustering/superworm_fastANI_output_sorted_1.txt # Sort the file by first column then by second column
sort -k2,2 -k1,1 MAG_clustering/superworm_fastANI_output.txt > MAG_clustering/superworm_fastANI_output_sorted_2.txt # Sort the file by second column then by first column
cut -f 1,2,3 MAG_clustering/superworm_fastANI_output_sorted_1.txt > temp.txt # Get the relevant columns of the first sorted file
cut -f 3 MAG_clustering/superworm_fastANI_output_sorted_2.txt > temp2.txt # Get the relevant columns of the second sorted file
paste temp.txt temp2.txt > MAG_clustering/superworm_fastANI_output_twoWay.txt # Combine the relevant columns
rm temp* # remove the temporary files
perl scripts/prepareANImatrix.pl MAG_clustering/superworm_fastANI_output_twoWay.txt MAG_clustering/genomePath_superworm.txt | sed 's/binning\/final_MAGs\/mealworm\///g' | sed 's/binning\/final_MAGs\/superworm\///g' | sed 's/\.fasta//g' > MAG_clustering/superworm_ANI_matrix.txt # make a two-way ANI matrix from the fastANI output

# Parse FastANI output for all MAGs
sort -k1,1 -k2,2 MAG_clustering/all_fastANI_output.txt > MAG_clustering/all_fastANI_output_sorted_1.txt # Sort the file by first column then by second column
sort -k2,2 -k1,1 MAG_clustering/all_fastANI_output.txt > MAG_clustering/all_fastANI_output_sorted_2.txt # Sort the file by second column then by first column
cut -f 1,2,3 MAG_clustering/all_fastANI_output_sorted_1.txt > temp.txt # Get the relevant columns of the first sorted file
cut -f 3 MAG_clustering/all_fastANI_output_sorted_2.txt > temp2.txt # Get the relevant columns of the second sorted file
paste temp.txt temp2.txt > MAG_clustering/all_fastANI_output_twoWay.txt # Combine the relevant columns
rm temp* # remove the temporary files
perl scripts/prepareANImatrix.pl MAG_clustering/all_fastANI_output_twoWay.txt MAG_clustering/genomePath_all.txt | sed 's/binning\/final_MAGs\/mealworm\///g' | sed 's/binning\/final_MAGs\/superworm\///g' | sed 's/\.fasta//g' > MAG_clustering/all_ANI_matrix.txt # make a two-way ANI matrix from the fastANI output

