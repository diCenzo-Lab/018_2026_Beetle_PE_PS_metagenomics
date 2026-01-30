# Get the eukaryotic MAGs for mealworm
mkdir binning/dRep_output/mealworm_euk/
dRep dereplicate binning/dRep_output/mealworm_euk/ -g binning/dRep_input/mealworm/*/*.fa -p 16 -comp 0 -con 0
mkdir binning/dRep_output/mealworm_euk/busco/
mkdir binning/dRep_output/mealworm_euk/final/
cd binning/dRep_output/mealworm_euk/dereplicated_genomes/
for i in *.fa
do
    busco -i "${i}" -o ../busco/"${i}" -l eukaryota_odb10 -c 16 -m genome -f
    perl ../../../../scripts/get_euk.pl "${i}"
done;
cd ../../../../
perl scripts/rename_MAGs.pl Mealworm_eMAG binning/dRep_output/mealworm_euk/final/ binning/final_MAGs/mealworm/ 0 > binning/final_MAGs/mealworm/log_euk.txt

# Get the eukaryotic MAGs for superworm
mkdir binning/dRep_output/superworm_euk/
dRep dereplicate binning/dRep_output/superworm_euk/ -g binning/dRep_input/superworm/*/*.fa -p 16 -comp 0 -con 0
mkdir binning/dRep_output/superworm_euk/busco/
mkdir binning/dRep_output/superworm_euk/final/
cd binning/dRep_output/superworm_euk/dereplicated_genomes/
for i in *.fa
do
    busco -i "${i}" -o ../busco/"${i}" -l eukaryota_odb10 -c 16 -m genome -f
    perl ../../../../scripts/get_euk.pl "${i}"
done;
cd ../../../../
perl scripts/rename_MAGs.pl Superworm_eMAG binning/dRep_output/superworm_euk/final/ binning/final_MAGs/superworm/ 0 > binning/final_MAGs/superworm/log_euk.txt
