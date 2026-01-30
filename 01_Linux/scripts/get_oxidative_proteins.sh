# Create the HMMs
hmmbuild secretome/oxidative_proteins/input_data/PF00394.hmm secretome/oxidative_proteins/input_data/PF00394_seed.txt # Create HMM
hmmbuild secretome/oxidative_proteins/input_data/PF07731.hmm secretome/oxidative_proteins/input_data/PF07731_seed.txt # Create HMM
hmmbuild secretome/oxidative_proteins/input_data/PF07732.hmm secretome/oxidative_proteins/input_data/PF07732_seed.txt # Create HMM
hmmbuild secretome/oxidative_proteins/input_data/PF00255.hmm secretome/oxidative_proteins/input_data/PF00255_seed.txt # Create HMM
hmmbuild secretome/oxidative_proteins/input_data/PF04261.hmm secretome/oxidative_proteins/input_data/PF04261_seed.txt # Create HMM
hmmbuild secretome/oxidative_proteins/input_data/PF21105.hmm secretome/oxidative_proteins/input_data/PF21105_seed.txt # Create HMM
hmmbuild secretome/oxidative_proteins/input_data/PF20628.hmm secretome/oxidative_proteins/input_data/PF20628_seed.txt # Create HMM
hmmbuild secretome/oxidative_proteins/input_data/PF13532.hmm secretome/oxidative_proteins/input_data/PF13532_seed.txt # Create HMM
hmmbuild secretome/oxidative_proteins/input_data/PF00487.hmm secretome/oxidative_proteins/input_data/PF00487_seed.txt # Create HMM

# Identify and extract the putative orthologs
mkdir secretome/oxidative_proteins/hmmsearch_output/ # Make directory to hold the output data
cat secretome/oxidative_proteins/input_data/reported_enzyme_sequences.faa secretome/full_secretome.fasta > secretome/oxidative_proteins/input_data/combined_proteins.fasta # Combine the secreted proteins with the proteins previously implicated in plastics modification
hmmsearch -o secretome/oxidative_proteins/hmmsearch_output/PF00394_hmmsearch.txt --noali --cpu 4 -E 0.01 secretome/oxidative_proteins/input_data/PF00394.hmm secretome/oxidative_proteins/input_data/combined_proteins.fasta # Run HMMSEARH with one of the HMMs
hmmsearch -o secretome/oxidative_proteins/hmmsearch_output/PF07731_hmmsearch.txt --noali --cpu 4 -E 0.01 secretome/oxidative_proteins/input_data/PF07731.hmm secretome/oxidative_proteins/input_data/combined_proteins.fasta # Run HMMSEARH with one of the HMMs
hmmsearch -o secretome/oxidative_proteins/hmmsearch_output/PF07732_hmmsearch.txt --noali --cpu 4 -E 0.01 secretome/oxidative_proteins/input_data/PF07732.hmm secretome/oxidative_proteins/input_data/combined_proteins.fasta # Run HMMSEARH with one of the HMMs
hmmsearch -o secretome/oxidative_proteins/hmmsearch_output/PF00255_hmmsearch.txt --noali --cpu 4 -E 0.01 secretome/oxidative_proteins/input_data/PF00255.hmm secretome/oxidative_proteins/input_data/combined_proteins.fasta # Run HMMSEARH with one of the HMMs
hmmsearch -o secretome/oxidative_proteins/hmmsearch_output/PF04261_hmmsearch.txt --noali --cpu 4 -E 0.01 secretome/oxidative_proteins/input_data/PF04261.hmm secretome/oxidative_proteins/input_data/combined_proteins.fasta # Run HMMSEARH with one of the HMMs
hmmsearch -o secretome/oxidative_proteins/hmmsearch_output/PF21105_hmmsearch.txt --noali --cpu 4 -E 0.01 secretome/oxidative_proteins/input_data/PF21105.hmm secretome/oxidative_proteins/input_data/combined_proteins.fasta # Run HMMSEARH with one of the HMMs
hmmsearch -o secretome/oxidative_proteins/hmmsearch_output/PF20628_hmmsearch.txt --noali --cpu 4 -E 0.01 secretome/oxidative_proteins/input_data/PF20628.hmm secretome/oxidative_proteins/input_data/combined_proteins.fasta # Run HMMSEARH with one of the HMMs
hmmsearch -o secretome/oxidative_proteins/hmmsearch_output/PF13532_hmmsearch.txt --noali --cpu 4 -E 0.01 secretome/oxidative_proteins/input_data/PF13532.hmm secretome/oxidative_proteins/input_data/combined_proteins.fasta # Run HMMSEARH with one of the HMMs
hmmsearch -o secretome/oxidative_proteins/hmmsearch_output/PF00487_hmmsearch.txt --noali --cpu 4 -E 0.01 secretome/oxidative_proteins/input_data/PF00487.hmm secretome/oxidative_proteins/input_data/combined_proteins.fasta # Run HMMSEARH with one of the HMMs
perl scripts/modifyFasta.pl secretome/oxidative_proteins/input_data/combined_proteins.fasta > secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta # Modify the format of the file for easy extraction
perl scripts/parseHMMsearch.pl secretome/oxidative_proteins/hmmsearch_output/PF00394_hmmsearch.txt > secretome/oxidative_proteins/hmmsearch_output/PF00394_hmmsearch.parsed.txt # Parse one of the hmmsearch output files
perl scripts/parseHMMsearch.pl secretome/oxidative_proteins/hmmsearch_output/PF07731_hmmsearch.txt > secretome/oxidative_proteins/hmmsearch_output/PF07731_hmmsearch.parsed.txt # Parse one of the hmmsearch output files
perl scripts/parseHMMsearch.pl secretome/oxidative_proteins/hmmsearch_output/PF07732_hmmsearch.txt > secretome/oxidative_proteins/hmmsearch_output/PF07732_hmmsearch.parsed.txt # Parse one of the hmmsearch output files
perl scripts/parseHMMsearch.pl secretome/oxidative_proteins/hmmsearch_output/PF00255_hmmsearch.txt > secretome/oxidative_proteins/hmmsearch_output/PF00255_hmmsearch.parsed.txt # Parse one of the hmmsearch output files
perl scripts/parseHMMsearch.pl secretome/oxidative_proteins/hmmsearch_output/PF04261_hmmsearch.txt > secretome/oxidative_proteins/hmmsearch_output/PF04261_hmmsearch.parsed.txt # Parse one of the hmmsearch output files
perl scripts/parseHMMsearch.pl secretome/oxidative_proteins/hmmsearch_output/PF21105_hmmsearch.txt > secretome/oxidative_proteins/hmmsearch_output/PF21105_hmmsearch.parsed.txt # Parse one of the hmmsearch output files
perl scripts/parseHMMsearch.pl secretome/oxidative_proteins/hmmsearch_output/PF20628_hmmsearch.txt > secretome/oxidative_proteins/hmmsearch_output/PF20628_hmmsearch.parsed.txt # Parse one of the hmmsearch output files
perl scripts/parseHMMsearch.pl secretome/oxidative_proteins/hmmsearch_output/PF13532_hmmsearch.txt > secretome/oxidative_proteins/hmmsearch_output/PF13532_hmmsearch.parsed.txt # Parse one of the hmmsearch output files
perl scripts/parseHMMsearch.pl secretome/oxidative_proteins/hmmsearch_output/PF00487_hmmsearch.txt > secretome/oxidative_proteins/hmmsearch_output/PF00487_hmmsearch.parsed.txt # Parse one of the hmmsearch output files
perl scripts/extractHMMsearchHits.pl secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta secretome/oxidative_proteins/hmmsearch_output/PF00394_hmmsearch.parsed.txt > secretome/oxidative_proteins/hmmsearch_output/PF00394.fasta # Extract the proteins corresponding to the hits
perl scripts/extractHMMsearchHits.pl secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta secretome/oxidative_proteins/hmmsearch_output/PF07731_hmmsearch.parsed.txt > secretome/oxidative_proteins/hmmsearch_output/PF07731.fasta # Extract the proteins corresponding to the hits
perl scripts/extractHMMsearchHits.pl secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta secretome/oxidative_proteins/hmmsearch_output/PF07732_hmmsearch.parsed.txt > secretome/oxidative_proteins/hmmsearch_output/PF07732.fasta # Extract the proteins corresponding to the hits
cat secretome/oxidative_proteins/hmmsearch_output/PF00394.fasta secretome/oxidative_proteins/hmmsearch_output/PF07731.fasta secretome/oxidative_proteins/hmmsearch_output/PF07732.fasta > secretome/oxidative_proteins/hmmsearch_output/Cu_oxidases.fasta # Merge potential laccases
perl scripts/modifyFasta.pl secretome/oxidative_proteins/hmmsearch_output/Cu_oxidases.fasta | sort -u | sed 's/\t/\n/' > secretome/oxidative_proteins/hmmsearch_output/Cu_oxidases.unique.fasta # Remove duplicates
perl scripts/extractHMMsearchHits.pl secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta secretome/oxidative_proteins/hmmsearch_output/PF00255_hmmsearch.parsed.txt > secretome/oxidative_proteins/hmmsearch_output/PF00255.fasta # Extract the proteins corresponding to the hits
perl scripts/extractHMMsearchHits.pl secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta secretome/oxidative_proteins/hmmsearch_output/PF04261_hmmsearch.parsed.txt > secretome/oxidative_proteins/hmmsearch_output/PF04261.fasta # Extract the proteins corresponding to the hits
perl scripts/extractHMMsearchHits.pl secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta secretome/oxidative_proteins/hmmsearch_output/PF21105_hmmsearch.parsed.txt > secretome/oxidative_proteins/hmmsearch_output/PF21105.fasta # Extract the proteins corresponding to the hits
perl scripts/extractHMMsearchHits.pl secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta secretome/oxidative_proteins/hmmsearch_output/PF20628_hmmsearch.parsed.txt > secretome/oxidative_proteins/hmmsearch_output/PF20628.fasta # Extract the proteins corresponding to the hits
cat secretome/oxidative_proteins/hmmsearch_output/PF04261.fasta secretome/oxidative_proteins/hmmsearch_output/PF21105.fasta secretome/oxidative_proteins/hmmsearch_output/PF20628.fasta > secretome/oxidative_proteins/hmmsearch_output/DyPs.fasta # Merge potential laccases
perl scripts/modifyFasta.pl secretome/oxidative_proteins/hmmsearch_output/DyPs.fasta | sort -u | sed 's/\t/\n/' > secretome/oxidative_proteins/hmmsearch_output/DyPs.unique.fasta # Remove duplicates
perl scripts/extractHMMsearchHits.pl secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta secretome/oxidative_proteins/hmmsearch_output/PF13532_hmmsearch.parsed.txt > secretome/oxidative_proteins/hmmsearch_output/PF13532.fasta # Extract the proteins corresponding to the hits
perl scripts/extractHMMsearchHits.pl secretome/oxidative_proteins/input_data/combined_proteins.modified.fasta secretome/oxidative_proteins/hmmsearch_output/PF00487_hmmsearch.parsed.txt > secretome/oxidative_proteins/hmmsearch_output/PF00487.fasta # Extract the proteins corresponding to the hits

# Get Pfam database
mkdir secretome/oxidative_proteins/Pfam_database/ # Make directory to hold the files
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz # get the Pfam v 38.0 HMM files
hmmconvert Pfam-A.hmm > Pfam-A_converted.hmm # convert the database to the necessary format
hmmpress Pfam-A_converted.hmm # prepare files for hmmscan searches
mv Pfam-A* secretome/oxidative_proteins/Pfam_database/ # Move files to the correct directory

# Run HMMscan and remove hits that are more closely related to other Pfam families
mkdir secretome/oxidative_proteins/hmmscan_output/ # Make directory to hold the files
hmmscan --cpu 16 secretome/oxidative_proteins/Pfam_database/Pfam-A_converted.hmm secretome/oxidative_proteins/hmmsearch_output/Cu_oxidases.unique.fasta > secretome/oxidative_proteins/hmmscan_output/Cu_oxidases.txt # Run HMMscan against the Pfam database
perl scripts/parseHMMscan.pl secretome/oxidative_proteins/hmmscan_output/Cu_oxidases.txt > secretome/oxidative_proteins/hmmscan_output/Cu_oxidases.parsed.csv
perl scripts/HMMscanTopHit.pl secretome/oxidative_proteins/hmmscan_output/Cu_oxidases.parsed.csv > secretome/oxidative_proteins/hmmscan_output/Cu_oxidases.top.csv
grep 'Cu-oxidase' secretome/oxidative_proteins/hmmscan_output/Cu_oxidases.top.csv | cut -f1 -d',' | pullseq -N -i secretome/oxidative_proteins/hmmsearch_output/Cu_oxidases.unique.fasta > secretome/oxidative_proteins/hmmscan_output/Cu_oxidases.fasta
hmmscan --cpu 16 secretome/oxidative_proteins/Pfam_database/Pfam-A_converted.hmm secretome/oxidative_proteins/hmmsearch_output/PF00255.fasta > secretome/oxidative_proteins/hmmscan_output/PF00255.txt # Run HMMscan against the Pfam database
perl scripts/parseHMMscan.pl secretome/oxidative_proteins/hmmscan_output/PF00255.txt > secretome/oxidative_proteins/hmmscan_output/PF00255.parsed.csv
perl scripts/HMMscanTopHit.pl secretome/oxidative_proteins/hmmscan_output/PF00255.parsed.csv > secretome/oxidative_proteins/hmmscan_output/PF00255.top.csv
grep 'GSHPx' secretome/oxidative_proteins/hmmscan_output/PF00255.top.csv | cut -f1 -d',' | pullseq -N -i secretome/oxidative_proteins/hmmsearch_output/PF00255.fasta > secretome/oxidative_proteins/hmmscan_output/PF00255.fasta
hmmscan --cpu 16 secretome/oxidative_proteins/Pfam_database/Pfam-A_converted.hmm secretome/oxidative_proteins/hmmsearch_output/DyPs.unique.fasta > secretome/oxidative_proteins/hmmscan_output/DyPs.txt # Run HMMscan against the Pfam database
perl scripts/parseHMMscan.pl secretome/oxidative_proteins/hmmscan_output/DyPs.txt > secretome/oxidative_proteins/hmmscan_output/DyPs.parsed.csv
perl scripts/HMMscanTopHit.pl secretome/oxidative_proteins/hmmscan_output/DyPs.parsed.csv > secretome/oxidative_proteins/hmmscan_output/DyPs.top.csv
grep 'Dy' secretome/oxidative_proteins/hmmscan_output/DyPs.top.csv | cut -f1 -d',' | pullseq -N -i secretome/oxidative_proteins/hmmsearch_output/DyPs.unique.fasta > secretome/oxidative_proteins/hmmscan_output/DyPs.fasta
hmmscan --cpu 16 secretome/oxidative_proteins/Pfam_database/Pfam-A_converted.hmm secretome/oxidative_proteins/hmmsearch_output/PF13532.fasta > secretome/oxidative_proteins/hmmscan_output/PF13532.txt # Run HMMscan against the Pfam database
perl scripts/parseHMMscan.pl secretome/oxidative_proteins/hmmscan_output/PF13532.txt > secretome/oxidative_proteins/hmmscan_output/PF13532.parsed.csv
perl scripts/HMMscanTopHit.pl secretome/oxidative_proteins/hmmscan_output/PF13532.parsed.csv > secretome/oxidative_proteins/hmmscan_output/PF13532.top.csv
grep '2OG-FeII_Oxy_2' secretome/oxidative_proteins/hmmscan_output/PF13532.top.csv | cut -f1 -d',' | pullseq -N -i secretome/oxidative_proteins/hmmsearch_output/PF13532.fasta > secretome/oxidative_proteins/hmmscan_output/PF13532.fasta
hmmscan --cpu 16 secretome/oxidative_proteins/Pfam_database/Pfam-A_converted.hmm secretome/oxidative_proteins/hmmsearch_output/PF00487.fasta > secretome/oxidative_proteins/hmmscan_output/PF00487.txt # Run HMMscan against the Pfam database
perl scripts/parseHMMscan.pl secretome/oxidative_proteins/hmmscan_output/PF00487.txt > secretome/oxidative_proteins/hmmscan_output/PF00487.parsed.csv
perl scripts/HMMscanTopHit.pl secretome/oxidative_proteins/hmmscan_output/PF00487.parsed.csv > secretome/oxidative_proteins/hmmscan_output/PF00487.top.csv
grep 'FA_desaturase' secretome/oxidative_proteins/hmmscan_output/PF00487.top.csv | cut -f1 -d',' | pullseq -N -i secretome/oxidative_proteins/hmmsearch_output/PF00487.fasta > secretome/oxidative_proteins/hmmscan_output/PF00487.fasta
