# Make new directories and copy all bins to them
mkdir binning/dRep_input/
mkdir binning/dRep_output/
mkdir binning/dRep_input/mealworm/
mkdir binning/dRep_input/superworm/
mkdir binning/dRep_output/mealworm/
mkdir binning/dRep_output/superworm/
cp -r binning/DAS_Tool_output/mealworm*/*DASTool_bins binning/dRep_input/mealworm/
cp -r binning/DAS_Tool_output/superworm*/*DASTool_bins binning/dRep_input/superworm/

# Rename the mealworm files
rename 's/concoct/mealworm_Standard_1a_concoct/' binning/dRep_input/mealworm/mealworm_Standard_1a_DASTool_bins/*
rename 's/metabat/mealworm_Standard_1a_metabat/' binning/dRep_input/mealworm/mealworm_Standard_1a_DASTool_bins/*
rename 's/maxbin/mealworm_Standard_1a_maxbin/' binning/dRep_input/mealworm/mealworm_Standard_1a_DASTool_bins/*
rename 's/concoct/mealworm_Control_1a_concoct/' binning/dRep_input/mealworm/mealworm_Control_1a_DASTool_bins/*
rename 's/metabat/mealworm_Control_1a_metabat/' binning/dRep_input/mealworm/mealworm_Control_1a_DASTool_bins/*
rename 's/maxbin/mealworm_Control_1a_maxbin/' binning/dRep_input/mealworm/mealworm_Control_1a_DASTool_bins/*
rename 's/concoct/mealworm_PB_1a_concoct/' binning/dRep_input/mealworm/mealworm_PB_1a_DASTool_bins/*
rename 's/metabat/mealworm_PB_1a_metabat/' binning/dRep_input/mealworm/mealworm_PB_1a_DASTool_bins/*
rename 's/maxbin/mealworm_PB_1a_maxbin/' binning/dRep_input/mealworm/mealworm_PB_1a_DASTool_bins/*
rename 's/concoct/mealworm_Styrofoam_1a_concoct/' binning/dRep_input/mealworm/mealworm_Styrofoam_1a_DASTool_bins/*
rename 's/metabat/mealworm_Styrofoam_1a_metabat/' binning/dRep_input/mealworm/mealworm_Styrofoam_1a_DASTool_bins/*
rename 's/maxbin/mealworm_Styrofoam_1a_maxbin/' binning/dRep_input/mealworm/mealworm_Styrofoam_1a_DASTool_bins/*
rename 's/concoct/mealworm_Control_1b_concoct/' binning/dRep_input/mealworm/mealworm_Control_1b_DASTool_bins/*
rename 's/metabat/mealworm_Control_1b_metabat/' binning/dRep_input/mealworm/mealworm_Control_1b_DASTool_bins/*
rename 's/maxbin/mealworm_Control_1b_maxbin/' binning/dRep_input/mealworm/mealworm_Control_1b_DASTool_bins/*
rename 's/concoct/mealworm_PB_1b_concoct/' binning/dRep_input/mealworm/mealworm_PB_1b_DASTool_bins/*
rename 's/metabat/mealworm_PB_1b_metabat/' binning/dRep_input/mealworm/mealworm_PB_1b_DASTool_bins/*
rename 's/maxbin/mealworm_PB_1b_maxbin/' binning/dRep_input/mealworm/mealworm_PB_1b_DASTool_bins/*
rename 's/concoct/mealworm_Styrofoam_1b_concoct/' binning/dRep_input/mealworm/mealworm_Styrofoam_1b_DASTool_bins/*
rename 's/metabat/mealworm_Styrofoam_1b_metabat/' binning/dRep_input/mealworm/mealworm_Styrofoam_1b_DASTool_bins/*
rename 's/maxbin/mealworm_Styrofoam_1b_maxbin/' binning/dRep_input/mealworm/mealworm_Styrofoam_1b_DASTool_bins/*
rename 's/concoct/mealworm_Standard_2a_concoct/' binning/dRep_input/mealworm/mealworm_Standard_2a_DASTool_bins/*
rename 's/metabat/mealworm_Standard_2a_metabat/' binning/dRep_input/mealworm/mealworm_Standard_2a_DASTool_bins/*
rename 's/maxbin/mealworm_Standard_2a_maxbin/' binning/dRep_input/mealworm/mealworm_Standard_2a_DASTool_bins/*
rename 's/concoct/mealworm_Control_2a_concoct/' binning/dRep_input/mealworm/mealworm_Control_2a_DASTool_bins/*
rename 's/metabat/mealworm_Control_2a_metabat/' binning/dRep_input/mealworm/mealworm_Control_2a_DASTool_bins/*
rename 's/maxbin/mealworm_Control_2a_maxbin/' binning/dRep_input/mealworm/mealworm_Control_2a_DASTool_bins/*
rename 's/concoct/mealworm_PB_2a_concoct/' binning/dRep_input/mealworm/mealworm_PB_2a_DASTool_bins/*
rename 's/metabat/mealworm_PB_2a_metabat/' binning/dRep_input/mealworm/mealworm_PB_2a_DASTool_bins/*
rename 's/maxbin/mealworm_PB_2a_maxbin/' binning/dRep_input/mealworm/mealworm_PB_2a_DASTool_bins/*
rename 's/concoct/mealworm_Styrofoam_2a_concoct/' binning/dRep_input/mealworm/mealworm_Styrofoam_2a_DASTool_bins/*
rename 's/metabat/mealworm_Styrofoam_2a_metabat/' binning/dRep_input/mealworm/mealworm_Styrofoam_2a_DASTool_bins/*
rename 's/maxbin/mealworm_Styrofoam_2a_maxbin/' binning/dRep_input/mealworm/mealworm_Styrofoam_2a_DASTool_bins/*
rename 's/concoct/mealworm_Standard_2b_concoct/' binning/dRep_input/mealworm/mealworm_Standard_2b_DASTool_bins/*
rename 's/metabat/mealworm_Standard_2b_metabat/' binning/dRep_input/mealworm/mealworm_Standard_2b_DASTool_bins/*
rename 's/maxbin/mealworm_Standard_2b_maxbin/' binning/dRep_input/mealworm/mealworm_Standard_2b_DASTool_bins/*
rename 's/concoct/mealworm_Control_2b_concoct/' binning/dRep_input/mealworm/mealworm_Control_2b_DASTool_bins/*
rename 's/metabat/mealworm_Control_2b_metabat/' binning/dRep_input/mealworm/mealworm_Control_2b_DASTool_bins/*
rename 's/maxbin/mealworm_Control_2b_maxbin/' binning/dRep_input/mealworm/mealworm_Control_2b_DASTool_bins/*
rename 's/concoct/mealworm_PB_2b_concoct/' binning/dRep_input/mealworm/mealworm_PB_2b_DASTool_bins/*
rename 's/metabat/mealworm_PB_2b_metabat/' binning/dRep_input/mealworm/mealworm_PB_2b_DASTool_bins/*
rename 's/maxbin/mealworm_PB_2b_maxbin/' binning/dRep_input/mealworm/mealworm_PB_2b_DASTool_bins/*
rename 's/concoct/mealworm_Styrofoam_2b_concoct/' binning/dRep_input/mealworm/mealworm_Styrofoam_2b_DASTool_bins/*
rename 's/metabat/mealworm_Styrofoam_2b_metabat/' binning/dRep_input/mealworm/mealworm_Styrofoam_2b_DASTool_bins/*
rename 's/maxbin/mealworm_Styrofoam_2b_maxbin/' binning/dRep_input/mealworm/mealworm_Styrofoam_2b_DASTool_bins/*
rename 's/concoct/mealworm_Standard_3a_concoct/' binning/dRep_input/mealworm/mealworm_Standard_3a_DASTool_bins/*
rename 's/metabat/mealworm_Standard_3a_metabat/' binning/dRep_input/mealworm/mealworm_Standard_3a_DASTool_bins/*
rename 's/maxbin/mealworm_Standard_3a_maxbin/' binning/dRep_input/mealworm/mealworm_Standard_3a_DASTool_bins/*
rename 's/concoct/mealworm_Control_3a_concoct/' binning/dRep_input/mealworm/mealworm_Control_3a_DASTool_bins/*
rename 's/metabat/mealworm_Control_3a_metabat/' binning/dRep_input/mealworm/mealworm_Control_3a_DASTool_bins/*
rename 's/maxbin/mealworm_Control_3a_maxbin/' binning/dRep_input/mealworm/mealworm_Control_3a_DASTool_bins/*
rename 's/concoct/mealworm_PB_3a_concoct/' binning/dRep_input/mealworm/mealworm_PB_3a_DASTool_bins/*
rename 's/metabat/mealworm_PB_3a_metabat/' binning/dRep_input/mealworm/mealworm_PB_3a_DASTool_bins/*
rename 's/maxbin/mealworm_PB_3a_maxbin/' binning/dRep_input/mealworm/mealworm_PB_3a_DASTool_bins/*
rename 's/concoct/mealworm_Styrofoam_3a_concoct/' binning/dRep_input/mealworm/mealworm_Styrofoam_3a_DASTool_bins/*
rename 's/metabat/mealworm_Styrofoam_3a_metabat/' binning/dRep_input/mealworm/mealworm_Styrofoam_3a_DASTool_bins/*
rename 's/maxbin/mealworm_Styrofoam_3a_maxbin/' binning/dRep_input/mealworm/mealworm_Styrofoam_3a_DASTool_bins/*
rename 's/concoct/mealworm_Control_3b_concoct/' binning/dRep_input/mealworm/mealworm_Control_3b_DASTool_bins/*
rename 's/metabat/mealworm_Control_3b_metabat/' binning/dRep_input/mealworm/mealworm_Control_3b_DASTool_bins/*
rename 's/maxbin/mealworm_Control_3b_maxbin/' binning/dRep_input/mealworm/mealworm_Control_3b_DASTool_bins/*
rename 's/concoct/mealworm_PB_3b_concoct/' binning/dRep_input/mealworm/mealworm_PB_3b_DASTool_bins/*
rename 's/metabat/mealworm_PB_3b_metabat/' binning/dRep_input/mealworm/mealworm_PB_3b_DASTool_bins/*
rename 's/maxbin/mealworm_PB_3b_maxbin/' binning/dRep_input/mealworm/mealworm_PB_3b_DASTool_bins/*
rename 's/concoct/mealworm_Styrofoam_3b_concoct/' binning/dRep_input/mealworm/mealworm_Styrofoam_3b_DASTool_bins/*
rename 's/metabat/mealworm_Styrofoam_3b_metabat/' binning/dRep_input/mealworm/mealworm_Styrofoam_3b_DASTool_bins/*
rename 's/maxbin/mealworm_Styrofoam_3b_maxbin/' binning/dRep_input/mealworm/mealworm_Styrofoam_3b_DASTool_bins/*

# Rename the superworm files
rename 's/concoct/superworm_Standard_1a_concoct/' binning/dRep_input/superworm/superworm_Standard_1a_DASTool_bins/*
rename 's/metabat/superworm_Standard_1a_metabat/' binning/dRep_input/superworm/superworm_Standard_1a_DASTool_bins/*
rename 's/maxbin/superworm_Standard_1a_maxbin/' binning/dRep_input/superworm/superworm_Standard_1a_DASTool_bins/*
rename 's/concoct/superworm_Control_1a_concoct/' binning/dRep_input/superworm/superworm_Control_1a_DASTool_bins/*
rename 's/metabat/superworm_Control_1a_metabat/' binning/dRep_input/superworm/superworm_Control_1a_DASTool_bins/*
rename 's/maxbin/superworm_Control_1a_maxbin/' binning/dRep_input/superworm/superworm_Control_1a_DASTool_bins/*
rename 's/concoct/superworm_PB_1a_concoct/' binning/dRep_input/superworm/superworm_PB_1a_DASTool_bins/*
rename 's/metabat/superworm_PB_1a_metabat/' binning/dRep_input/superworm/superworm_PB_1a_DASTool_bins/*
rename 's/maxbin/superworm_PB_1a_maxbin/' binning/dRep_input/superworm/superworm_PB_1a_DASTool_bins/*
rename 's/concoct/superworm_Styrofoam_1a_concoct/' binning/dRep_input/superworm/superworm_Styrofoam_1a_DASTool_bins/*
rename 's/metabat/superworm_Styrofoam_1a_metabat/' binning/dRep_input/superworm/superworm_Styrofoam_1a_DASTool_bins/*
rename 's/maxbin/superworm_Styrofoam_1a_maxbin/' binning/dRep_input/superworm/superworm_Styrofoam_1a_DASTool_bins/*
rename 's/concoct/superworm_Standard_1b_concoct/' binning/dRep_input/superworm/superworm_Standard_1b_DASTool_bins/*
rename 's/metabat/superworm_Standard_1b_metabat/' binning/dRep_input/superworm/superworm_Standard_1b_DASTool_bins/*
rename 's/maxbin/superworm_Standard_1b_maxbin/' binning/dRep_input/superworm/superworm_Standard_1b_DASTool_bins/*
rename 's/concoct/superworm_Control_1b_concoct/' binning/dRep_input/superworm/superworm_Control_1b_DASTool_bins/*
rename 's/metabat/superworm_Control_1b_metabat/' binning/dRep_input/superworm/superworm_Control_1b_DASTool_bins/*
rename 's/maxbin/superworm_Control_1b_maxbin/' binning/dRep_input/superworm/superworm_Control_1b_DASTool_bins/*
rename 's/concoct/superworm_PB_1b_concoct/' binning/dRep_input/superworm/superworm_PB_1b_DASTool_bins/*
rename 's/metabat/superworm_PB_1b_metabat/' binning/dRep_input/superworm/superworm_PB_1b_DASTool_bins/*
rename 's/maxbin/superworm_PB_1b_maxbin/' binning/dRep_input/superworm/superworm_PB_1b_DASTool_bins/*
rename 's/concoct/superworm_Styrofoam_1b_concoct/' binning/dRep_input/superworm/superworm_Styrofoam_1b_DASTool_bins/*
rename 's/metabat/superworm_Styrofoam_1b_metabat/' binning/dRep_input/superworm/superworm_Styrofoam_1b_DASTool_bins/*
rename 's/maxbin/superworm_Styrofoam_1b_maxbin/' binning/dRep_input/superworm/superworm_Styrofoam_1b_DASTool_bins/*
rename 's/concoct/superworm_Standard_2a_concoct/' binning/dRep_input/superworm/superworm_Standard_2a_DASTool_bins/*
rename 's/metabat/superworm_Standard_2a_metabat/' binning/dRep_input/superworm/superworm_Standard_2a_DASTool_bins/*
rename 's/maxbin/superworm_Standard_2a_maxbin/' binning/dRep_input/superworm/superworm_Standard_2a_DASTool_bins/*
rename 's/concoct/superworm_Control_2a_concoct/' binning/dRep_input/superworm/superworm_Control_2a_DASTool_bins/*
rename 's/metabat/superworm_Control_2a_metabat/' binning/dRep_input/superworm/superworm_Control_2a_DASTool_bins/*
rename 's/maxbin/superworm_Control_2a_maxbin/' binning/dRep_input/superworm/superworm_Control_2a_DASTool_bins/*
rename 's/concoct/superworm_PB_2a_concoct/' binning/dRep_input/superworm/superworm_PB_2a_DASTool_bins/*
rename 's/metabat/superworm_PB_2a_metabat/' binning/dRep_input/superworm/superworm_PB_2a_DASTool_bins/*
rename 's/maxbin/superworm_PB_2a_maxbin/' binning/dRep_input/superworm/superworm_PB_2a_DASTool_bins/*
rename 's/concoct/superworm_Styrofoam_2a_concoct/' binning/dRep_input/superworm/superworm_Styrofoam_2a_DASTool_bins/*
rename 's/metabat/superworm_Styrofoam_2a_metabat/' binning/dRep_input/superworm/superworm_Styrofoam_2a_DASTool_bins/*
rename 's/maxbin/superworm_Styrofoam_2a_maxbin/' binning/dRep_input/superworm/superworm_Styrofoam_2a_DASTool_bins/*
rename 's/concoct/superworm_Standard_2b_concoct/' binning/dRep_input/superworm/superworm_Standard_2b_DASTool_bins/*
rename 's/metabat/superworm_Standard_2b_metabat/' binning/dRep_input/superworm/superworm_Standard_2b_DASTool_bins/*
rename 's/maxbin/superworm_Standard_2b_maxbin/' binning/dRep_input/superworm/superworm_Standard_2b_DASTool_bins/*
rename 's/concoct/superworm_Control_2b_concoct/' binning/dRep_input/superworm/superworm_Control_2b_DASTool_bins/*
rename 's/metabat/superworm_Control_2b_metabat/' binning/dRep_input/superworm/superworm_Control_2b_DASTool_bins/*
rename 's/maxbin/superworm_Control_2b_maxbin/' binning/dRep_input/superworm/superworm_Control_2b_DASTool_bins/*
rename 's/concoct/superworm_PB_2b_concoct/' binning/dRep_input/superworm/superworm_PB_2b_DASTool_bins/*
rename 's/metabat/superworm_PB_2b_metabat/' binning/dRep_input/superworm/superworm_PB_2b_DASTool_bins/*
rename 's/maxbin/superworm_PB_2b_maxbin/' binning/dRep_input/superworm/superworm_PB_2b_DASTool_bins/*
rename 's/concoct/superworm_Styrofoam_2b_concoct/' binning/dRep_input/superworm/superworm_Styrofoam_2b_DASTool_bins/*
rename 's/metabat/superworm_Styrofoam_2b_metabat/' binning/dRep_input/superworm/superworm_Styrofoam_2b_DASTool_bins/*
rename 's/maxbin/superworm_Styrofoam_2b_maxbin/' binning/dRep_input/superworm/superworm_Styrofoam_2b_DASTool_bins/*
rename 's/concoct/superworm_Standard_3a_concoct/' binning/dRep_input/superworm/superworm_Standard_3a_DASTool_bins/*
rename 's/metabat/superworm_Standard_3a_metabat/' binning/dRep_input/superworm/superworm_Standard_3a_DASTool_bins/*
rename 's/maxbin/superworm_Standard_3a_maxbin/' binning/dRep_input/superworm/superworm_Standard_3a_DASTool_bins/*
rename 's/concoct/superworm_Control_3a_concoct/' binning/dRep_input/superworm/superworm_Control_3a_DASTool_bins/*
rename 's/metabat/superworm_Control_3a_metabat/' binning/dRep_input/superworm/superworm_Control_3a_DASTool_bins/*
rename 's/maxbin/superworm_Control_3a_maxbin/' binning/dRep_input/superworm/superworm_Control_3a_DASTool_bins/*
rename 's/concoct/superworm_PB_3a_concoct/' binning/dRep_input/superworm/superworm_PB_3a_DASTool_bins/*
rename 's/metabat/superworm_PB_3a_metabat/' binning/dRep_input/superworm/superworm_PB_3a_DASTool_bins/*
rename 's/maxbin/superworm_PB_3a_maxbin/' binning/dRep_input/superworm/superworm_PB_3a_DASTool_bins/*
rename 's/concoct/superworm_Styrofoam_3a_concoct/' binning/dRep_input/superworm/superworm_Styrofoam_3a_DASTool_bins/*
rename 's/metabat/superworm_Styrofoam_3a_metabat/' binning/dRep_input/superworm/superworm_Styrofoam_3a_DASTool_bins/*
rename 's/maxbin/superworm_Styrofoam_3a_maxbin/' binning/dRep_input/superworm/superworm_Styrofoam_3a_DASTool_bins/*
rename 's/concoct/superworm_Standard_3b_concoct/' binning/dRep_input/superworm/superworm_Standard_3b_DASTool_bins/*
rename 's/metabat/superworm_Standard_3b_metabat/' binning/dRep_input/superworm/superworm_Standard_3b_DASTool_bins/*
rename 's/maxbin/superworm_Standard_3b_maxbin/' binning/dRep_input/superworm/superworm_Standard_3b_DASTool_bins/*
rename 's/concoct/superworm_Control_3b_concoct/' binning/dRep_input/superworm/superworm_Control_3b_DASTool_bins/*
rename 's/metabat/superworm_Control_3b_metabat/' binning/dRep_input/superworm/superworm_Control_3b_DASTool_bins/*
rename 's/maxbin/superworm_Control_3b_maxbin/' binning/dRep_input/superworm/superworm_Control_3b_DASTool_bins/*
rename 's/concoct/superworm_PB_3b_concoct/' binning/dRep_input/superworm/superworm_PB_3b_DASTool_bins/*
rename 's/metabat/superworm_PB_3b_metabat/' binning/dRep_input/superworm/superworm_PB_3b_DASTool_bins/*
rename 's/maxbin/superworm_PB_3b_maxbin/' binning/dRep_input/superworm/superworm_PB_3b_DASTool_bins/*
rename 's/concoct/superworm_Styrofoam_3b_concoct/' binning/dRep_input/superworm/superworm_Styrofoam_3b_DASTool_bins/*
rename 's/metabat/superworm_Styrofoam_3b_metabat/' binning/dRep_input/superworm/superworm_Styrofoam_3b_DASTool_bins/*
rename 's/maxbin/superworm_Styrofoam_3b_maxbin/' binning/dRep_input/superworm/superworm_Styrofoam_3b_DASTool_bins/*

# Run dRep
dRep dereplicate binning/dRep_output/mealworm/ -g binning/dRep_input/mealworm/*/*.fa -p 16 -comp 70 -con 10
dRep dereplicate binning/dRep_output/superworm/ -g binning/dRep_input/superworm/*/*.fa -p 16 -comp 70 -con 10

# Rename MAGs
mkdir binning/final_MAGs/
mkdir binning/final_MAGs/mealworm/
mkdir binning/final_MAGs/superworm/
perl scripts/rename_MAGs.pl Mealworm_MAG binning/dRep_output/mealworm/dereplicated_genomes/ binning/final_MAGs/mealworm/ 0 > binning/final_MAGs/mealworm/log.txt
perl scripts/rename_MAGs.pl Superworm_MAG binning/dRep_output/superworm/dereplicated_genomes/ binning/final_MAGs/superworm/ 0 > binning/final_MAGs/superworm/log.txt
