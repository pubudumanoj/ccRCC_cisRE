File running order (needs modifications in the code)

1. H3K27ac_aggregate_heatmap.promo_gbody_3FR.pbs
2. H3K27ac.job.diff.norm.col.mean.pbs
3. supp_fig_1.Rmd

#################################################################################################
#################################################################################################
#################################################################################################
Description about H3K27ac_aggregate_heatmap.promo_gbody_3FR.pbs

you have to create script for H3K4me3 and H3K4me1 using this script
You only need to do is change the histone mark in the script and create tables for file names

then run the script and create the files



#################################################################################################
#################################################################################################
#################################################################################################
Description about H3K27ac.job.diff.norm.col.mean.pbs

you have to create script for H3K4me3 and H3K4me1 using this script
You only need to do is change the histone mark in the script and create tables for file names

then run the script and create the files



#################################################################################################
#################################################################################################
#################################################################################################


Description about supp_fig_1.Rmd

This script creates graphs in supplementary figure 1

supplementary figure 1A
inputs: 1
1. hg19ToHg38.over.chain (downloaded from encode)

It contains chr and it size (2 columns) 

Since reference files are mapped in to hg19 we have to convert those coordinates to HG38 to align with our samples

2. Downloaded reference histone mark files

We only need 1st three columns

Chr, start and end position

3. our peak files stored in makepeak folder

These chip peak files were made from original peak files obtained from MACS2. Just select the chromosome coordinates from the files and create new files.

Three columns

Chr, start and end coordinates

#################################################################################################
#################################################################################################
#################################################################################################

supplementary figure 1B

Inputs: chip pileups (bedgraphs)

They are in chip_pileups folder

first we have to normalized them
Used the code in the file and create normalized bedgraph files for all the histone modifications and all the samples

e.g: R354_N_1_Kidney_ChIP_H3K27ac_treat_pileup.bdg

columns
chr, start, end, value

ChIP pileup values are calculated in this way: 1) Sequencing reads are extended into N bps fragments according to fragment size prediction step; 2) For each position, MACS2 computes how many fragments can be found; 3) Only if ChIP has more reads than control, these values will be multiplied by control_depth/ChIP_depth (default behavior). The 'pileup' value you saw in  file is the value at peak summit, which means how many fragments you can find covering the peak summit. 

chr8    145073482       145074204       0.00000
chr8    145074204       145074388       0.60602
chr8    145074388       145074467       0.00000
chr8    145074467       145074615       0.60602
chr8    145074615       145074635       1.21204
chr8    145074635       145074651       1.81807










