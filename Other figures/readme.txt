Before running the scripts in this folder you should run all the scripts in sup fig 1 and fig 1 folders and create all the files

After that start from step 1 file and then run the step 2 file

They are related to most of the figures

main input file for step 1 is 
"resulted enhancers from binomial test_p_value_combined_promoter5kb_new.csv" created from previous steps

It has many columns (based on GRCh37)
This is the description of each column 

	Peak = ID of the enhancer
	chr	=  chromosome of where the enhancer locatedstart
	start	= start coordinates of 1kb enhancer
	end	= end coordinates of 1kb enhancer
	Width	= width of the enhancer (1kb)
	Methylation enrichment on 10kb region	= whether enriched in hypo- or hypomethylated CpGs (WGBS)
	WGBS hypomethylated CpGs count on the enhancer 10kb region
	WGBS hypermethylated CpGs count on the enhancer 10kb region
	Total WGBS CpGs on the enhancer 10kb region
	P value of 10 kb region
	FDR of 10kb region
	RCC4_450k_array hypermethylated CpGs count on the enhancer 10kb region
	RCC4_450k_array hypomethylated CpGs count on the enhancer 10kb region
	RCC4_450k_array total CpGs count on the enhancer 10kb region
	Number of Overlaps of 10kb region with reference-lost enhancers in tumor(N)
	Number of overlaps of 10kb region with reference-gained enhancers in tumor(T)
	Overlaps of 10kb region with query-enhancers from tumor samples
	Overlaps of 10kb region with query-enhancers from normal samples
	closest_gene from the centre of enhancer
	Gene name
	distance_to_the_gene_TSS_from_centre_of_enhancer
	log2FoldChange_VHL_vs_EV_H	= gene expression change when compared VHL+ve relative to mock cells in hypoxic condition
	padj_VHL_vs_EV_H	= FDR of log2foldchange when compared VHL+ve relative to mock cells in hypoxic condition
	log2FoldChange_VHL_vs_EV_N	= gene expression change when compared VHL+ve relative to mock cells in normoxic condition
	padj_VHL_vs_EV_N	= FDR of log2foldchange when compared VHL+ve relative to mock cells in normoxic condition 
	log2FoldChange_patients_T_vs_N	= gene expression change when compared tumor relative to normal cells (29 patients)
	padj_patients_T_vs_N	= DR of log2foldchange when compared tumor relative to normal cells (29 patients)
	Methylation enrichment on 1kb region = whether enriched in hypo- or hypomethylated CpGs (WGBS)
	WGBS hypomethylated CpGs count on the enhancer 1kb region
	WGBS hypermethylated CpGs count on the enhancer 1kb region
	Total WGBS CpGs on the enhancer 1kb region
	P value of 1kb region
	FDR of 1kb region
	RCC4_450k_array hypermethylated CpGs count on the enhancer 1kb region
	RCC4_450k_array hypomethylated CpGs count on the enhancer 1kb region
	RCC4_450k_array total CpGs count on the enhancer 1kb region
	Number of Overlaps of 1kb region with reference-lost enhancers in tumor(N)
	Number of overlaps of 1kb region with reference-gained enhancers in tumor(T)
	Overlaps of 1kb region with query-enhancers from tumor samples
	Overlaps of 1kb region with query-enhancers from normal samples
	log2FoldChange of diff expressed enhancer RNA VHL/EV in hypoxia
	FDR of diff expressed enhancer RNA VHL/EV in hypoxia
	log2FoldChange of diff expressed enhancer RNA VHL/EV in Normoxia
	FDR of diff expressed enhancer RNA VHL/EV in Normoxia
	hypo_count_5kb_downstream
	hyper_count_5kb_downstream
	p_value_5kb_downstream
	methylation_enrichment_5kb_downstream
	total_cpg_5kb_downstream
