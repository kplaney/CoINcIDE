These scripts use the FULL package Coincide version Coincide_0.99.2.tar.gz.  "Full" simply means extra helper functions 
specific to the breast and ovarian cancer analyses produced for the March 2016 publication of Coincide, and my closely
related doctoral dissertation from the Stanford Biomedical Informatics program, were included in this package.

I will be the first to admit that this full package is not documented in terms of R package documentation files,
but the Coincide_0.99.2.tar.gz code in the R directory of this tarbell contains copious (perhaps too many for 
the average user) comments.  The analyses scripts can also be followed to see how the Genome Medicine and my 
dissertation analyes were run from start to finish.

Note that I attempt to write path names at the top of the scripts, but given this code was written over my 
entire PhD, I sometimes did hard code in paths, especially earlier in my PhD.  The user can search for "/home/ywrfc09" 
(my user directory) to replace these paths with their own paths. However, I made a good-faith report before posting 
these scripts to add path variables at the top of each script. Use the same path variables for all of your breast, and 
all of your ovarian analyses, as I assume this as you "link" analyses between the different scripts.

You'll note that many of these scripts seem repetitive, and could have been run in a loop. I coded them this way because I 
usually just ran them in a naively parallel fashion by kicking off different batch runs or screen sessions (the joys of 
computing in academia!)

The breast and ovarian scripts were run in a fairly similar fashion, using curatedBreastData (a package I created)
and curatedOvarianData from Bioconductor.  The processing scripts for both diseases respectively are:

-breastProcessAndGeneFeatures_script.R
-ovarianProcessAndGeneFeatures_script.R

First things first: my research involved clustering, which means we have to chose both the baseline clustering method 
and the way we select the number of clusters. I used the breast cancer dataset to help select these parameters, using
the script:

-Concide_selectK_cluster_script.R

I slightly tweaked the Bioconductor ConsensusClusterPlus code to not have it spit out lots of plots on each loop,
and I also added an input variable to choose how many starts for k-means to use. Other than that, the code is the 
same from the 2015 version of the ConsensusClusterPlus package. I use the Coincide clustMatrixListWrapper() mainly here, 
and in the cluster scripts below, to run all of this on all breast (or ovarian) datasets at once.

The results (see the code at the bottom of script) are reported in Supplementary Table 2 in the Genome Medicine publication.

Now we can cluster! The above results helped me to decide to pick Hartigan-Wong's k-means consensus clustering with 
one start, so that's what you'll see used in subsequent scripts, except for the centroid clustering, described below.

Note that clustering ~20 datasets takes a few hours to run both the breast and then the ovarian datasets, using 
the scripts:

-Coincide_breastCluster_script.R
-Coincide_ovarianCluster_script.R

The initial analyses in this script are for semi-supervised centroid clustering using the PAM50 centroids. This 
was done to create very clear, obvious clusters (to make sure Coincide could find well-defined clusters on real data 
as a gut check.) 

Some small RData objects are used here to cluster on a well-known breast cancer feature set, the PAM50 genes:
-pam50_centroids_updatedSymbols.RData: pam50 gene symbol centroids (see reference in online publication) updated using HGNChelper package.
-pam50Short_genes.RData: the 35 pam50 genes taken from pam50_centroids_updatedSymbols.RData that intersect all datasets.

For such small datasets included in the repo, I just load them and assume they're in your current working directory.

I also included in the breast clustering script my work to cluster the entire concatenated/merged matrix, using a few      
different normalization schemes, as this is the current way to cluster across multiple datasets.

OK we have all of our clusters for each dataset. Now let's run the core functions of Coincide. I did this by using two Rscript 
command line files, and a third file that called them for each of the clustering analyses from the two files above. 

-Coincide_baseAnalysisNullFDRScript.R
-Coincide_baseAnalysisScript.R
-Coincide_breastAndOvarian_script.R


Next, as a pseudo basline comparison to later Coincide analyses, I analyzed the breast concatenated/merged 
matrix clusterings using the script:

-Coincide_breastMergedAnalysis_script.R

Next, a lot of up-front work was done on simulated data to prove Coincide behaves as "expected".  Functions to 
both create the simulated clusters and run Coincide on it can be found in the scripts below:

(a quick note: below this line, I'm still working on not having paths hard-coded to make it 
easier for other folks to copy and paste/run the code):
-Coincide_tissueSimulation_script.R
-Coincide_simulationVisualizations_script.R

Finally, on to the visualization and interpretation of the results on real data. I used a wrapper function that is more closely tailored to my specific breast and ovarian analyses that uses the broader 
Coincide package functions to analyze survival results for the breast and ovarian data. This wrapper script is:

-Coincide_metaFeaturesAnalysisWrapper.R

Scripts that use this wrapper function to interpret the output for the breast and ovarian functions can be found at:

-Coincide_breastMergedAnalysis_script.R
-Coincide_ovarianMetaFeatures_analysis_withTop20Genes_script.R

Two additional scripts to dig into further detail on specific aspects of the ovarian analyses are:

-Coincide_ovarianMetaFeaturesAnalysis_withTop20Genes_serousOnly_script.R
-Coincide_ovarianWithTop20Genes_MM7_min30OutcomesPerSubtypes_script.R 

And a script showing how I inspected druggable genes that had a high meta-ranked effect size is:

-Coincide_findDruggableTopESGenes_script. R

Note that I used GSEA on the genes with a high meta-ranked effect size in these metaFeatures scripts, 
and that requires gene lists to be inputted (this is the refGeneListDir="~/GSEA_base_MSigDB_lists_merged.RData.gzip" 
keyword in the metaFeaturesAnalysisWrapper function.)  I specifiy the Broad institute gene lists 
in the supplemental sections of the Genome Medicine publication, but if the RData object is not in the Github repo,
it's because I was nervous about its size.  However, I do have these gene lists in an RData object for complete 
reproducibility, as the Broad Institute may change the gene lists over time. You can reach out to the email below 
if you wish to obtain this RData object. The script I used to format the Broad GMT files is:

Coincide_readIn_MSigDB_GMTfiles.R

Questions? Contact Katie Planey at katie.planey@gmail.com
