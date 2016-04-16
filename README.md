# Coincide

This directory/repo contains both a library and extensive scripts for clustering across multiple datasets. 
The method is called CoINcIDE, or Clustering INtra and Inter DatasEts.

The Coincide_fullPackage/ directory contains the "full" package in that it contains all main and helper functions 
used in my dissertation and the recent publication https://www.ncbi.nlm.nih.gov/pubmed/26961683 in Genome Medicine. 

The functions are not fully documented in this package yet, but Coincide_GenomeMedicine_dissertation_scripts/ 
contains all scripts run for the Genome Medicine analysis, and thus provides actual examples of all of the functiosn 
in the package. The README file in this directory contains many more details on the order of scripts I ran, and 
all of the scripts, from pre-processing to analysis and visualization of results.

I am working on "stripping down" the full package into a more manageable and well-documented package to submit to 
Bioconductor, and then slowly add on some of the more data-type-specific functions I used that are tailored to the 
datasets I used in the publication, curatedBreastData and curatedOvarianData, which are on Biocondcutor.

I am not currently in academia anymore (i.e. I defended my thesis and graduated), but I do my best to respond to emails quickly, and am still developing 
extensions to Coincide in my spare time. For example, Coincide can also be used as a new consensus clustering 
framework to select the number of clusters on a single dataset. If you are interested in collaborating on this and 
other work, please do let me know.

email: catherine.planey@gmail.com


