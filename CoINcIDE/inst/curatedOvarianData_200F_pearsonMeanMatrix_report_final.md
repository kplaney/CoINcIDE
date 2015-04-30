curatedOvarianData meta-clustering with 200 (meta-rank) features using kmeans consensus and pearson mean matrix correlation
========================================================
**Data loading**


This is a report for the curatedOvarianData meta-clustering. Initial clustering was done using k-means and selecting k via the the rounded PACR metric with nstart=1.

First, I loaded up my data just past the CoINcIDE" "getAdjMatrices() function (that's the function that takes  a really long time to run and is hogging up server space, as this is the part that computes the edge p-values):


```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     Filter, Find, Map, Position, Reduce, anyDuplicated, append,
##     as.data.frame, as.vector, cbind, colnames, do.call,
##     duplicated, eval, evalq, get, intersect, is.unsorted, lapply,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, rank, rbind, rep.int, rownames, sapply, setdiff,
##     sort, table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
## 
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
## 
## Loading required package: annotate
## Loading required package: AnnotationDbi
## Loading required package: GenomeInfoDb
## Loading required package: graph
## 
## Attaching package: 'graph'
## 
## The following objects are masked from 'package:igraph':
## 
##     degree, edges
## 
## 
## Attaching package: 'plyr'
## 
## The following object is masked from 'package:graph':
## 
##     join
## 
## Loading required package: bitops
## Loading required package: GOstats
## Loading required package: Category
## Loading required package: Matrix
## 
## Attaching package: 'Matrix'
## 
## The following objects are masked from 'package:base':
## 
##     crossprod, tcrossprod
## 
## Loading required package: GO.db
## Loading required package: DBI
## 
## 
## Attaching package: 'GOstats'
## 
## The following object is masked from 'package:AnnotationDbi':
## 
##     makeGOGraph
## 
## 
## Attaching package: 'RDAVIDWebService'
## 
## The following object is masked from 'package:GSEABase':
## 
##     ids
## 
## The following object is masked from 'package:AnnotationDbi':
## 
##     species
## 
## The following objects are masked from 'package:igraph':
## 
##     is.connected, membership
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     counts
```

**Input variables used to derive adjacency matrix**



```
## Input variables used to derive the adjacency matrix:
```

```
##                  date edgeMethod numParallelCores minTrueSimilThresh
## 1 2015-04-29 07:07:36    pearson                3                0.3
##   maxTrueSimilThresh  sigMethod maxNullFractSize numSims
## 1                Inf meanMatrix              0.2     500
##   includeRefClustInNull fractFeatIntersectThresh numFeatIntersectThresh
## 1                  TRUE                      0.8                    150
##   clustSizeThresh clustSizeFractThresh
## 1               5                 0.05
```

```
## There were 92 total input clusters from 24 studies
## The total number of input features was 240
## Across the entire square (nonsymmetric) p-value matrix, there are 1610 pvalues less than or equal to .1
## Across the entire square (nonsymmetric) p-value matrix, there are 1429 pvalues less than or equal to .05
## Across the entire square (symmetric) similarity matrix, there are 1426 similarities greater than or equal to 0.5
```

**Network Analysis**

As we can see in these plots, 4 meta-clusters remained after edge filtering. 

```
## 17 clusters dropped because they were below the clust size thresh clustSizeThresh
##           threshold of 5
## A total of 45 clusters removed because they have no significant edges.
## A total of 47 clusters have significant edges.
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-3.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-4.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-5.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-6.png) 

```
## Overall network stats (I believe the origNumClusters variable may be off-need to debug:
```

```
##  numCommunities     numClusters origNumClusters        numEdges 
##               5              38              38             155 
##      numStudies 
##              18
```

**Gene meta-rank Analysis**


I ranked genes within each meta-cluster for all samples, and then ran a Kruskal test to see which genes significantly stratified/differentiated patients across the 4 meta-clusters. I still need to implement GSEA; it turns out there's a base GSEA package in Biocondcutor so I've decided to just adapt my code and use their baseline functions.


In the heatmap: red means that gene was ranked high in terms of expression level for patients in that meta-cluster (I took the median rank across all samples in a meta-cluster to create the heatmap. Only significant genes are shown in the heatmap but at over 80 significant genes, of course the gene names are illegible...I print out the top 20 genes below.) The actual meta-cluster numbers are not 1:4 because meta-clusters with only 2 studies were thresholded out earlier.


```
## There are 83 genes with an FDR corrected p-value below .05 whose ranks significantly differed among all of the meta-clusters (using Kruskal's test)
## Red in heatmap means genes were ranked higher across all samples in that meta-cluster.
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

```
## Running effect size analysis:
```

```
## Error in colnames(computeMetaclustEffectSizesOutput$summWilcoxon_qvalue): error in evaluating the argument 'x' in selecting a method for function 'colnames': Error: object 'computeMetaclustEffectSizesOutput' not found
```

**GSEA Analysis**

It turns out "GSEABase" isn't all that great..so just using some Broad gene sets I already downloaded and the standard hypergeometric tests:


```
## 
##  163  genes for commmunity number  1 (name)  1  have a positive effect size above .2 with a wilcoxon rank q-value below .1
##  index 9  gene  CDKN2A  had multiple records returned.
## 
##  index 24  gene  GLDC  had multiple records returned.
## 
##  index 33  gene  COL3A1  had multiple records returned.
## 
##  index 38  gene  MMP11  had multiple records returned.
## 
##  index 45  gene  CFB  had multiple records returned.
## 
##  index 62  gene  PDGFRA  had multiple records returned.
## 
##  index 72  gene  TDO2  had multiple records returned.
## 
##  index 83  gene  CRIP1  had multiple records returned.
## 
##  index 89  gene  COL6A2  had multiple records returned.
## 
##  index 93  gene  COL6A3  had multiple records returned.
## 
##  index 97  gene  WT1  had multiple records returned.
## 
##  index 98  gene  C1S  had multiple records returned.
## 
##  index 99  gene  C1QB  had multiple records returned.
## 
##  index 105  gene  HLA-DPB1  had multiple records returned.
## 
##  index 106  gene  IFI27  had multiple records returned.
## 
##  index 108  gene  CCL5  had multiple records returned.
## 
##  index 116  gene  HLA-DPA1  had multiple records returned.
## 
##  index 125  gene  PMP22  had multiple records returned.
## 
##  index 126  gene  LY6E  had multiple records returned.
## 
##  index 127  gene  CBS  had multiple records returned.
## 
##  index 144  gene  COL1A2  had multiple records returned.
## 
##  index 155  gene  COL6A1  had multiple records returned.
## 
##  index 159  gene  IL10RB  had multiple records returned.
```

```
## Connecting to DAVID web service. Occassionally it times out and you need to re-run this function.
## Significantly enriched DAVID functional groups for overexpressed genes in meta-cluster 1 with soft q-value thresholds:
```

```
## Secreted signal signal peptide GO:0044421~extracellular region part GO:0005576~extracellular region extracellular matrix GO:0005578~proteinaceous extracellular matrix GO:0031012~extracellular matrix GO:0009611~response to wounding disulfide bond GO:0005615~extracellular space disulfide bond glycoprotein GO:0005201~extracellular matrix structural constituent GO:0044420~extracellular matrix part GO:0006954~inflammatory response GO:0007155~cell adhesion GO:0022610~biological adhesion glycosylation site:N-linked (GlcNAc...) hydroxylation cell adhesion GO:0005581~collagen hydroxylysine triple helix GO:0001501~skeletal system development hydroxyproline calcium binding hsa04512:ECM-receptor interaction GO:0001568~blood vessel development GO:0001944~vasculature development IPR001811:Small chemokine, interleukin-8-like pyroglutamic acid chemotaxis GO:0006952~defense response SM00199:SCY GO:0042330~taxis GO:0006935~chemotaxis GO:0008009~chemokine activity region of interest:Triple-helical region GO:0048545~response to steroid hormone stimulus cytokine GO:0042379~chemokine receptor binding trimer IPR008160:Collagen triple helix repeat collagen acute phase GO:0006955~immune response heterodimer GO:0042060~wound healing GO:0030198~extracellular matrix organization inflammatory response GO:0043062~extracellular structure organization IPR018048:Small chemokine, C-X-C, conserved site GO:0007568~aging GO:0007626~locomotory behavior GO:0016477~cell migration GO:0032963~collagen metabolic process short sequence motif:Cell attachment site GO:0030247~polysaccharide binding GO:0001871~pattern binding GO:0044259~multicellular organismal macromolecule metabolic process plasma GO:0007610~behavior GO:0005125~cytokine activity GO:0010033~response to organic substance GO:0048870~cell motility GO:0051674~localization of cell GO:0019838~growth factor binding GO:0051384~response to glucocorticoid stimulus GO:0005539~glycosaminoglycan binding GO:0044236~multicellular organismal metabolic process hsa04510:Focal adhesion duplication GO:0048407~platelet-derived growth factor binding GO:0048514~blood vessel morphogenesis GO:0031960~response to corticosteroid stimulus IPR001089:Small chemokine, C-X-C GO:0006928~cell motion GO:0005583~fibrillar collagen cell binding GO:0007267~cell-cell signaling GO:0009719~response to endogenous stimulus GO:0004175~endopeptidase activity Protease calcium-binding region:1; low affinity GO:0009725~response to hormone stimulus GO:0030199~collagen fibril organization cleavage on pair of basic residues Serine protease GO:0005604~basement membrane Pyrrolidone carboxylic acid 109.Chemokine_families calcium-binding region:2; high affinity domain:Peptidase S1 region of interest:Nonhelical region propeptide:C-terminal propeptide propeptide:Activation peptide GO:0040007~growth calcium GO:0005518~collagen binding GO:0004252~serine-type endopeptidase activity GO:0005198~structural molecule activity inflammation IPR001314:Peptidase S1A, chymotrypsin PIRSF002522:CXC chemokine hsa04060:Cytokine-cytokine receptor interaction GO:0042127~regulation of cell proliferation zymogen IPR018114:Peptidase S1/S6, chymotrypsin/Hap, active site GO:0030574~collagen catabolic process metal ion-binding site:Zinc 2; in inhibited form proteoglycan domain:Fibrillar collagen NC1 IPR013032:EGF-like region, conserved site GO:0010035~response to inorganic substance IPR001254:Peptidase S1 and S6, chymotrypsin/Hap GO:0002526~acute inflammatory response IPR013787:S100/CaBP-9k-type, calcium binding, subdomain GO:0031099~regeneration IPR000885:Fibrillar collagen, C-terminal IPR001751:S100/CaBP-9k-type, calcium binding GO:0008236~serine-type peptidase activity egf-like domain GO:0008233~peptidase activity GO:0017171~serine hydrolase activity EF hand IPR017891:Insulin-like growth factor binding protein, N-terminal GO:0009991~response to extracellular stimulus GO:0044243~multicellular organismal catabolic process GO:0000302~response to reactive oxygen species GO:0050840~extracellular matrix binding IPR016060:Complement control module PIRSF002353:S-100 protein metal ion-binding site:Zinc 2; catalytic heparin-binding active site:Charge relay system GO:0048705~skeletal system morphogenesis GO:0070011~peptidase activity, acting on L-amino acid peptides SM00038:COLFI IPR002473:Small chemokine, C-X-C/Interleukin 8 GO:0001503~ossification transmembrane protein GO:0001558~regulation of cell growth GO:0032355~response to estradiol stimulus collagen degradation hsa04610:Complement and coagulation cascades GO:0032101~regulation of response to external stimulus smooth muscle GO:0048589~developmental growth GO:0005178~integrin binding GO:0030246~carbohydrate binding GO:0060348~bone development SM00020:Tryp_SPc GO:0006979~response to oxidative stress IPR000716:Thyroglobulin type-1 domain:IGFBP N-terminal disease mutation Plasminogen activation propeptide:Removed in mature form domain:Fibronectin type-II 2 domain:Fibronectin type-II 1 GO:0008285~negative regulation of cell proliferation Bethlem myopathy Ullrich congenital muscular dystrophy PIRSF002259:collagen VI GO:0005509~calcium ion binding IPR000867:Insulin-like growth factor-binding protein, IGFBP GO:0008201~heparin binding GO:0007584~response to nutrient GO:0007596~blood coagulation GO:0050817~coagulation GO:0050878~regulation of body fluid levels metalloproteinase blocked amino end propeptide:N-terminal propeptide SM00211:TY serine proteinase lipid moiety-binding region:GPI-anchor amidated serine GO:0043627~response to estrogen stimulus copper GO:0007599~hemostasis GO:0009615~response to virus GO:0007565~female pregnancy domain:VWFA 3 GO:0040008~regulation of growth GO:0032989~cellular component morphogenesis GO:0046658~anchored to plasma membrane GO:0031667~response to nutrient levels GO:0000904~cell morphogenesis involved in differentiation SM00121:IB GO:0034097~response to cytokine stimulus GO:0006959~humoral immune response bone IPR001881:EGF-like calcium-binding IPR018097:EGF-like calcium-binding, conserved site GO:0031175~neuron projection development IPR000742:EGF-like, type 3 IPR000837:Fos transforming protein GO:0005520~insulin-like growth factor binding glycosylation site:O-linked (Gal...) PIRSF001719:fos transforming protein IPR006026:Peptidase, metallopeptidases GO:0007517~muscle organ development GO:0009628~response to abiotic stimulus hsa04620:Toll-like receptor signaling pathway GO:0030334~regulation of cell migration IPR006210:EGF-like GO:0004867~serine-type endopeptidase inhibitor activity GO:0030335~positive regulation of cell migration GO:0009612~response to mechanical stimulus Ehlers-Danlos syndrome gpi-anchor hsa04062:Chemokine signaling pathway Serine protease inhibitor GO:0048666~neuron development GO:0051272~positive regulation of cell motion GO:0031589~cell-substrate adhesion GO:0040017~positive regulation of locomotion GO:0004866~endopeptidase inhibitor activity GO:0008015~blood circulation GO:0003013~circulatory system process GO:0048584~positive regulation of response to stimulus SM00235:ZnMc GO:0002684~positive regulation of immune system process domain:Thyroglobulin type-1 GO:0042246~tissue regeneration short sequence motif:Cysteine switch GO:0040012~regulation of locomotion GO:0051270~regulation of cell motion pharmaceutical SM00179:EGF_CA GO:0030414~peptidase inhibitor activity GO:0001525~angiogenesis GO:0033273~response to vitamin GO:0000902~cell morphogenesis GO:0005507~copper ion binding homodimer IPR006209:EGF site:Susceptible to oxidation sequence variant IPR000562:Type II fibronectin, collagen-binding IPR001818:Peptidase M10A and M12B, matrixin and adamalysin Growth factor binding domain:VWFA 1 GO:0042981~regulation of apoptosis GO:0042552~myelination GO:0004857~enzyme inhibitor activity GO:0043067~regulation of programmed cell death GO:0032844~regulation of homeostatic process GO:0010941~regulation of cell death GO:0051216~cartilage development GO:0048812~neuron projection morphogenesis GO:0051241~negative regulation of multicellular organismal process GO:0050727~regulation of inflammatory response GO:0042493~response to drug SM00181:EGF GO:0034329~cell junction assembly hsa04621:NOD-like receptor signaling pathway sushi GO:0006956~complement activation GO:0007272~ensheathment of neurons GO:0051591~response to cAMP GO:0008366~axon ensheathment GO:0019748~secondary metabolic process microsome SM00059:FN2 GO:0016337~cell-cell adhesion GO:0002541~activation of plasma proteins involved in acute inflammatory response protease inhibitor IPR000152:EGF-type aspartate/asparagine hydroxylation conserved site GO:0051412~response to corticosterone stimulus IPR000233:Cadherin cytoplasmic region GO:0050865~regulation of cell activation muscle protein GO:0045137~development of primary sexual characteristics GO:0002683~negative regulation of immune system process IPR002477:Peptidoglycan binding-like hsa05200:Pathways in cancer IPR000436:Sushi/SCR/CCP GO:0016638~oxidoreductase activity, acting on the CH-NH2 group of donors IPR018247:EF-HAND 1 angiogenesis IPR003129:Laminin G, thrombospondin-type, N-terminal GO:0032403~protein complex binding IPR018486:Hemopexin/matrixin, conserved site IPR018487:Hemopexin/matrixin, repeat IPR000585:Hemopexin/matrixin GO:0004833~tryptophan 2,3-dioxygenase activity complement pathway polymorphism quinoprotein topaquinone innate immunity extracellular protein Signal transduction inhibitor immune response sialoglycoprotein skin ltq aortic aneurysm lipoprotein basement membrane 
## 
## 
##  44  genes for commmunity number  2 (name)  9  have a positive effect size above .2 with a wilcoxon rank q-value below .1
##  index 2  gene  CDKN2A  had multiple records returned.
## 
##  index 9  gene  GLDC  had multiple records returned.
## 
##  index 14  gene  CFB  had multiple records returned.
## 
##  index 29  gene  IFI27  had multiple records returned.
## 
##  index 36  gene  LY6E  had multiple records returned.
```

```
## Connecting to DAVID web service. Occassionally it times out and you need to re-run this function.
## Significantly enriched DAVID functional groups for overexpressed genes in meta-cluster 2 with soft q-value thresholds:
```

```
## disulfide bond signal signal peptide Secreted disulfide bond lipid moiety-binding region:GPI-anchor amidated serine cleavage on pair of basic residues GO:0031225~anchored to membrane GO:0005576~extracellular region glycoprotein microsome propeptide:Removed in mature form glycosylation site:N-linked (GlcNAc...) gpi-anchor lipoprotein oxidoreductase chloride transmembrane protein 
## 
## 
##  53  genes for commmunity number  3 (name)  10  have a positive effect size above .2 with a wilcoxon rank q-value below .1
##  index 1  gene  C7  had multiple records returned.
## 
##  index 10  gene  PDGFRA  had multiple records returned.
## 
##  index 33  gene  PMP22  had multiple records returned.
## 
##  index 41  gene  MYH11  had multiple records returned.
## 
##  index 45  gene  COL6A1  had multiple records returned.
## 
##  index 50  gene  IL10RB  had multiple records returned.
```

```
## Connecting to DAVID web service. Occassionally it times out and you need to re-run this function.
## Significantly enriched DAVID functional groups for overexpressed genes in meta-cluster 3 with soft q-value thresholds:
```

```
## signal signal peptide Secreted GO:0044421~extracellular region part transmembrane protein glycosylation site:N-linked (GlcNAc...) GO:0048514~blood vessel morphogenesis GO:0005578~proteinaceous extracellular matrix GO:0001501~skeletal system development GO:0007155~cell adhesion GO:0022610~biological adhesion glycoprotein GO:0031012~extracellular matrix calcium binding GO:0001568~blood vessel development GO:0005576~extracellular region GO:0001944~vasculature development GO:0003006~reproductive developmental process disulfide bond extracellular matrix GO:0044420~extracellular matrix part GO:0007507~heart development cell adhesion GO:0046546~development of primary male sexual characteristics GO:0046661~male sex differentiation GO:0005604~basement membrane deafness GO:0030334~regulation of cell migration basement membrane GO:0010035~response to inorganic substance GO:0005615~extracellular space disease mutation smooth muscle GO:0030485~smooth muscle contractile fiber 
## 
## 
##  41  genes for commmunity number  4 (name)  12  have a positive effect size above .2 with a wilcoxon rank q-value below .1
##  index 3  gene  C7  had multiple records returned.
## 
##  index 10  gene  SERPINA1  had multiple records returned.
## 
##  index 14  gene  PDGFRA  had multiple records returned.
## 
##  index 35  gene  MYH11  had multiple records returned.
```

```
## Connecting to DAVID web service. Occassionally it times out and you need to re-run this function.
## Significantly enriched DAVID functional groups for overexpressed genes in meta-cluster 4 with soft q-value thresholds:
```

```
## Secreted signal signal peptide GO:0044421~extracellular region part GO:0005576~extracellular region GO:0009611~response to wounding GO:0005615~extracellular space disulfide bond GO:0048545~response to steroid hormone stimulus extracellular matrix GO:0006954~inflammatory response GO:0005578~proteinaceous extracellular matrix GO:0009725~response to hormone stimulus disulfide bond GO:0031012~extracellular matrix GO:0009719~response to endogenous stimulus Serine protease inhibitor calcium-binding region:1; low affinity GO:0042127~regulation of cell proliferation GO:0006952~defense response calcium-binding region:2; high affinity serine proteinase inhibitor calcium binding protease inhibitor GO:0042060~wound healing GO:0005201~extracellular matrix structural constituent PIRSF002353:S-100 protein GO:0004867~serine-type endopeptidase inhibitor activity IPR013787:S100/CaBP-9k-type, calcium binding, subdomain GO:0010035~response to inorganic substance IPR001751:S100/CaBP-9k-type, calcium binding GO:0010033~response to organic substance glycosylation site:N-linked (GlcNAc...) growth factor GO:0045768~positive regulation of anti-apoptosis glycoprotein site:Reactive bond IPR013032:EGF-like region, conserved site GO:0004857~enzyme inhibitor activity cytokine GO:0004866~endopeptidase inhibitor activity EF hand GO:0030414~peptidase inhibitor activity GO:0008083~growth factor activity Antimicrobial transmembrane protein egf-like domain 
## 
## 
##  49  genes for commmunity number  5 (name)  13  have a positive effect size above .2 with a wilcoxon rank q-value below .1
##  index 2  gene  C7  had multiple records returned.
## 
##  index 9  gene  GLDC  had multiple records returned.
## 
##  index 11  gene  APOA1  had multiple records returned.
## 
##  index 16  gene  CFB  had multiple records returned.
## 
##  index 28  gene  WT1  had multiple records returned.
## 
##  index 36  gene  LY6E  had multiple records returned.
```

```
## Connecting to DAVID web service. Occassionally it times out and you need to re-run this function.
## Significantly enriched DAVID functional groups for overexpressed genes in meta-cluster 5 with soft q-value thresholds:
```

```
## disulfide bond signal signal peptide Secreted disulfide bond glycoprotein glycosylation site:N-linked (GlcNAc...) GO:0005576~extracellular region GO:0048732~gland development lipoprotein GO:0050767~regulation of neurogenesis lipid moiety-binding region:GPI-anchor amidated serine GO:0051960~regulation of nervous system development GO:0044057~regulation of system process GO:0060284~regulation of cell development GO:0031225~anchored to membrane lipid binding GO:0007584~response to nutrient GO:0035270~endocrine system development plasma propeptide:Removed in mature form GO:0044421~extracellular region part gpi-anchor transmembrane protein GO:0031226~intrinsic to plasma membrane
```
**Survival Analyses**
I'm still working on determing which long-term and binary variables have the least amount of NAs across the samples in these meta-clusters, but it does look like the binary vital_status variable (alive or dead) and continous days to death variables provide survival curves that significantly stratify patients. It was less significant when I used a 5-year cutoff.


```
## 
## Attaching package: 'survival'
## 
## The following object is masked from 'package:RDAVIDWebService':
## 
##     cluster
## 
## This function assumes samples/patient clinical data rows are not duplicated
## coxfit summary for overall survival.
```

```
## Call:
## coxph(formula = OverallSurvival ~ groupings, data = Survival)
## 
##   n= 1316, number of events= 686 
##    (567 observations deleted due to missingness)
## 
##                coef exp(coef)  se(coef)      z Pr(>|z|)  
## groupings -0.012824  0.987258  0.006937 -1.849   0.0645 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##           exp(coef) exp(-coef) lower .95 upper .95
## groupings    0.9873      1.013    0.9739     1.001
## 
## Concordance= 0.526  (se = 0.011 )
## Rsquare= 0.003   (max possible= 0.999 )
## Likelihood ratio test= 3.47  on 1 df,   p=0.06246
## Wald test            = 3.42  on 1 df,   p=0.0645
## Score (logrank) test = 3.42  on 1 df,   p=0.06426
```

```
## kaplan meier p-value for overall survival:
```

```
## [1] 0.07964773
```

```
## calculating the sign of the survival relationship
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-2.png) 

```
## coxfit summary for survival cutoff at 5years:
```

```
## Call:
## coxph(formula = OverallSurvivalCutoff ~ groupings, data = SurvivalCutoff)
## 
##   n= 1316, number of events= 602 
##    (567 observations deleted due to missingness)
## 
##                coef exp(coef)  se(coef)      z Pr(>|z|)  
## groupings -0.013306  0.986783  0.007382 -1.802   0.0715 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##           exp(coef) exp(-coef) lower .95 upper .95
## groupings    0.9868      1.013    0.9726     1.001
## 
## Concordance= 0.526  (se = 0.011 )
## Rsquare= 0.003   (max possible= 0.998 )
## Likelihood ratio test= 3.3  on 1 df,   p=0.06924
## Wald test            = 3.25  on 1 df,   p=0.07149
## Score (logrank) test = 3.25  on 1 df,   p=0.07121
```

```
## kaplan meier p-value for survival cutoff:
```

```
## [1] 0.07995688
```

```
##       
##        clearcell endo  mix mucinous other  ser undifferentiated <NA>
##   1           13   32    1        3    25 1003                2   35
##   9            1    0    3        0     0   40                2   10
##   10           6   27    0       25    43   21                1   51
##   12          20   34    0       16     2    8                0    0
##   13           2   17    0        0     2  437                0    1
##   <NA>         0    0    0        0     0    0                0    0
```

```
##       
##          1   2   3   4 <NA>
##   1     38  23 802 156   95
##   9      0   0  20   5   31
##   10    22  13  31  14   94
##   12    48  13  15   4    0
##   13    33  35 333  53    5
##   <NA>   0   0   0   0    0
```

```
##       
##        norecurrence recurrence <NA>
##   1             208        417  489
##   9               0          0   56
##   10             12         10  152
##   12             20         14   46
##   13            191        257   11
##   <NA>            0          0    0
```

```
##       
##          1   2   3   4 <NA>
##   1     18 240 545   6  305
##   9      0  11  34   1   10
##   10    14  26  36   2   96
##   12    18  18  18   0   26
##   13    14  99 327   1   18
##   <NA>   0   0   0   0    0
```

```
## groupings
##    1    9   10   12   13 
## 1114   56  174   80  459
```
