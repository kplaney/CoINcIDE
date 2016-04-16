
#NOTE: I found wget plus the broad http link gave garbage. I had to download to my local computer and then scp the files over.
#
#the .gmt files I downloaded  from http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C7 were:
#"MSigDB_onco_symbols","MSigDB_CanPath_symbols",
                                         #"MSigDB_TFT_symbols","MSigDB_immun_symbols",
                                         #"MSigDB_cancerNeigh_symbols"

#NOTE: I downloaded these in 2014. newer versions may exist now,
#but I do have the final output from here saved.
library("GSA");
library("limma");

#NOTE: need to change the paths here to where you'd like to save these lists.
MSigDBPath <- "/home/data/MSigDB/"

#set to directory where downloaded files:
setwd(MSigDBPath);

#look for .gmt files:
gmtFiles_symb <- list.files("./",pattern="symbols.gmt");
gmtNames_symb <- strsplit2(gmtFiles_symb,split="\\.")[,1];
gmtFiles_entrez <- list.files("./",pattern="entrez.gmt");
gmtNames_entrez <- strsplit2(gmtFiles_entrez,split="\\.")[,1];
MSigDB_lists_symb <- list();
MSigDB_lists_merged_symb <- list();
MSigDB_lists_entrez <- list();
MSigDB_lists_merged_entrez <- list();

count <- 0;
for(g in 1:length(gmtFiles_symb)){
  
  tmp <- GSA.read.gmt(gmtFiles_symb[g]);
  MSigDB_lists_symb[[g]] <- tmp$genesets;
  names(MSigDB_lists_symb[[g]]) <- tmp$geneset.names;
  #count <- count + length(MSigDB_lists[[g]]);
  #MSigDB_lists_merged[c((count-length(MSigDB_lists[[g]])+1):count)] <- MSigDB_lists[[g]];
  MSigDB_lists_merged_symb <- append(MSigDB_lists_merged_symb,MSigDB_lists_symb[[g]]);
  tmp <- GSA.read.gmt(gmtFiles_entrez[g]);
  MSigDB_lists_entrez[[g]] <- tmp$genesets;
  names(MSigDB_lists_entrez[[g]]) <- tmp$geneset.names;
  MSigDB_lists_merged_entrez <- append(MSigDB_lists_merged_entrez,MSigDB_lists_entrez[[g]]);
  
}

names(MSigDB_lists_symb) <- gmtNames_symb;
names(MSigDB_lists_entrez) <- gmtNames_entrez;

#REACTOME pathways are already included in some of the cancer pathway lists. so don't re-included the REACTOME here!
#many of these immun symbol lists are derived from GSE studies - interesting!
GSEA_base_MSigDB_lists_symb <- MSigDB_lists_symb[c("MSigDB_onco_symbols","MSigDB_CanPath_symbols",
                                         "MSigDB_TFT_symbols","MSigDB_immun_symbols",
                                         "MSigDB_cancerNeigh_symbols")];

GSEA_base_MSigDB_lists_merged <- list();

for(e in 1:length(GSEA_base_MSigDB_lists_symb)){
  
  GSEA_base_MSigDB_lists_merged <- append(GSEA_base_MSigDB_lists_merged,GSEA_base_MSigDB_lists_symb[[e]]);
  if(any(duplicated(names(GSEA_base_MSigDB_lists_merged)))){
    cat("loop:",e,"\n")
    stop("error! Some pathway names are duplicated")
  }
}


save(MSigDB_lists_symb,file=paste0(MSigDBPath ,"/MSigDB_lists_symb.RData.gzip"),compress="gzip");
save(GSEA_base_MSigDB_lists_symb,file=paste0(MSigDBPath ,"/GSEA_base_MSigDB_lists_symb.RData.gzip",)compress="gzip");
save(MSigDB_lists_entrez,file=paste0(MSigDBPath ,"/MSigDB_lists_entrez.RData.gzip"),compress="gzip");
save(MSigDB_lists_merged_entrez,file=paste0(MSigDBPath ,"/MSigDB_lists_merged_entrez.RData.gzip"),compress="gzip");
save(MSigDB_lists_merged_symb,file=paste0(MSigDBPath ,"/MSigDB_lists_merged_symb.RData.gzip"),compress="gzip");
save(GSEA_base_MSigDB_lists_merged,file=paste0(MSigDBPath ,"/GSEA_base_MSigDB_lists_merged.RData.gzip"),compress="gzip");
