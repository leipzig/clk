#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(stringr)))

df<-read.table(args[1],header=TRUE)

if(args[2]=='omit'){
    df$IncLevel1<-str_replace_all(pattern = ',?NA,?',replacement = '',string = df$IncLevel1)
    df$IncLevel2<-str_replace_all(pattern = ',?NA,?',replacement = '',string = df$IncLevel2)
}
if(args[2]=='zero'){
        df$IncLevel1<-str_replace_all(pattern = 'NA',replacement = '0.0',string = df$IncLevel1)
    df$IncLevel2<-str_replace_all(pattern = 'NA',replacement = '0.0',string = df$IncLevel2)
}
if(args[2]=='skip'){
    df %>% dplyr::filter(!str_detect(IncLevel1,pattern='NA')) %>% dplyr::filter(!str_detect(IncLevel2,pattern='NA')) -> df
}
df %>% dplyr::filter(FDR<=0.05) -> df

df$GeneID<-paste0('"',df$GeneID,'"')
df$geneSymbol<-paste0('"',df$geneSymbol,'"')
#ID      GeneID  geneSymbol      chr     strand  exonStart_0base exonEnd upstreamES      upstreamEE      downstreamES    downstreamEE    ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2IncFormLen       SkipFormLen     PValue  FDR     IncLevel1       IncLevel2       IncLevelDifference
#0       "ENSG00000096696.13"    "DSP"   chr6    +       7555717 7555820 7541574 7542085 7558115 7558264 0       0,427,695,889   0,0,0,0 256,0,638,593   0,0,2,4 -2      -1      1       1.0     NA,1.0,1.0,1.0       1.0,NA,0.994,0.987      0.006

if(nrow(df)>0){
    write.table(df,file=str_replace(args[1],'.txt','.filtered.txt'),quote=FALSE,row.names = FALSE,sep="\t")
}else{
    cat("No clusters passed FDR<=0.05\n")
}