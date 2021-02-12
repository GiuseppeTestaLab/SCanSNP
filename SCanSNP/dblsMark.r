
#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)

scoreTablePath= args[1]
outDir = args[2]




library(ggplot2)
library(mixtools)
library(dplyr)
library(robustbase)
library(cluster)
library(gridExtra)
library(grid)


##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
#---Functions
#---Functions
#---Functions
#---Functions
#---Functions
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
IterateFitting <- function(DBLmetricsDF) {
  FittedIterations <- data.frame()
  #Iterating over ranger of Second ID / Total droplet Private reads contributions
  #To identify peaks in Doublets detection and best thresholding
    for (SecondIDfraction in seq(0.10,.60, .01)){
    ReadSignals <- data.frame(row.names = rownames(DBLmetricsDF))
    Outliers="Trimmed"
    ID_specificThresholds <- list()
    #Defining "negative" droplets for each ID
    for (id in unique(DBLmetricsDF$FirstID)){
      IDNoise <- DBLmetricsDF[(DBLmetricsDF$FirstID != id) & !((DBLmetricsDF$SecondID == id) & (DBLmetricsDF[paste0("ID_",id)] >= rowSums(DBLmetricsDF[paste0("ID_",unique(DBLmetricsDF$FirstID))])*SecondIDfraction )),paste0("ID_",id)]
      #We try to trim top 5% of droplet signal before fitting negbin in order to remove outliers
      CroppedNoise <- IDNoise[IDNoise < quantile(IDNoise, probs =.95)[[1]]]
      if (sum(CroppedNoise) == 0) {
        #In case trimming fails (result in all 0 values) we fall back on untrimmed values
        #storing the result in $Outliers variable for subsequent assessment
        print(".5 outliers trimming was not possible")
        Outliers="NotTrimmed"
        CroppedNoise <- IDNoise
      }
        distFitt <- fitdistrplus::fitdist(CroppedNoise, distr = "nbinom")
        Threshold <- quantile(distFitt, probs = .99)
        ID_specificThresholds[[paste0(id,"_R.Threshold")]] <- Threshold$quantiles[[1]]
        Positiveid <- rownames(DBLmetricsDF[DBLmetricsDF[,paste0("ID_",id)] > Threshold$quantiles[[1]],])
        ReadSignals[,id] <- ifelse(rownames(ReadSignals) %in% Positiveid, 1, 0)
    }
    #Subsetting for double (or more) positive droplets
    MultiSIg <- ReadSignals[rowSums(ReadSignals) > 1,]

    # We store the number of droplets in which each ID appear as multiplet component
    idSpecComp <- apply(MultiSIg, 1, function(x) colnames(MultiSIg)[which(x == 1)])
    idSpecComp <- table(unlist(apply(MultiSIg, 1, function(x) colnames(MultiSIg)[which(x == 1)])))

    tempDF <- data.frame("SecondID_fraction"=SecondIDfraction, "DBLs"=rep(length(which(rowSums(ReadSignals) > 1)), ))

    for (id in names(idSpecComp)){
      tempDF[paste0("DBLloss_",id)]  <- idSpecComp[[id]]
    }

    tempDF <-cbind(tempDF, as.data.frame(ID_specificThresholds))

    FittedIterations <- data.table::.rbind.data.table(FittedIterations, tempDF, fill = T)

    print(paste0("Completed ",SecondIDfraction, " ","ID2 fraction iteration"))
    }
    FittedIterations <- as.data.frame(FittedIterations)

    retList <- list("FittedIterations" = FittedIterations, "Outliers" = Outliers)
    return(retList)
}


ID2_RatioSelect <- function(FittedIterations,Outliers ){
  #For each iteration we store the number of lost doublets with respect of previous one to spot the highest "jump"
  FittedIterationsLoss <- FittedIterations
  FittedIterationsLoss$DBLsLoss <- abs(c(0,(FittedIterationsLoss[-1,"DBLs"] - FittedIterationsLoss[-nrow(FittedIterations),"DBLs"])))


  #If 5% outliers trimming failed we ignore the first peak because it will likely be due to non trimmed outliers retention
  if (Outliers == "NotTrimmed"){
  FittedIterationsLoss <-  FittedIterationsLoss[3:length(FittedIterationsLoss$DBLsLoss),]
  ID2_Ratio = FittedIterationsLoss[which.max(FittedIterationsLoss$DBLsLoss)-1,"SecondID_fraction"]

  } else {
  ID2_Ratio = FittedIterationsLoss[which.max(FittedIterationsLoss$DBLsLoss)-1,"SecondID_fraction"]
  }

  return(ID2_Ratio)
}


DBLsMark <- function(DBLmetricsDF, ID2_Ratio) {
  #Perform actual DBLs marking based on selected ID2 ratio negbin fitting
  ReadSignals <- data.frame(row.names = rownames(DBLmetricsDF))
  for (id in unique(DBLmetricsDF$FirstID)){
    IDNoise <- DBLmetricsDF[(DBLmetricsDF$FirstID != id) & !((DBLmetricsDF$SecondID == id) & (DBLmetricsDF[paste0("ID_",id)] >= rowSums(DBLmetricsDF[paste0("ID_",unique(DBLmetricsDF$FirstID))])*ID2_Ratio )),paste0("ID_",id)]
    #We try to trim top 5% of droplet signal before fitting negbin in order to remove outliers
    CroppedNoise <- IDNoise[IDNoise < quantile(IDNoise, probs =.95)[[1]]]
    if (sum(CroppedNoise) == 0) {
      #In case trimming fails (result in all 0 values) we fall back on untrimmed values
      #storing the result in $Outliers variable for subsequent assessment
      print(".5 outliers trimming was not possible")
      Outliers="NotTrimmed"
      CroppedNoise <- IDNoise
    }
      distFitt <- fitdistrplus::fitdist(CroppedNoise, distr = "nbinom")
      Threshold <- quantile(distFitt, probs = .99)
      Positiveid <- rownames(DBLmetricsDF[DBLmetricsDF[,paste0("ID_",id)] > Threshold$quantiles[[1]],])
      ReadSignals[,id] <- ifelse(rownames(ReadSignals) %in% Positiveid, 1, 0)
  }
  #Subsetting for double (or more) positive droplets
  MultiSIg <- ReadSignals[rowSums(ReadSignals) > 1,]
  DBLbarcodes <- names(rowSums(ReadSignals)[rowSums(ReadSignals) > 1])
  DBLsMarked <- data.frame(row.names = rownames(DBLmetricsDF), DBLstatus = ifelse(rownames(DBLmetricsDF) %in% DBLbarcodes, "Doublet", "Singlet"))

  return(DBLsMarked)
}



PlotDBLsLoss <- function(FittedIterations, DBLmetricsDF){
  IDSpecDBLloss <- abs(FittedIterations[-1,paste0("DBLloss_",unique(DBLmetricsDF$FirstID))] - FittedIterations[-nrow(FittedIterations),paste0("DBLloss_",unique(DBLmetricsDF$FirstID))])
  #Add a 0 first row
  firstrow <- as.data.frame(matrix(nrow = 1,ncol = ncol(IDSpecDBLloss)))
  colnames(firstrow) <- colnames(IDSpecDBLloss)
  firstrow[ is.na(firstrow)] = 0
  IDSpecDBLloss <- rbind(firstrow,IDSpecDBLloss )
  #Parallely create same dataframe but with ID losses normalized by total ID of same kind
  IDsScaling <- table(DBLmetricsDF$FirstID)
  names(IDsScaling) <- paste0("DBLloss_", names(IDsScaling))
  IDSpecDBLlossNorm <- as.data.frame(t(t(IDSpecDBLloss)/as.vector(IDsScaling[colnames(IDSpecDBLloss)])))


  #Re-add info about SeconID fraction
  IDSpecDBLloss$SecondID_fraction <- as.factor(FittedIterations$SecondID_fraction)
  IDSpecDBLlossNorm$SecondID_fraction <- as.factor(FittedIterations$SecondID_fraction)
  #And Total DBL loss
  IDSpecDBLloss$OverallDBLloss <- abs(c(0,(FittedIterations[-1,"DBLs"] - FittedIterations[-nrow(FittedIterations),"DBLs"])))

  #NotNormPlot
  MeltedLoss <- reshape2::melt(IDSpecDBLloss)
  colnames(MeltedLoss) <- c("SecondID_fraction","Lost_ID","Lost_DBLs")
  DeltaLossPlot <- ggplot(data=MeltedLoss, aes(x=SecondID_fraction, y=Lost_DBLs, fill=Lost_ID)) +
  geom_bar(stat="identity", position=position_dodge())+
  ylab("Lost DBLs")+
  theme(axis.text.x = element_text(angle = 45, size = 1))+
    theme_bw()+
    theme(plot.title = element_text(size = 12),
    legend.title=element_text(size=7),
    legend.text=element_text(size=5),
    legend.key.size = unit(.5,"line"),
    panel.grid.major = element_blank(),
    axis.title=element_text(size=5),
    axis.text.y = element_text(size = 4),
    axis.text.x = element_text(angle = 45, size = 4))

    #NormPlot
    MeltedLossNorm <- reshape2::melt(IDSpecDBLlossNorm)
    colnames(MeltedLossNorm) <- c("SecondID_fraction","Lost_ID","Lost_DBLs")
    DeltaLossPlotNorm <- ggplot(data=MeltedLossNorm, aes(x=SecondID_fraction, y=Lost_DBLs, fill=Lost_ID)) +
    geom_bar(stat="identity", position=position_dodge())+
    ylab("Lost DBLs normalized by FirstID counts")+
    theme(axis.text.x = element_text(angle = 45, size = 1))+
      theme_bw()+
      theme(plot.title = element_text(size = 12),
      legend.title=element_text(size=7),
      legend.text=element_text(size=5),
      panel.grid.major = element_blank(),
      legend.key.size = unit(.5,"line"),
      axis.title=element_text(size=5),
      axis.text.y = element_text(size = 4),
      axis.text.x = element_text(angle = 45, size = 4))

  GridDeltaPlot <- gridExtra::arrangeGrob(DeltaLossPlot, DeltaLossPlotNorm, nrow = 2,top = textGrob("Delta DBLs loss per ID2 fraction",gp=gpar(fontsize=7)))

  ggsave(GridDeltaPlot, file=paste0(outDir, "/DeltaLossPlot.png") , width = 21, height = 10, units = "cm")

  #Prepare linespecific DBLs loss
  LinesLossDF <- FittedIterations[paste0("DBLloss_",unique(DBLmetricsDF$FirstID))]
  #Re-add infos
  LinesLossDF$SecondID_fraction <- as.factor(FittedIterations$SecondID_fraction)
  LinesLossDF$OverallDBLloss <- FittedIterations$DBLs

  MeltedLinesLoss <- reshape2::melt(LinesLossDF)
  colnames(MeltedLinesLoss) <- c("SecondID_fraction","Identity","Doublets_Containing_ID")

  LinesLossPlot <- ggplot(data=MeltedLinesLoss, aes(x=SecondID_fraction, y=Doublets_Containing_ID, group=Identity)) +
  geom_line(aes(color=Identity))+
  geom_point(aes(color=Identity), size = 0.05)+
  theme(axis.text.x = element_text(angle = 45, size = 1))+
    theme_bw()+
    geom_segment(aes(x = as.factor(ID2_Ratio), y = min(Doublets_Containing_ID), xend = as.factor(ID2_Ratio),
    yend = max(MeltedLinesLoss$Doublets_Containing_ID)),
        size = 0.5, color = "red",linetype="dotted")+
    theme(plot.title = element_text(size = 12),
    legend.title=element_text(size=7),
    legend.text=element_text(size=5),
    legend.key.size = unit(.5,"line"),
    axis.title=element_text(size=5),
    axis.text.y = element_text(size = 4),
    axis.text.x = element_text(angle = 45, size = 4))

  ggsave(LinesLossPlot, file=paste0(outDir, "/LinesLossPlot.png") , width = 21, height = 10, units = "cm")

}




AddPseudoCounts <- function(scoreTableSNGs) {

  NoiseDF <- data.frame()

  NoisedscoreTableSNGs <- scoreTableSNGs

  TempNoisedscoreTableSNGs <- NoisedscoreTableSNGs
  TempNoisedscoreTableSNGs$Paste <- paste0(NoisedscoreTableSNGs$FirstID,"@",NoisedscoreTableSNGs$SecondID)
  for (idPair in unique(TempNoisedscoreTableSNGs$Paste)){

    AllIDs <- unique(c(unique(TempNoisedscoreTableSNGs$FirstID),unique(TempNoisedscoreTableSNGs$SecondID)))
    NotDBLids <- AllIDs[! AllIDs %in% c(unlist(strsplit(idPair, "@"))[1]  ,unlist(strsplit(idPair, "@"))[2])]

    noiseSlice <- TempNoisedscoreTableSNGs[TempNoisedscoreTableSNGs$Paste == idPair,paste0("ID_",NotDBLids)]
    if( length(AllIDs == 3)){
      noiseSlice <- as.data.frame(noiseSlice, stringsAsFactors=F)
      colnames(noiseSlice) <- paste0("ID_",NotDBLids)
    }
    NoiseDF <- bind_rows(NoiseDF,noiseSlice)

  }


  NoiseDFNorm <- NoiseDF[, !names(NoiseDF) %in% c("FirstID_Counts")]/TempNoisedscoreTableSNGs$DBL_FirstID_Score
  NoiseDFNorm[is.na(NoiseDFNorm)] <- 0
  NoiseDFNorm <- NoiseDFNorm[!is.infinite(rowSums(NoiseDFNorm)),]
  NoisedscoreTableSNGs$NoiseFactor_ID1 <- paste0("ID_",NoisedscoreTableSNGs$FirstID)
  NoisedscoreTableSNGs$NoiseFactor_ID2 <- paste0("ID_",NoisedscoreTableSNGs$SecondID)



  NoiseFactor <- colMeans(NoiseDFNorm)

  NoiseIDs <- unique(NoisedscoreTableSNGs$NoiseFactor_ID1)
  NoisedscoreTableSNGs$NoiseFactor_ID1 <- as.numeric(c(as.vector(NoiseFactor[NoiseIDs]), NoisedscoreTableSNGs$NoiseFactor_ID1)[match(NoisedscoreTableSNGs$NoiseFactor_ID1, c(NoiseIDs, NoisedscoreTableSNGs$NoiseFactor_ID1))])

  NoiseIDsII <- unique(NoisedscoreTableSNGs$NoiseFactor_ID2)
  NoisedscoreTableSNGs$NoiseFactor_ID2 <- as.numeric(c(as.vector(NoiseFactor[NoiseIDsII]), NoisedscoreTableSNGs$NoiseFactor_ID2)[match(NoisedscoreTableSNGs$NoiseFactor_ID2, c(NoiseIDsII, NoisedscoreTableSNGs$NoiseFactor_ID2))])


  NoisedscoreTableSNGs$DBL_FirstIDScore_Noised <- (NoisedscoreTableSNGs$NoiseFactor_ID1*NoisedscoreTableSNGs$DBL_FirstID_Score)+NoisedscoreTableSNGs$DBL_FirstID_Score
  NoisedscoreTableSNGs$DBL_SecondIDScore_Noised <- (NoisedscoreTableSNGs$NoiseFactor_ID2*NoisedscoreTableSNGs$DBL_FirstID_Score)+NoisedscoreTableSNGs$DBL_SecondID_Score

  return(NoisedscoreTableSNGs)
}







ModelFitting <- function(NoisedScoreTableSNGs,outDir) {
  set.seed(43)


  LogTable <- NoisedScoreTableSNGs[complete.cases(NoisedScoreTableSNGs[c("DBL_FirstIDScore_Noised","DBL_SecondIDScore_Noised")]),c("DBL_FirstIDScore_Noised","DBL_SecondIDScore_Noised") ]
  LogTable$logFC <- log(LogTable$DBL_FirstIDScore_Noised/LogTable$DBL_SecondIDScore_Noised)
  LogTable <- LogTable[complete.cases(LogTable),]

  #We trim top and bottom 1% of the values assuming they are arising by a) Adding Pcounts to 0 counts background when firstID is highly represented or b) Absence in counts of the firstID (Should be very rare cases)
  # "a" Cases can be safely considered good quality, , while we consider "b" cases as "Low quality droplets"
  GoodQualBarcodes <- rownames(LogTable[LogTable$logFC >= quantile(LogTable$logFC, .99),])
  LowQualBarcodes <- rownames(LogTable[LogTable$logFC <= quantile(LogTable$logFC, .01),])
  LogTable <- LogTable[LogTable$logFC < quantile(LogTable$logFC, .99) & (LogTable$logFC > quantile(LogTable$logFC, .01)),]



  plotDsist = function(x, mu=0, s=1, lambda=1){lambda*dnorm(x, mean=mu, sd=s)}

  models <- list()
  fittings <- list()
  MuRatiRange <- seq(0.10, 0.31, 0.01)
  for (MuRatio in seq(length(MuRatiRange))){
    LocalRatio <- MuRatiRange[MuRatio]
    models[[MuRatio]] <- normalmixEM(LogTable$logFC , k = 2, mean.constr = c(paste0(LocalRatio, "a"),"a"), epsilon = 1e-08, maxit = 1000, maxrestarts=100, arbmean = TRUE, arbvar = TRUE)
    fittings[[MuRatio]] <- models[[MuRatio]]$loglik
  }


  bestFit <- which.max(unlist(fittings))
  FittedModel <- models[[bestFit]]

  FittedMixture <- ggplot(data.frame(x=LogTable$LogFC)) +
    geom_histogram(aes(x=LogTable$logFC,y=..density..),fill="white",color="black", binwidth = .028) +
    stat_function(fun=plotDsist,
                  args=list(mu=FittedModel$mu[1],
                            s=FittedModel$sigma[1],
                            lambda=FittedModel$lambda[1]),
                  fill="#FF3030",geom="polygon",alpha = .4)+
    stat_function(fun=plotDsist,
                  args=list(mu=FittedModel$mu[2],
                            s=FittedModel$sigma[2],
                            lambda=FittedModel$lambda[2]),
                  fill="#00FA9A",geom="polygon", alpha = .4) +
    geom_density(aes(x=LogTable$logFC,y=..density..), linetype="dashed", size=.8)+
    labs(title="FirstID-SecondID LogFC mixture distributions",x="FirstID-SecondID logFC", y = "Density")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5))

  rownames(FittedModel$posterior) <- rownames(LogTable)


  ggsave(FittedMixture, file=paste0(outDir, "/logFC_FittedMixture.png") , width = 14, height = 10, units = "cm")

  ReturnList <- list("FittedModel" = FittedModel, "GoodQualBarcodes" = GoodQualBarcodes, "LowQualBarcodes"= LowQualBarcodes)
  return(ReturnList)


}




##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
#---MAIN
#---MAIN
#---MAIN
#---MAIN
#---MAIN
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################



DBLmetricsDF <- read.delim(scoreTablePath,header = T,stringsAsFactors=FALSE, row.names = 1,check.names = FALSE)

#Threshold iterations
IterateFittingRetList <- IterateFitting(DBLmetricsDF)
FittedIterations <- IterateFittingRetList$FittedIterations
Outliers <- IterateFittingRetList$Outliers





ID2_Ratio <- ID2_RatioSelect(FittedIterations,Outliers )
DBLsMarked <- DBLsMark(DBLmetricsDF, ID2_Ratio)

#DBLlossPlots
PlotDBLsLoss(FittedIterations, DBLmetricsDF)

#Removing spotted doublets before GMM fitting
DBLsBarcodes<- rownames(DBLsMarked)[(DBLsMarked$DBLstatus == "Doublet")]
scoreTableSNGs <- DBLmetricsDF[! rownames(DBLmetricsDF) %in% DBLsBarcodes,]


#Adding PScounts that mimic ambient RNA (inspired by https://rdrr.io/github/MarioniLab/DropletUtils/man/hashedDrops.html)
NoisedScoreTableSNGs <- AddPseudoCounts(scoreTableSNGs)

#GMM fitting
ModelFittingRetList <- ModelFitting(NoisedScoreTableSNGs,outDir)
FittedModel <- ModelFittingRetList$FittedModel


SafeGoodQualBarcodes <- ModelFittingRetList$GoodQualBarcodes
if (length(SafeGoodQualBarcodes) > 0) {
  SafeGoodQualBarcodes <- data.frame(row.names = SafeGoodQualBarcodes, Qual = rep("GoodQuality", length(SafeGoodQualBarcodes)))
} else {SafeGoodQualBarcodes <- data.frame()}


SafeLowQualBarcodes <- ModelFittingRetList$LowQualBarcodes
if (length(SafeLowQualBarcodes) > 0 ){
  SafeLowQualBarcodes <- data.frame(row.names = SafeLowQualBarcodes, Qual = rep("LowQuality", length(SafeLowQualBarcodes)))
} else {SafeLowQualBarcodes <- data.frame()}





#Retrieving Quality Label from fitted mixture
LowQualComp <- paste0("comp.",which.min(FittedModel$mu))
GoodQualComp <- paste0("comp.",which.max(FittedModel$mu))
Posteriors <- as.data.frame(FittedModel$posterior)

#LowQual component
LowQual <- data.frame(row.names  = rownames(Posteriors[Posteriors[,LowQualComp] - Posteriors[,GoodQualComp] > 0,]),
                      Qual = rep("LowQuality",length(rownames(Posteriors[Posteriors[,LowQualComp] - Posteriors[,GoodQualComp] > 0,])), stringsAsFactors = F))
LowQual <- rbind(LowQual, SafeLowQualBarcodes)

#GoodQual component
GoodQual <- data.frame(row.names  = rownames(Posteriors[Posteriors[,LowQualComp] - Posteriors[,GoodQualComp] < 0,]),
                       Qual = rep("GoodQuality",length(rownames(Posteriors[Posteriors[,LowQualComp] - Posteriors[,GoodQualComp] < 0,])), stringsAsFactors = F))
GoodQual <- rbind(GoodQual, SafeGoodQualBarcodes)

#Re-add DBLs that were removed before computing quality
FinalDBLsDF <- rbind(LowQual,GoodQual, data.frame(Qual = rep("Doublet",length(DBLsBarcodes)), row.names = DBLsBarcodes, stringsAsFactors = F))
FinalDBLsDF$Qual <- as.vector(FinalDBLsDF$Qual)

#Building the final DF
FinalDBLsDF$FirstID <- DBLmetricsDF[rownames(FinalDBLsDF),"FirstID"]
FinalDBLsDF$SecondID <- DBLmetricsDF[rownames(FinalDBLsDF),"SecondID"]
FinalDBLsDF$Type <- as.vector(DBLsMarked[rownames(FinalDBLsDF),"DBLstatus"])

#ID-QUAL column
FinalDBLsDF$ID_Qual <- ifelse(FinalDBLsDF$Qual == "GoodQuality", FinalDBLsDF$FirstID, FinalDBLsDF$Qual)
#ID-TYPE column
FinalDBLsDF$ID <- ifelse(FinalDBLsDF$Type == "Singlet", FinalDBLsDF$FirstID, FinalDBLsDF$Type)
#Barcode column
FinalDBLsDF$Barcode <- rownames(FinalDBLsDF)


#Arrangin cols
FinalDBLsDF <- FinalDBLsDF[,c("Barcode","ID","Type","Qual","FirstID","SecondID","ID_Qual")]


write.table(FinalDBLsDF, paste0(outDir, "/doubletsMarked.tsv"), sep = "\t", row.names = F, col.names = T,quote = F)
