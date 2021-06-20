
#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)

scoreTablePath= args[1]
outDir = args[2]


# .libPaths("/usr/local/lib/R/site-library")
.libPaths("/usr/lib/R/library")
# .libPaths("/usr/local/lib/R/site-library")


library(ggplot2)
library(mixtools)
library(dplyr)
library(robustbase)
library(cluster)
library(gridExtra)
library(grid)





ModelFitting <- function(NoisedScoreTableSNGs,outDir, HighOutlierThreshold, LowOutlierThreshold) {
  set.seed(43)


   LogTable <- NoisedScoreTableSNGs[complete.cases(NoisedScoreTableSNGs[c("Noised_FirstID_Score","Noised_SecondID_Score")]),c("Noised_FirstID_Score","Noised_SecondID_Score") ]
   LogTable$logFC <- log(LogTable$Noised_FirstID_Score/LogTable$Noised_SecondID_Score)
   LogTable <- LogTable[complete.cases(LogTable),]
   #We trim top and bottom 1% of the values assuming they are arising by a) Adding Pcounts to 0 counts background when firstID is highly represented or b) Absence in counts of the firstID (Should be very rare cases)
   # "a" Cases can be safely considered good quality, , while we consider "b" cases as "Low quality droplets"
   GoodQualBarcodes <- rownames(LogTable[LogTable$logFC >= quantile(LogTable$logFC, HighOutlierThreshold),])
   LowQualBarcodes <- rownames(LogTable[LogTable$logFC <= quantile(LogTable$logFC, LowOutlierThreshold),])
   LogTable <- LogTable[LogTable$logFC < quantile(LogTable$logFC, HighOutlierThreshold) & (LogTable$logFC > quantile(LogTable$logFC, LowOutlierThreshold)),]
   
   
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
   ModelFittingRetList <- list("FittedModel" = FittedModel, "GoodQualBarcodes" = GoodQualBarcodes, "LowQualBarcodes"= LowQualBarcodes)


   SafeGoodQualBarcodes <- ModelFittingRetList$GoodQualBarcodes
if (length(SafeGoodQualBarcodes) > 0) {
  SafeGoodQualBarcodes <- data.frame(row.names = SafeGoodQualBarcodes, Quality= rep("GoodQuality", length(SafeGoodQualBarcodes)))
} else {SafeGoodQualBarcodes <- data.frame()}


SafeLowQualBarcodes <- ModelFittingRetList$LowQualBarcodes
if (length(SafeLowQualBarcodes) > 0 ){
  SafeLowQualBarcodes <- data.frame(row.names = SafeLowQualBarcodes, Quality= rep("LowQuality", length(SafeLowQualBarcodes)))
} else {SafeLowQualBarcodes <- data.frame()}
   


#Retrieving Quality Label from fitted mixture
LowQualComp <- paste0("comp.",which.min(FittedModel$mu))
GoodQualComp <- paste0("comp.",which.max(FittedModel$mu))
Posteriors <- as.data.frame(FittedModel$posterior)

#LowQual component
LowQual <- data.frame(row.names  = rownames(Posteriors[Posteriors[,LowQualComp] - Posteriors[,GoodQualComp] > 0,]),
                      Quality= rep("LowQuality",length(rownames(Posteriors[Posteriors[,LowQualComp] - Posteriors[,GoodQualComp] > 0,])), stringsAsFactors = F))
LowQual <- rbind(LowQual, SafeLowQualBarcodes)

#GoodQual component
GoodQual <- data.frame(row.names  = rownames(Posteriors[Posteriors[,LowQualComp] - Posteriors[,GoodQualComp] < 0,]),
                       Quality= rep("GoodQuality",length(rownames(Posteriors[Posteriors[,LowQualComp] - Posteriors[,GoodQualComp] < 0,])), stringsAsFactors = F))
GoodQual <- rbind(GoodQual, SafeGoodQualBarcodes)


QualDF=rbind(LowQual,GoodQual)

return(QualDF)

}






