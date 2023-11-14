###############################################################################
# TITLE: CODE festschrift Stull et al
# PURPOSE: initial data analysis
# LAST UPDATED ON 2023/09/02 (JTH)
# UPDATES: 
# 1. 2023/09/06:  added in ida  code from course
# NOTES: 
################################################################################
# START WITH A CLEAN WORKING DIRECTORY
rm(list=ls())
################################################################################
# OPEN USEFUL LIBRARIES 
################################################################################
library(mice)
library(naniar)
library(kableExtra) #display tables
library(forcats)  #in tidyverse, to work with categorical variables
library(GGally) #includes ggpairs
library(pheatmap) #to plot heatmaps
library(gt) #for tables that use gt
library(dendextend) #coloring dendrograms
library(viridis) #colors
library(RColorBrewer) #colors
library(plotly)
library(reshape2)
library(gridExtra)
library(descr) #set descr package option to skip plots by default when using freq() & crossbar().
  options(descr.plot=FALSE)
################################################################################
# set seed for repeatability
set.seed(1234)
################################################################################
# load data
data <- read.csv("data.csv")
################################################################################
# FUNCTIONS for various
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0
trim <- function (x) gsub("^\\s+|\\s+$", "", x) # cleans up string data
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}
#################################
# color - can be used inline to change the color of the text, default is red
color <- function(x, color="red") {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
            x)
  } else x
}
################################################################################
#initial data analysis
# descriptive stats
des.stat<-psych::describeBy(data,data$Sex)
################################################################################
# Item (per variable) missingness
tab1 <- data  %>%
  dplyr::select(Sex, ANS,INA,IOB,MT,NAS,NAW,NBC,NBS,NO,OBS,PBD,PZT,GOL,NOL,BNL,
                BBH,XCB,WFB,ZYB,AUB,ASB,BPL,NPH,NLH,JUB,MDH,
                NLB,OBH,OBB,DKB,FMB,EKB,FRC,PAC,OCC,FOL,FOB) %>%    miss_var_summary() %>%
  gt::gt() %>%
  gt::cols_label(
    variable = "Variable",
    n_miss = "Missing (count)",
    pct_miss = "Missing (%)"
  ) %>%
  gt::fmt_number(
    columns = c(pct_miss),
    decimals = 2
  )
################################################################################
# to get scrollable gt tables
tab1 %>% gt::tab_options(container.height = px(2000))
################################################################################
#visualizations (mms and sex)

a<-ggplot(data, aes(x=ANS)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Anterior Nasal Spine")
b<-ggplot(data, aes(x=INA)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Inferior Nasal Aperture")
c<-ggplot(data, aes(x=IOB)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Interorbital Breadth")
d<-ggplot(data, aes(x=MT)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Malar Tubercle")
e<-ggplot(data, aes(x=NAS)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Nasal Aperture Shape")
f<-ggplot(data, aes(x=NAW)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Nasal Aperture Width")
g<-ggplot(data, aes(x=NBC)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Nasal Bone Contour")
h<-ggplot(data, aes(x=NBS)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Nasal Bone Shape")
i<-ggplot(data, aes(x=NO)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Nasal Overgrowth")
j<-ggplot(data, aes(x=OBS)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Orbit Shape")
k<-ggplot(data, aes(x=PBD)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Postbregmatic Depression")
l<-ggplot(data, aes(x=PZT)) + geom_bar(fill="white", color="black") + theme_bw()+labs(x="Posterior Zygomatic Tubercle")

grid.arrange(a,b,c,d,e,f,g,h,i,j,k,l,nrow=4)
# craniometric - bleh...don't care for this
 data.1<-data[15:25]
 data.2<-data[26:36]
 data.3<-data[36:39]
 pairs(data.1)
 pairs(data.2)
 pairs(data.3)
################################################################################
#set  breaks for heatmaps, to make colors comparable
Breaks <- seq(-1, 1, by=.01)
# select a color palette
my_palette <- c( colorRampPalette( rev(RColorBrewer::brewer.pal(n = 11,
              name ="Spectral")))(length(Breaks)-2) , "grey80", "grey80")
################################################################################
# heat maps for variable exploration
which.variables <- c(3:39)

pheatmap(cor(data[,which.variables], use="p", method="p"), cluster_cols=TRUE,
         cluster_rows=TRUE, display_numbers = TRUE,
         main="Correlation between variables", breaks=Breaks,
         #labels_row = paste("Wave", 1:7),  labels_col = paste("Wave", 1:7),
         color= my_palette)
################################################################################
# hierarchical clustering
plot(hclust(as.dist(1-cor(data[3:39], method="p"))), xlab="Distance: 1-correlation")

# correlation/covariation
my.mat <- cor(data[3:39], use="p", method="s")
diag(my.mat) <- apply(data[3:39], 2, sd, na.rm=TRUE)
my.mat[lower.tri(my.mat)] <- cov(data[3:39], use="p")[lower.tri(cov(data[3:39], use="p"))]
write.csv(my.mat,"cor_sd_covary.csv")
################################################################################
# table the above
kable(my.mat, digits=2, caption ="Correlation/SD/covariances between variables") %>% kable_styling()
################################################################################
# heatmap
pheatmap(data[,c(3:39)], scale="row", color= my_palette,  show_rownames = TRUE)
clust.compl <- hclust(dist(data[,3:39]), method="complete")
plot(clust.compl)
table(cutree(clust.compl, k=150))
nodePar <- list(lab.cex = 0.6, pch = c(NA, 5),
                cex = 0.7, col = my_palette)
################################################################################
# dendrogram
#save the results in the format dealt with by ibrary(dendextend)
dend <- as.dendrogram(hclust(dist(data[,3:39]), method="complete"),
                      nodePar = nodePar, leaflab = "none")
dend1 <- color_branches(dend, k=150)
plot(dend1, type = "rectangle",nodePar = nodePar, leaflab = "none",  ylab= "Height",
     xlab = "Sex")
par(mfrow=c(1,1))
par(mar = c(12,4,1,1))
################################################################################
# variable exploration with medians, IQR
# color selection
cols3 <- my_palette
################################################################################
# set cuts for variables
ANS <- (cols3[cut(data$ANS, breaks=quantile(data$ANS, c(0,.25,1)), include.lowest=TRUE)])
INA <- (cols3[cut(data$INA, breaks=quantile(data$INA, c(0,.5,1)), include.lowest=TRUE)])
IOB <- (cols3[cut(data$IOB, breaks=quantile(data$IOB, c(0,.5,1)), include.lowest=TRUE)])
MT <- (cols3[cut(data$MT, breaks=quantile(data$MT, c(0,.5,1)), include.lowest=TRUE)])
NAS <- (cols3[cut(data$NAS, breaks=quantile(data$NAS, c(0,1)), include.lowest=TRUE)])
NAW <- (cols3[cut(data$NAW, breaks=quantile(data$NAW, c(0,.5,1)), include.lowest=TRUE)])
NBC <- (cols3[cut(data$NBC, breaks=quantile(data$NBC, c(0,.5,1)), include.lowest=TRUE)])
NBS <- (cols3[cut(data$NBS, breaks=quantile(data$NBS, c(0,.5,1)), include.lowest=TRUE)])
NO <- (cols3[cut(data$NO, breaks=quantile(data$NO, c(0,1)), include.lowest=TRUE)])
OBS <- (cols3[cut(data$OBS, breaks=quantile(data$OBS, c(0,1)), include.lowest=TRUE)])
PBD <- (cols3[cut(data$PBD, breaks=quantile(data$PBD, c(0,1)), include.lowest=TRUE)])
PZT <- (cols3[cut(data$PZT, breaks=quantile(data$PZT, c(0,.5,1)), include.lowest=TRUE)])
Sex <- brewer.pal(11, 'Spectral') [factor(data$Sex)]
plot(dend1, nodePar = nodePar, leaflab = "none", ylab="Height")
#mms medians
my.medians.max <- data %>% group_by(Group ) %>% summarise(ANS=median(ANS), INA=median(INA), IOB=median(IOB),
                                                      MT=median(MT), NAW=median(NAW),NBC=median(NBC),
                                                       NO=median(NO),PBD=median(PBD),PZT=median(PZT), ZS=median(ZS))
my.medians.max.all <- data %>% summarise(ANS=median(ANS), INA=median(INA), IOB=median(IOB),
                                         MT=median(MT), NAW=median(NAW),NBC=median(NBC),
                                         NO=median(NO),PBD=median(PBD),PZT=median(PZT), ZS=median(ZS))
my.IQR.max.all <- data %>% summarise(ANS=median(ANS), INA=median(INA), IOB=median(IOB),
                                              MT=median(MT), NAW=median(NAW),NBC=median(NBC),
                                              NO=median(NO),PBD=median(PBD),PZT=median(PZT), ZS=median(ZS))
################################################################################

# transforming in a data frame with types as names
my.medians.max <- as.data.frame(my.medians.max)
dimnames(my.medians.max)[[1]] <- my.medians.max[,1]
my.medians.max <- my.medians.max[,-1]
pheatmap(my.medians.max,main="MMS Traits - Median Distribution by Population Affinity")
#mms means
my.means.max <- data %>% group_by(Group) %>% summarise(ANS=mean(ANS), INA=mean(INA), IOB=mean(IOB),
                                                       MT=mean(MT), NAW=mean(NAW),NBC=mean(NBC),
                                                       NO=mean(NO),PBD=mean(PBD),PZT=mean(PZT), SPS=mean(SPS),ZS=mean(ZS))

my.means.max.all <- data %>% summarise(ANS=mean(ANS), INA=mean(INA), IOB=mean(IOB),
                                                MT=mean(MT), NAW=mean(NAW),NBC=mean(NBC),
                                                NO=mean(NO),PBD=mean(PBD),PZT=mean(PZT),SPS=mean(SPS), ZS=mean(ZS))

my.IQR.max.all <- data %>% summarise(ANS=mean(ANS), INA=mean(INA), IOB=mean(IOB),
                                              MT=mean(MT), NAW=mean(NAW),NBC=mean(NBC),
                                              NO=mean(NO),PBD=mean(PBD),PZT=mean(PZT),SPS=mean(SPS), ZS=mean(ZS))

################################################################################
# transforming in a data frame with types as names
my.means.max <- as.data.frame(my.means.max)
dimnames(my.means.max)[[1]] <- my.means.max[,1]
my.means.max <- my.means.max[,-1]
pheatmap(my.means.max,main="MMS Traits - Distribution by Population Affinity",
         angle_col = 0,  # Specify the rotation angle (in degrees)
         fontsize_col = 10,
         fontsize_row = 12# Adjust the font size as needed
)
################################################################################
# End of program
################################################################################