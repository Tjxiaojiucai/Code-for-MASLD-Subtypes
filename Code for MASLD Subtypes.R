#########------------ Feature selection-------------------------
library(reshape2)
library(car) 
library(randomForest) 
library(leaps) #best subsets regression
library(glmnet) 
library(caret) 
library(pROC)
library(ggplot2)
library(rms)
library(varSelRF) 
library(InformationValue)
library(reportROC) 
library(Hmisc) 
library(tidyverse)

rm(list=ls())

data <- read.csv("All data statistic_Logistic.csv", header = TRUE, row.name = 1) #Group:“Control" and "Case"
Markers <- read.csv("All data statistic_a_select markers.csv", header = TRUE) 

outputDir <- "C:/HQX/Code/R code/Statisics/All data statistic_Logis_output" 
if(!dir.exists(outputDir)){
  dir.create(outputDir) 
} else {
  print("Dir already exists!")
}

#-----data structure------
categ <- 15 #Number of clinical category variables
e <- 42 #Number of clinical continuous variables

a <- 1 + categ #1 is Group
b <- 1 + categ + e # 2:b is clinical
c <- a + 1 #The begining of continuous variables
d <- b + 1 #The begining of metabolites

#-----set for data pretreatment-------
Trans2 = T #Whether to transform the data; #T or F

HighCorr = T #Whether to remove high correlation; #T or F
r <- 0.8 #correlation coefficient

Train.data = T #Whether to randomize training and testing sets; #T or F
t <- 0.7 #ratio for training and testing sets

Metdata = 4 #Designation of candidate markers
#1:features, Clinical and Metabolites (p<0.05)  
#2:features, only Metabolites (p<0.05);
#3:features, only clinical parameters (p<0.05)  
#4:parameters from "All data statistic_a_select markers.csv"

############################################################
#----------------chisq.test--------------------
Chi_data <- data[1:a] 
if(categ>1){
  name <- names(Chi_data[,2:ncol(Chi_data)])
}else{
  name <- names(Chi_data[2])
}

Chi <- function(x){
  data1 <- factor(x == 0, levels = c(TRUE, FALSE), 
                  labels = c("Missing","Real")) 
  data2 <- table(Chi_data$Group,data1,
                 dnn = c("Group","Metabolite")) 
  data3 <- as.data.frame(chisq.test(data2)$expected) 
  all <- data3[1,1] + data3[1,2] + data3[2,1] + data3[2,2] 
  
  if(all >= 40 && data3[1,1] >= 1 && data3[1,2] >= 1 && data3[2,1] >= 1 && data3[2,2] >= 1){  #N >= 40, T >= 1
    if(data3[1,1] >= 5 && data3[1,2] >= 5 && data3[2,1] >= 5 && data3[2,2] >= 5){ 
      p <- chisq.test(data2, correct = FALSE)$p.value #N >= 40, T >= 5
    }else{ 
      p <- chisq.test(data2, correct = TRUE)$p.value #N >= 40, 1 =< T < 5
    }
  }else{ 
    p <- fisher.test(data2)$p.value #N<40, or T < 1
  }
  
  data.T <- as.data.frame(data2)
  N1 <- data.T[3,3] 
  N1p <- data.T[3,3]/(data.T[1,3]+data.T[3,3])*100 
  N2 <- data.T[4,3]
  N2p <- data.T[4,3]/(data.T[2,3]+data.T[4,3])*100 
  out1 <- c(p, N1, N1p, N2, N2p)
  names(out1) <- c("p", "Case_1N", "Case_1per", "Control_1N", "Control_1per")
  return(out1)
}

Chi_p <- c()
for(i in 2:ncol(Chi_data)){
  x <- Chi_data[i]
  Chi_p <- rbind(Chi_p, Chi(x))
}
Univar_chi <- data.frame('Metabolites' = name[1:length(name)], Chi_p)
write.csv(Univar_chi, paste0(outputDir, "/", "All statistic_Chisquare.csv"), 
          row.names = F, quote = F)

#------T-test & Wilcox test----------------

Ttest_data1 <- cbind(data[,1], data[, c:ncol(data)]) 
colnames(Ttest_data1)[1] <- "Group" 
Ttest_data1$Group <- as.numeric(ifelse(Ttest_data1$Group == "Control",'1','2')) 
Ttest_data2 <- as.data.frame(t(Ttest_data1))
Col1 <- as.data.frame(rownames(Ttest_data2))
Ttest_data <- data.frame(Col1, Ttest_data2)
colnames(Ttest_data)[1] <- "ID" 

pcal2 <- function(datafile,group){
  data.p <- datafile
  groups <- c(1:group)
  comb <- combn(groups,2)
  num <- ncol(comb)
  out <- data.frame(name = data.p[-1,1])
  
  for (k in 1:num){
    m <- comb[,k]
    g1 <- data.p[,data.p[1,] == m[1]]
    group1 <- ncol(g1)
    g2 <- data.p[,data.p[1,] == m[2]]
    group2 <- ncol(g2)
    data.l <- data.frame(g1[-1,],g2[-1,])
    
    tag <- c(rep(0,nrow(data.l)))
    normality1.l <- c(rep(0, nrow(data.l)))
    normality2.l <- c(rep(0, nrow(data.l)))
    homogeneity.l <- c(rep(0, nrow(data.l)))
    nonparametric.l <- c(rep(0, nrow(data.l)))
    pvalue.l <- c(rep(0, nrow(data.l)))
    fdr <- c(rep(0, nrow(data.l)))
    mean1 <- c(rep(0, nrow(data.l)))
    sd1 <- c(rep(0, nrow(data.l)))
    mean2 <- c(rep(0, nrow(data.l)))
    sd2 <- c(rep(0, nrow(data.l)))
    median1 <- c(rep(0, nrow(data.l)))
    median2 <- c(rep(0, nrow(data.l)))
    FC.mea <- c(rep(0, nrow(data.l)))
    FC.med <- c(rep(0, nrow(data.l)))
    #  shapiro.test for group1
    for(i in 1:nrow(data.l)){
      y = shapiro.test(as.numeric(data.l[i,1:group1]))
      normality1.l[i] <- y$p.value
    }
    #   shapiro.test for group2
    for(i in 1:nrow(data.l)){
      y = shapiro.test(as.numeric(data.l[i,(group1+1):(group1+group2)])) 
      normality2.l[i] <- y$p.value
    }
    #-----F-test----
    for(i in 1:nrow(data.l)){
      y = var.test(as.numeric(data.l[i,1:group1]),as.numeric(data.l[i,(group1+1):(group1+group2)])) 
      homogeneity.l[i] <- y$p.value
    }
    #-----p valuue------
    for(i in 1:nrow(data.l)){
      if(normality1.l[i]>0.05&&normality2.l[i]>0.05){ 
        if(homogeneity.l[i] > 0.05){ 
          y = t.test(as.numeric(data.l[i,1:group1]),as.numeric(data.l[i,(group1+1):(group1+group2)]),
                     var.equal = TRUE)
          pvalue.l[i] <- y$p.value
          tag[i] <- "t.test"
        }else{ 
          y = t.test(as.numeric(data.l[i,1:group1]),as.numeric(data.l[i,(group1+1):(group1+group2)]),
                     var.equal = FALSE)
          pvalue.l[i] <- y$p.value
          tag[i] <- "t.test_non"
        }
      }else{  
        y = wilcox.test(as.numeric(data.l[i,1:group1]),as.numeric(data.l[i,(group1+1):(group1+group2)]),
                        exact=FALSE)
        pvalue.l[i] <- y$p.value
        tag[i] <- "wilcox"
      }
    }
    #-----mean and sd----
    for(i in 1:nrow(data.l)){
      mean2[i] <- mean(as.numeric(data.l[i,(group1+1):(group1+group2)])) #Case
      mean1[i] <- mean(as.numeric(data.l[i,1:group1])) #Control
      
      sd2[i] <- sd(as.numeric(data.l[i,(group1+1):(group1+group2)]))
      sd1[i] <- sd(as.numeric(data.l[i,1:group1]))  
      
      median2[i] <- median(as.numeric(data.l[i,(group1+1):(group1+group2)])) #Case
      median1[i] <- median(as.numeric(data.l[i,1:group1])) #Control
      
      FC.mea[i] <- (mean(as.numeric(data.l[i,(group1+1):(group1+group2)]))/mean(as.numeric(data.l[i,1:group1])))
      FC.med[i] <- (median(as.numeric(data.l[i,(group1+1):(group1+group2)]))/median(as.numeric(data.l[i,1:group1])))
    }
    
    #-----p value correction
    fdr = p.adjust(pvalue.l,"BH")
    
    name.l1 <- paste0(m[1],"&",m[2],"_P")
    name.l2 <- paste0(m[1],"&",m[2],"_fdr")
    name.l3 <- paste0(m[1],"&",m[2],"_mean1")
    name.l4 <- paste0(m[1],"&",m[2],"_sd1")
    name.l5 <- paste0(m[1],"&",m[2],"_mean2")
    name.l6 <- paste0(m[1],"&",m[2],"_sd2")    
    name.l7 <- paste0(m[1],"&",m[2],"_FC.mea")
    name.l8 <- paste0(m[1],"&",m[2],"_median1")
    name.l9 <- paste0(m[1],"&",m[2],"_median2")    
    name.l10 <- paste0(m[1],"&",m[2],"_FC.med")
    name.l11 <- paste0(m[1],"&",m[2],"_tag") 
    
    out <- cbind(out,pvalue.l)
    names(out)[length(names(out))] <- name.l1
    
    out <- cbind(out,fdr)
    names(out)[length(names(out))] <- name.l2
    
    out <- cbind(out,mean1)
    names(out)[length(names(out))] <- name.l3
    
    out <- cbind(out,sd1)
    names(out)[length(names(out))] <- name.l4
    
    out <- cbind(out,mean2)
    names(out)[length(names(out))] <- name.l5
    
    out <- cbind(out,sd2)
    names(out)[length(names(out))] <- name.l6   
    
    out <- cbind(out,FC.mea)
    names(out)[length(names(out))] <- name.l7
    
    out <- cbind(out,median1)
    names(out)[length(names(out))] <- name.l8
    
    out <- cbind(out,median2)
    names(out)[length(names(out))] <- name.l9   
    
    out <- cbind(out,FC.med)
    names(out)[length(names(out))] <- name.l10
    
    out <- cbind(out, tag)
    names(out)[length(names(out))] <- name.l11
  }
  return(out)
}

Ttest_result <- pcal2(Ttest_data, group = 2)
write.csv(Ttest_result, paste0(outputDir, "/", file = "All statistic_ttest.csv"), 
          row.names = F, quote = F)

##########################################################
#########------------ Feature selection-------------------------
#-------data pretreatment------------
if(Metdata < 4){
  if(Metdata > 1){
    if(Metdata == 3){
      #3: features, only clinical parameters (p<0.05)
      ttest_Result <- Ttest_result
      ttest_Result_CL <- ttest_Result[1:e,] 
      colnames(ttest_Result_CL)[1:2] <- c("Metabolites", "p") 
      
      Chi_result <- Univar_chi #imput Chisquare results
      colnames(Chi_result)[1:2] <- c("Metabolites", "p") 
      
      Metaname7 <- rbind(Chi_result[1:2], ttest_Result_CL[1:2])
      Metaname8 <- Metaname7[Metaname7[2] < 0.05,] 
      Metaname <- data.frame(Metaname8[1])
      colnames(Metaname)[1] <- "Metabolites" 
    }else{
      #2: features, only Metabolites (p<0.05)
      ttest_Result <- Ttest_result
      f <- e+1 
      Metaname4 <- ttest_Result[f:nrow(ttest_Result),] 
      colnames(Metaname4)[1:2] <- c("Metabolites", "p") 
      
      Metaname5 <- Metaname4[Metaname4[2] < 0.05,] 
      Metaname <- data.frame(Metaname5[1])
      colnames(Metaname)[1] <- "Metabolites" 
    }
    
  }else{
    #1: features, Clinical and Metabolites (p<0.05)
    ttest_Result <- Ttest_result #imput ttest results
    colnames(ttest_Result)[1:2] <- c("Metabolites", "p") 
    Chi_result <- Univar_chi #imput Chisquare results
    colnames(Chi_result)[1:2] <- c("Metabolites", "p") 
    Metaname1 <- rbind(Chi_result[1:2], ttest_Result[1:2])
    Metaname2 <- Metaname1[Metaname1[2] < 0.05,] 
    Metaname <- data.frame(Metaname2[1])
    colnames(Metaname)[1] <- "Metabolites" 
  }
  
}else{
  #4: parameters from "All data statistic_a_select markers.csv"
  Metaname <- Markers
}

cols_remain <- factor(Metaname[,1]) 
data2 <- data[, colnames(data) %in% cols_remain]

#--------data transform--------
if(Trans2 == T){ 
  int <- function(x){
    qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x)))
  }
  data2 <- apply(data2, 2, int)#inverse normal transformation
  data2 <- scale(data2) #Scale
}else{ 
  data2 <- data2 
}

All_data <- data.frame(data[1], data2)
All_data$Group <- factor(All_data$Group,levels = c('Control','Case'))

#-------train and test data------------
if(Train.data == T){ 
  set.seed(1234)
  train <- sample(nrow(All_data), t*nrow(All_data)) #70% train.data
  train.data <- All_data[train,]
  test.data <- All_data[-train,] #30% test.data
}else{ 
  train.data <- All_data 
  test.data <- All_data 
}

#remove high correlation--------------------------------------------------
if(HighCorr == T){ 
  highCorr <- findCorrelation(cor(train.data[,-1], use = 'na.or.complete'), r) 
  if(length(highCorr) > 0){               
    train.data <- train.data[,-(highCorr+1)]
  }                                        
  write.csv(train.data, paste0(outputDir, "/", "All statistic_Highcorr_0.8.csv"))
  #train.data$Group <- factor(train.data$Group,levels = c('Control','Case'))
}else{ 
  train.data <- train.data 
}

#--------------------LASSO---------------------
x <- as.matrix(train.data[, 2:ncol(train.data)])
y <- train.data[, 1]
set.seed(3)
lasso.CV <- cv.glmnet(x, y, family = "binomial", 
                      glmnet(x, y, family = "gaussian", alpha = 0)
                      type.measure = "auc", #"auc","mse"
                      standardize = T,
                      nfolds = 10) 
#plot(lasso.CV) 
plot(lasso.CV, xvar = "lambda", label = TRUE)

#1
lasso.CV$lambda.1se 
LassoVIP1 <- as.matrix(coef(lasso.CV, s = "lambda.1se")) 
LassoVIP2 <- data.frame(LassoVIP1)
LassoVIP2$ABS <- abs(LassoVIP2[,1]) 
LassoVIP2 <- LassoVIP2[-1,] 
LassoVIP2 <- LassoVIP2[LassoVIP2$ABS > 0,] 
LassoVIP2 <- LassoVIP2[order(LassoVIP2$ABS, decreasing = TRUE), ] 
write.csv(LassoVIP2, paste0(outputDir, "/", "All statistic_Features_LASSOlse.csv"))
#2
lasso.CV$lambda.min #lambda.min 
LassoVIP3 <- as.matrix(coef(lasso.CV, s = "lambda.min")) 
LassoVIP4 <- data.frame(LassoVIP2)
LassoVIP4$ABS <- abs(LassoVIP4[,1]) 
LassoVIP4 <- LassoVIP4[-1,] 
LassoVIP4 <- LassoVIP4[LassoVIP4$ABS > 0,] 
LassoVIP4 <- LassoVIP4[order(LassoVIP4$ABS, decreasing = TRUE), ] 
write.csv(LassoVIP4, paste0(outputDir, "/", "All statistic_Features_LASSOmin.csv"))

#################################################################
#--------AUC--------
Model_Markers <- read.csv("All data statistic_AUC markers.csv", header = TRUE) #
Trans2 = F #T or F
Train.data = F  #T or F
t <- 0.7 

Metaname <- Model_Markers
cols_remain <- factor(Metaname[,1]) 
Model_data2 <- data[, colnames(data) %in% cols_remain]

if(Trans2 == T){ 
  int <- function(x){
    qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x)))
  }
  Model_data2 <- apply(Model_data2, 2, int)#inverse normal transformation
  Model_data2 <- scale(Model_data2) #Scale
}else{ 
  Model_data2 <- Model_data2 
}

Model_data <- data.frame(data[1], Model_data2)
Model_data$Group <- as.numeric(ifelse(Model_data$Group == "Control",'0','1')) 

if(Train.data == T){ 
  set.seed(1234)
  train <- sample(nrow(Model_data), t*nrow(Model_data))
  train.data <- Model_data[train,]
  test.data <- Model_data[-train,]
}else{ ）
  train.data <- Model_data 
  test.data <- Model_data 
}

#------logistic model-------
full.fit <- glm(Group ~ ., family = binomial(), data = train.data)
summary(full.fit)
vif(full.fit)

logt1 <- coef(summary(full.fit))
OR <- exp(coef(full.fit))
CI <- exp(confint(full.fit))
logt3 <- cbind(logt1, OR, CI)
names(logt3) <- c("Estimate", "Std. Error", "z value", "pvalue", "OR", 
                  "OR.lower", "OR.upper")
write.csv(logt3, paste0(outputDir, "/", "All statistic_predict_Logis_OR_train data.csv"))

#-------------train.data--------------------
model_pre1 <- predict(full.fit, newdata = train.data, type = "response")
predout1 <- cbind(train.data$Group, model_pre1)
#write.csv(predout1, paste0(outputDir, "/", "All statistic_prediction_Logis_pre.csv"))

report1 <- as.data.frame(reportROC(gold = train.data$Group, predictor = model_pre1))
write.csv(report1, paste0(outputDir, "/", "All statistic_predict_Logis_AUC_train data.csv"))

#-------------test.data--------------------
model_pre2 <- predict(full.fit, newdata = test.data, type = "response")
predout2 <- cbind(test.data$Group, model_pre2)
write.csv(predout2, paste0(outputDir, "/", "All statistic_predict_Logis_pre_test data.csv"))

report2 <- as.data.frame(reportROC(gold = test.data$Group, predictor = model_pre2))
write.csv(report2, paste0(outputDir, "/", "All statistic_predict_Logis_AUC_test data.csv"))

#-------------ROC-----------------------
roc1 <- roc((test.data$Group), model_pre2, ci = T)
round(auc(roc1),3)
ci(roc1)

plot(roc1, print.auc = T,  auc.polygon = T, 
     auc.polygon.col = "white", print.thres = T)

#############################################################
#------Pathway scores Calculation----------
rm(list=ls())
data <- read.csv("Calculate Score data.csv", header = TRUE, row.name = 1) 

#LASSO-New Pathway-All features, mg/dL
Pathway1 <- -2.017 +0.138*data$AST
Pathway2 <- 4.941 -0.420*data$BMI +0.116*data$H4LP #mg/dL
Pathway3 <- -0.442 +0.067*data$V1TG #mg/dL
Pathway4 <- 2.607 +7.530*data$Gly -2.836*data$Ala -20.624*data$AcOH -15.421*data$Val
Pathway5 <- 6.535 +73.370*data$Acaa -96.054*data$Tyr -25.585*data$Ile -24.670*data$Acto

#mmol/L for UKB and Finnish
Pathway2_F <- 4.941 -0.420*data$BMI +0.116*55.1*data$H4LP #mmol/L
Pathway3_F <- -0.442 +0.067*88.5*data$V1TG #mmol/L

name <- rownames(data)
data_out <- data.frame('Metabolites' = name[1:length(name)], 
                       Pathway1, Pathway2, Pathway3,  
                       Pathway4, Pathway5)
write.csv(data_out,"Calculate Score_results.csv", row.names = F, quote = F)

###################################################################
#----K-Means cluster--------------------------
library(cluster)
library(factoextra) 
library(dplyr)
library(pacman)
library(NbClust) 
library(fMultivar)
library(ggplot2)

rm(list=ls())
data1 <- read.csv("Cluster data.csv", header = TRUE, row.name = 1) 

Metaname <- read.csv("Cluster_markers.csv", header = TRUE) #c("Pathway1"..."Pathway5")
cols_remain <- factor(Metaname[,1]) 
data <- data1[, colnames(data1) %in% cols_remain]

data <- na.omit(data) 
data.scale <- data.frame(scale(data)) 

d <- dist(data.scale, method = "euclidean") 

#-----number of clusters-----
set.seed(123)
fviz_nbclust(data.scale, kmeans, k.max = 15)

set.seed(666) 
kmeans1 <- kmeans(data.scale, centers = 3, nstart = 25) 

#PCA polt
fviz_cluster(object = kmeans1, data = data.scale,
             ellipse.type = "euclid", star.plot = T,
             repel = T, geom = ("point"), palette = c("#FF6B6B", "#4ECDC4", "#556270"),
             main = "Cluster plot with k = 3 (K-means)", ggtheme = theme_minimal()) +
  theme(panel.background = element_rect(fill ="white", colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "white", colour = "black"))

out <- data.frame(data, kmeans1$cluster)
write.csv(out,"Kmeans_3Clusters.csv")

# average pathway scores for each cluster
names(out) <- c("P1", "P2", "P3", "P4", "P5", "Cluster")
P1 <- aggregate(P1 ~ Cluster, data = out, FUN = mean)
P2 <- aggregate(P2 ~ Cluster, data = out, FUN = mean)
P3 <- aggregate(P3 ~ Cluster, data = out, FUN = mean)
P4 <- aggregate(P4 ~ Cluster, data = out, FUN = mean)
P5 <- aggregate(P5 ~ Cluster, data = out, FUN = mean)

data_out <- data.frame(P1, P2[,2], P3[,2], P4[,2], P5[,2])
write.csv(data_out,"Kmeans_3Clusters_NMR_CN Biops-V4_average score.csv", row.names = F, quote = F)

#--------Jaccard------
library(fpc)
cb <- clusterboot(data.scale, B = 1000, 
                  clustermethod = kmeansCBI, 
                  k = 3, seed = 123, 
                  bootmethod = "boot") 
cb$bootmean 
mean(cb$bootmean) 
sd(cb$bootmean)

################################################################
#1: Multi-group spider chart, no options
library(fmsb)
library(RColorBrewer)
library(scales)

rm(list = ls())
data <- read.csv("SpiderChart data.csv", header = T, row.names = 1)

coul <- c("#FF6B6B","#4ECDC4","#556270","#2171B5","green","#E7B800") # Set graphic colors
colors_border <- coul 
colors_in <- alpha(coul, 0.2) 

radarchart(data, axistype = 0 , maxmin = F,
           pcol = colors_border, pfcol = colors_in, 
           plty = 1, plwd = 1.5, 
           cglcol = "grey", cglty = 1, axislabcol = "black", cglwd = 1, 
           caxislabels = seq(0, 20, 5), 
           vlcex = 1) 
# Add a legend
legend(x = 0.85, y = 1.3, legend = rownames(data), 
       bty = "n", pch = 20 , col = colors_in , text.col = "black", 
       cex = 0.7, pt.cex = 2)
title(main = "Treat1-3Cluster", #sub = "Geeksforgeeks article",
      #xlab = "X axis", ylab = "Y axis",
      cex.main = 1,   font.main = 1, col.main = "black",
      #cex.sub = 0.8, font.sub = 1, col.sub = "darkgreen",
      col.lab ="black")

#2: Multi-group spider chart, with options
library(fmsb) 
library(RColorBrewer)
library(scales)
library(dplyr)

rm(list = ls())
data <- read.csv("SpiderChart data.csv", header = T)

mat <- data %>% 
  select(-Group) %>% 
  as.matrix() %>% 
  rbind(rep(2, 5), rep(-0.8, 5), .)   # Max / Min
rownames(mat) <- c("Max", "Min", data$Group)
mat <- data.frame(mat)

# Set graphic colors
coul <- c("#FF6B6B", "#4ECDC4", "#556270","#2171B5", "#E7B800", "green")
colors_border <- coul 
colors_in <- alpha(coul, 0.2) 

radarchart(mat, axistype = 1, 
           pcol = colors_border, pfcol = colors_in, 
           plty = 1, plwd = 1.5, 
           cglcol = "grey", cglty = 1, axislabcol = "black", cglwd = 1, 
           vlcex = 1)
# Add a legend
legend(x = 0.7, y = 0.95, legend = rownames(mat[-c(1,2),]), 
       bty = "n", pch = 20, col = colors_in, text.col = "black", 
       cex = 0.8, pt.cex = 2)
title(main = "Changfeng and UKB-3Cluster", cex.main = 1,
      font.main = 1, col.main = "black", col.lab ="black")

##################################################################
#---------circlize heatmap
library(circlize)
library(ComplexHeatmap)
library(magrittr)
library(dendsort)
library(dplyr)

rm(list=ls())
data <- read.csv("Circos data_Group.csv", header = T, row.names = 1)

mat1 <- data %>% select(-Group)
mat2 <- t(scale(t(mat1)))
data$Group <- factor(data$Group,levels = data$Group %>% unique())
split <- data$Group 

circos.clear() 
range(mat2) 
color <- colorRamp2(c(-1.5, 0, 1.5), c("#4ECDC4", "white", "#FF6B6B"))

circos.par(start.degree = 90, 
           gap.after = c(5,5,5,30)) 

circos.heatmap(mat2, split = split,
               col = color, 
               dend.side = "inside", 
               dend.track.height = 0.1, 
               rownames.side = "outside", 
               rownames.col = "black", 
               rownames.cex = 0.6, 
               rownames.font = 1, 
               track.height = 0.3, 
               #clustering.method = "complete", 
               #distance.method = "euclidean",
               cluster = T) 

circos.trackPlotRegion(ylim = c(0, 0.1), track.height = 0.05, 
                       bg.col = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1" ,"#6BAED6", "#4292C6" ,"#2171B5" ,"#084594"), 
                       panel.fun = function(x, y) {
                         chr = CELL_META$sector.index
                         xlim = CELL_META$xlim
                         ylim = CELL_META$ylim
                         circos.text(mean(xlim), mean(ylim),
                                     chr, cex = 0.4, col = "black",
                                     facing = "bending.inside", 
                                     niceFacing = F)
                       })

lgd = Legend(title = "Z-score", col_fun = color, direction = "horizontal",
             legend_width = unit(3,"cm"))
draw(lgd, x = unit(0.3,"npc"), y = unit(0.9,"npc"))

############################################################
#------pie plot
library(ggplot2)

rm(list=ls())
data <- read.csv("pie plot data.csv", header = T)

myLabel = as.vector(data$category) 
data$fraction = data$count/sum(data$count)
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n = -1))

# Make the plot
ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, 
                 fill = category)) +
  geom_rect(alpha = 0.8) + 
  coord_polar(theta = "y") + 
  scale_fill_manual(values = c("#FF6B6B", "#4ECDC4", "#556270", "#E7B800")) +
  xlim(c(2, 4)) + 
  labs(x = "", y = "", title = "Treat1_2 cohorts") + 
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.border = element_blank(), 
    legend.title = element_blank(),
    legend.position = c(0.5, 0.5)
  )  +
  geom_text(aes(x = 3.5, y = ((ymin + ymax)/2), label = count)) 

######################################################################
###---------------Cox analysis for multi-group----------------------------
library(survival) 
library(survminer) 
library(devtools) 
library(jstable) 
library(cli) 
library(forestploter) 
library(grid) 

rm(list=ls())
data <- read.csv("Cox HR data_Multigroup.csv", header = TRUE, row.name = 1) 

cox_data <- data
data$Cluster <- as.factor(data$Cluster) 

#--------non-adjust
res.cox1 <- coxph(Surv(time = cox_data$time, event = cox_data$Group == 1) ~ Cluster, 
                  data = cox_data)  
summary(res.cox1)
#pairwise log-rank
library(multcomp)
contr1 <- glht(res.cox1, linfct = mcp(Cluster = "Tukey"))
summary(contr1)

#----adjustment for age, sex and smoke
res.cox2 <- coxph(Surv(time = cox_data$time, event = cox_data$Group == 1) ~ Cluster+Sex+Age+Smoke, 
                  data = cox_data)  
summary(res.cox2)
#pairwise log-rank
contr2 <- glht(res.cox2, linfct = mcp(Cluster = "Tukey"))
summary(contr2)

#------------survival plot---------
library(survival) 
library(survminer) 

rm(list=ls())
data1 <- read.csv("Coxplot data.csv", header = TRUE)
coxdata <- data1

fit <- survfit(Surv(time = coxdata$time, event = coxdata$Group == 1) ~ Cluster, 
               data = coxdata)

A <- ggsurvplot(fit, data = coxdata, 
                ggtheme = theme_bw(),
                conf.int = F, 
                pval = T, 
                pval.coord = c(0, 0.02), 
                legend.labs = c("Cluster1", "Cluster2", "Cluster3"), 
                legend.title = "", 
                ylim = c(0, 0.03), 
                #surv.median.line = "hv", 
                censor = FALSE, 
                fun = "cumhaz"  
) +
  ggtitle("Changfeng and UK Biobank") +
  xlab("Follow-up time (years)") + ylab("Cumulative mortality")
A$plot <- A$plot + 
  scale_color_manual(values = c("#FF6B6B", "#4ECDC4", "#556270", "#E7B800")) +
  theme(
    legend.background = element_rect(fill="white", colour = "white"),
    legend.position = "right", 
    panel.grid = element_blank(), 
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 14, colour = "black", vjust = 0.3, hjust = 0.3),
    axis.title = element_text(size = 14, colour = "black"),
  ) 
A

#############################################################
#----------ANOVA analysis, (including pairwise comparisons)

rm(list=ls())
data <- read.csv("ANOVA data.csv", header = TRUE, row.name = 1)
Anov_data <- data

categ <- 0 
e <- 0 

a <- 1 + categ 
b <- 1 + categ + e 
c <- a + 1 
d <- b + 1 

#inverse normal transformation
int <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  #qnorm(rank(x)/length(x))
}
Anov_data[, d:ncol(Anov_data)] <- apply(Anov_data[, d:ncol(Anov_data)], 2, int)
write.csv(Anov_data,'Anov_data_transformed.csv', row.names = T, quote = F)

#------non adjustment-------
Ano1 <- aov(unlist(Anov_data$FFAll_Delta) ~ Group, data = Anov_data)
summary(Ano1)
#pairwise comparisons
TukeyHSD(Ano1)

#----------adjustment for age, sex and baseline FFA----------
Ano2 <- aov(unlist(Anov_data$FFAll_Delta) ~ Group+Age+Sex+FFAll_Base, 
            data = Anov_data)
summary(Ano2)
#pairwise comparisons
library(emmeans) # Adjusted mean change
emmeans_result <- emmeans(Ano2, pairwise ~ Group, adjust = "none")
summary(emmeans_result, infer = TRUE)

################################################################
##-------------Proteome enrichment
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(DOSE)
library(ggnewscale)
library(circlize)
library(ComplexHeatmap)

rm(list=ls())
info <- read.csv("Proteome_enrich data.csv", header = TRUE)

#GO_database
GO_database <- 'org.Hs.eg.db'  # http://bioconductor.org/packages/release/BiocViews.html#___OrgDb

gene <- bitr(info$gene_symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

GO <- enrichGO(gene$ENTREZID, 
               OrgDb = GO_database,
               keyType = "ENTREZID", 
               ont = "ALL", #(ALL: Biological Process,Cellular Component,Mollecular Function）
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05, 
               readable = T)
write.csv(GO,'Proteome_enrichment_GO.csv', row.names = T, quote = F)

barplot(GO, showCategory = 10, split = "ONTOLOGY") + facet_grid(ONTOLOGY~., scale = "free")

enrichplot::cnetplot(GO, showCategory = 5, color_category = "blue", size_category = 1)


################################################################
##-------------Multi-group: bar plot
library(ggpubr)
library(reshape2)
library(ggplot2)
library(dplyr)

rm(list=ls())

Mean = read.csv("barplot_Mean.csv", header = T)
data1 <- melt(Mean,id = "Metabolites") 
names(data1) <- c("Metabolites", "Group", "Mean")

SD = read.csv("barplot_SD.csv", header = T)
data2 <- melt(SD,id = "Metabolites") 
names(data2) <- c("Metabolites", "Group", "SD")

data <- data.frame(data1, data2$SD)
names(data) <- c("Metabolites", "Group", "Mean", "SD")

#bar plot
ggplot(data, aes(x = Metabolites, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(0.8)) + 
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, position = position_dodge(0.8)) +
  labs(title = "Lipoproteins",
       x = "Groups",
       y = "Relative Intensity",
       fill = "Cluster") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_fill_manual(values = c("#FF6B6B", "#4ECDC4", "#556270", "#E7B800")) +
  #scale_color_manual(values = c("#FF6B6B", "#4ECDC4", "#556270", "#E7B800")) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.line = element_line(linewidth = 0.5, colour = 'black')
  )


