#小周的bulkRNA数据

setwd("/path/to/your/own/Input")

# Import packages ---------------------------------------------------------
library(edgeR)
library(readr)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)

#数据导入
PL4 <- read.table("PL4_featureCounts_output.txt", sep = "", header = T)
PL5 <- read.table("PL5_featureCounts_output.txt", sep = "", header = T)
PL11 <- read.table("PL11_featureCounts_output.txt", sep = "", header = T)
PL12 <- read.table("PL12_featureCounts_output.txt", sep = "", header = T)
V4 <- read.table("V4_featureCounts_output.txt", sep = "", header = T)
V5 <- read.table("V5_featureCounts_output.txt", sep = "", header = T)
V10 <- read.table("V10_featureCounts_output.txt", sep = "", header = T)
V12 <- read.table("V12_featureCounts_output.txt", sep = "", header = T)

#列名を変更する
colnames(PL4) <- c("Gene_ID", "PL4")
colnames(PL5) <- c("Gene_ID", "PL5")
colnames(PL11) <- c("Gene_ID", "PL11")
colnames(PL12) <- c("Gene_ID", "PL12")
colnames(V4) <- c("Gene_ID", "V4")
colnames(V5) <- c("Gene_ID", "V5")
colnames(V10) <- c("Gene_ID", "V10")
colnames(V12) <- c("Gene_ID", "V12")

#这b还轮流merge也是真爱护环境
Data <- merge(PL4, PL5, all = T, by = "Gene_ID")
Data <- merge(Data, PL11, all = T, by = "Gene_ID")
Data <- merge(Data, PL12, all = T, by = "Gene_ID")
Data <- merge(Data, V4, all = T, by = "Gene_ID")
Data <- merge(Data, V5, all = T, by = "Gene_ID")
Data <- merge(Data, V10, all = T, by = "Gene_ID")
Data <- merge(Data, V12, all = T, by = "Gene_ID")

#这一步是把第一列提出来作为行名，然后把行名赋值给Data这个数据后把第一列删掉
row_names <- Data[,1]
rownames(Data) <- row_names
Data <- Data[,-1]

#她想把这个表现矩阵存下来
#write.csv(Data,  "./231110_Data integration/Output/231110 read count.csv") ##save csv file

#cut-offf(read count >10) ---------------------------------过滤了一下
data <- filter_all(Data, all_vars(. >=10))

#她又想把过滤完的矩阵存下来
#write.csv(data,  "./231110_Data integration/Output/231110 read count>10.csv") ##save csv file

# Calculate logCPM ---好像是把原来的raw counts给标准化为一个计数形式-------------------------------------------
logcpm <- cpm(Data, log=TRUE)
logcpm

#write.csv(Data,  "./231110_Data integration/Output/231110 logCPM.csv")

# Hierarchical clustering --计算相关性然后聚类------------------------------------
#Hierarchical clustering between samples
#use Spearman's rank correlation coefficient

rho <- cor(logcpm, method = "spearman")

#Ward-based clustering
d <- as.dist(1 - rho)
tree <- hclust(d, method = "ward.D")
plot(tree) #看聚类的树

# Sample PCA ----------主成分分析------------------------------------------
pca <- prcomp(t(logcpm), scale = FALSE)
summary(pca)

#Importance of components:
#                            PC1     PC2     PC3     PC4     PC5      PC6      PC7
# Standard deviation     37.6624 29.1613 28.6927 25.9188 23.4865 22.41439 22.15855
# Proportion of Variance  0.2672  0.1602  0.1551  0.1265  0.1039  0.09463  0.09249
# Cumulative Proportion   0.2672  0.4274  0.5824  0.7090  0.8129  0.90751  1.00000
#                              PC8
# Standard deviation     4.401e-13
# Proportion of Variance 0.000e+00
# Cumulative Proportion  1.000e+00

#plot ---可视化PCA-------------------------------------------------
color <- c("#97ffb8","#ff7b7b", "#FFC354","#85fbff", "#509DC3","#be8fff", "#fe75a2","#dcff7b")
par(mar=c(5,5,5,10))
# 有边框版本
plot(pca$x[, 1], pca$x[,2], 
     col = "black",  # 先设定边框颜色
     bg = color,     # 再填充里面的颜色
     pch = 21,       # 要边框是21，不要就是19
     cex = 2,        # 点的大小
     main = "PCA", 
     xlab = "PC1 (37.7%)", 
     ylab = "PC2 (29.2%)")
par(xpd=TRUE)
legend(x=par()$usr[2], y=par()$usr[4], legend = colnames(logcpm), bty = "n", pch = 19,  col  = color, pt.cex =1.5)

# 无边框版本
plot(pca$x[, 1], pca$x[,2], 
     col = color,    #颜色就是设定好的颜色
     pch = 19,       # 要边框是21，不要就是19
     cex = 2,        # 点的大小
     main = "PCA", 
     xlab = "PC1 (37.7%)", 
     ylab = "PC2 (29.2%)")
par(xpd=TRUE)
legend(x=par()$usr[2], y=par()$usr[4], legend = colnames(logcpm), bty = "n", pch = 19,  col  = color, pt.cex =1.5)

#DEG
#你那个代码她把保存了的过滤后矩阵又读进来弄的，我是直接环境里还有，就没重新读，继续

# 把data变成矩阵
data <- as.matrix(data)
dim(data)

# 定义一下组别，就是在这里把两个V合成一个Vehicle的
Treatment <- factor(c("aPD-L1", "aPD-L1", "KO_aPD-L1", "KO_aPD-L1", "Vehicle", "Vehicle", "KO_Vehicle", "KO_Vehicle"))
data <- DGEList(data, group = Treatment)

# 过滤一下低表达基因
keep <- filterByExpr(data)
data <- data[keep, , keep.lib.sizes=FALSE]

# 标准化
data <- calcNormFactors(data)
data$samples
design <- model.matrix(~Treatment)
colnames(design) <- levels(Treatment)   # 将列名设置为分组名
data <- estimateCommonDisp(data)

# 提出差异表达基因DEG-------------------------------------------------------------
Treatment
# [1] aPD-L1     aPD-L1     KO_aPD-L1  KO_aPD-L1  Vehicle    Vehicle   
# [7] KO_Vehicle KO_Vehicle
# Levels: aPD-L1 KO_aPD-L1 KO_Vehicle Vehicle

#切FC是1.5的话那logFC切到0.585

### aPD-L1 vs Vehicle ###
et <- exactTest(data, pair = c("aPD-L1", "Vehicle"))
topTags(et)
FC1 <- as.data.frame(topTags(et, n = nrow(count)))
FC1 <- FC1 %>% 
  mutate(Expression = case_when(logFC >= 0.585 & logCPM >= 0.301 & FDR < 0.05 ~ "Up-regulated",           
                                logFC <= -0.585 & logCPM >= 0.301 & FDR < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))##FC>1.5, logCPM >1, FDR <0.05
# write.csv(FC1, "/Users/asakawayuuki/Desktop/がん研/RNAseq/小周的结果保存/DEG_aPD-L1 vs Vehicle.csv")



### KO_aPD-L1 vs KO_Vehicle ###
et <- exactTest(data, pair = c("KO_aPD-L1", "KO_Vehicle"))
topTags(et)
FC2 <- as.data.frame(topTags(et, n = nrow(count)))
FC2 <- FC2 %>% 
  mutate(Expression = case_when(logFC >= 0.585 & logCPM >= 0.301 & FDR < 0.05 ~ "Up-regulated",           
                                logFC <= -0.585 & logCPM >= 0.301 & FDR < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))##FC>1.5, logCPM >1, FDR <0.05
# write.csv(FC2, "/Users/asakawayuuki/Desktop/がん研/RNAseq/小周的结果保存/DEG_KO_aPD-L1 vs KO_Vehicle.csv")


### KO_Vehicle vs Vehicle ###
et <- exactTest(data, pair = c("KO_Vehicle", "Vehicle"))
topTags(et)
FC3 <- as.data.frame(topTags(et, n = nrow(count)))
FC3 <- FC3 %>% 
  mutate(Expression = case_when(logFC >= 0.585 & logCPM >= 0.301 & FDR < 0.05 ~ "Up-regulated",           
                                logFC <= -0.585 & logCPM >= 0.301 & FDR < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))##FC>1.5, logCPM >1, FDR <0.05
# write.csv(FC3, "/Users/asakawayuuki/Desktop/がん研/RNAseq/小周的结果保存/DEG_KO_Vehicle vs Vehicle.csv")


### KO_aPD-L1 vs aPD-L1 ###
et <- exactTest(data, pair = c("KO_aPD-L1", "aPD-L1"))
topTags(et)
FC4 <- as.data.frame(topTags(et, n = nrow(count)))
FC4 <- FC4 %>% 
  mutate(Expression = case_when(logFC >= 0.585 & logCPM >= 0.301 & FDR < 0.05 ~ "Up-regulated",           
                                logFC <= -0.585 & logCPM >= 0.301 & FDR < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))##FC>1.5, logCPM >1, FDR <0.05
# write.csv(FC4, "/Users/asakawayuuki/Desktop/がん研/RNAseq/小周的结果保存/DEG_KO_aPD-L1 vs PD-L1.csv")


#火山图
### aPD-L1 vs Vehicle ###
options(ggrepel.max.overlaps = Inf)
MA1 <- ggplot(FC1, aes(x=logFC, y=logCPM, colour=Expression)) + 
  geom_point(alpha=0.8) + 
  ylim(-3,20) +
  theme_bw() + 
  theme(
    legend.position = "right",            
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("#2166ac", "grey","#b2182b")) +
  labs(colour="Differential Expression") + 
  ggtitle(paste0("MA Plot (aPD-L1 vs Vehicle)")) +
  geom_vline(xintercept=c(-0.585,0.585),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept=0.301,lty=2,col="black",lwd=0.4) +  
  labs(x="log2(FC)",y="log2CPM")+guides(colour = guide_legend(override.aes = list(size = 4)))

MA1


### KO_aPD-L1 vs KO_Vehicle ###
options(ggrepel.max.overlaps = Inf)
MA2 <- ggplot(FC2, aes(x=logFC, y=logCPM, colour=Expression)) + 
  geom_point(alpha=0.8) + 
  ylim(-3,20) +
  theme_bw() + 
  theme(
    legend.position = "right",            
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("#2166ac", "grey","#b2182b")) +
  labs(colour="Differential Expression") + 
  ggtitle(paste0("MA Plot (KO_aPD-L1 vs KO_Vehicle)")) +
  geom_vline(xintercept=c(-0.585,0.585),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept=0.301,lty=2,col="black",lwd=0.4) +  
  labs(x="log2(FC)",y="log2CPM")+guides(colour = guide_legend(override.aes = list(size = 4)))

MA2

### KO_Vehicle vs Vehicle ###
options(ggrepel.max.overlaps = Inf)
MA3 <- ggplot(FC3, aes(x=logFC, y=logCPM, colour=Expression)) + 
  geom_point(alpha=0.8) + 
  ylim(-3,20) +
  theme_bw() + 
  theme(
    legend.position = "right",            
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("#2166ac", "grey","#b2182b")) +
  labs(colour="Differential Expression") + 
  ggtitle(paste0("MA Plot (KO_Vehicle vs Vehicle)")) +
  geom_vline(xintercept=c(-0.585,0.585),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept=0.301,lty=2,col="black",lwd=0.4) +  
  labs(x="log2(FC)",y="log2CPM")+guides(colour = guide_legend(override.aes = list(size = 4)))

MA3

### KO_aPD-L1 vs aPD-L1 ###
options(ggrepel.max.overlaps = Inf)
MA4 <- ggplot(FC4, aes(x=logFC, y=logCPM, colour=Expression)) + 
  geom_point(alpha=0.8) + 
  ylim(-3,20) +
  theme_bw() + 
  theme(
    legend.position = "right",            
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("#2166ac", "grey","#b2182b")) +
  labs(colour="Differential Expression") + 
  ggtitle(paste0("MA Plot (KO_aPD-L1 vs aPD-L1)")) +
  geom_vline(xintercept=c(-0.585,0.585),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept=0.301,lty=2,col="black",lwd=0.4) +  
  labs(x="log2(FC)",y="log2CPM")+guides(colour = guide_legend(override.aes = list(size = 4)))

MA4

## MA plot with label(top10 gene)带标识的-------------------------------

### aPD-L1 vs Vehicle ###
top <- 10
top_genes_1 <- bind_rows(
  FC1 %>% 
    filter(Expression == 'Up-regulated' ) %>% 
    arrange(desc(logFC)) %>%  head(top), FC1 %>% 
    filter(Expression == 'Down-regulated' ) %>% 
    arrange(logFC) %>%  head(top))

top_genes_1

# 得把这个行名转换为一列单独存在，取名 "Gene"
top_genes_1 <- top_genes_1 %>%
  rownames_to_column(var = "Gene")

# write.csv(top_genes_1,"./231110_Data integration/Output/Top genes/Top genes_RK-582 vs Vehicle.csv" )
# top_genes_1 <- read.csv("./231110_Data integration/Output/Top genes/Top genes_RK-582 vs Vehicle.csv" )
# colnames(top_genes_1) <- c("Gene", "logFC" , "logCPM", "PValue", "Expression")

ma1 <- ggplot(FC1, aes(x=logFC, y=logCPM, colour=Expression)) + 
  geom_point(alpha=0.8) + 
  ylim(-3,20) +
  theme_bw() + 
  theme(
    legend.position = "right",            
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("#2166ac", "grey","#b2182b")) +
  labs(colour="Differential Expression") + 
  ggtitle(paste0("MA Plot (aPD-L1 vs Vehicle)")) +
  geom_vline(xintercept=c(-0.585,0.585),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept=0.301,lty=2,col="black",lwd=0.4) +  
  labs(x="log2(FC)",y="log2CPM")+guides(colour = guide_legend(override.aes = list(size = 4)))+ geom_text_repel(data = top_genes_1, label = (top_genes_1$Gene), box.padding = unit(.100, "lines"),hjust= 0.30, size = 3, force = 1, nudge_y = 1, color ="black") 


ma1


### KO_aPD-L1 vs KO_Vehicle ###
top_genes_2 <- bind_rows(
  FC2 %>% 
    filter(Expression == 'Up-regulated' ) %>% 
    arrange(desc(logFC)) %>%  head(top), FC2 %>% 
    filter(Expression == 'Down-regulated' ) %>% 
    arrange(logFC) %>%  head(top))

top_genes_2

top_genes_2 <- top_genes_2 %>%
  rownames_to_column(var = "Gene")

ma2 <- ggplot(FC2, aes(x=logFC, y=logCPM, colour=Expression)) + 
  geom_point(alpha=0.8) + 
  ylim(-3,20) +
  theme_bw() + 
  theme(
    legend.position = "right",            
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("#2166ac", "grey","#b2182b")) +
  labs(colour="Differential Expression") + 
  ggtitle(paste0("MA Plot (KO_aPD-L1 vs KO_Vehicle)")) +
  geom_vline(xintercept=c(-0.585,0.585),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept=0.301,lty=2,col="black",lwd=0.4) +  
  labs(x="log2(FC)",y="log2CPM")+guides(colour = guide_legend(override.aes = list(size = 4))) +
  geom_text_repel(data = top_genes_2, label = (top_genes_2$Gene), box.padding = unit(.100, "lines"),hjust= 0.30, size = 3, force = 1, nudge_y = 1, color ="black") 

ma2


### KO_Vehicle vs Vehicle ###
top_genes_3 <- bind_rows(
  FC3 %>% 
    filter(Expression == 'Up-regulated' ) %>% 
    arrange(desc(logFC)) %>%  head(top), FC3 %>% 
    filter(Expression == 'Down-regulated' ) %>% 
    arrange(logFC) %>%  head(top))

top_genes_3

top_genes_3 <- top_genes_3 %>%
  rownames_to_column(var = "Gene")

ma3 <- ggplot(FC3, aes(x=logFC, y=logCPM, colour=Expression)) + 
  geom_point(alpha=0.8) + 
  ylim(-3,20) +
  theme_bw() + 
  theme(
    legend.position = "right",            
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("#2166ac", "grey","#b2182b")) +
  labs(colour="Differential Expression") + 
  ggtitle(paste0("MA Plot (KO_Vehicle vs Vehicle)")) +
  geom_vline(xintercept=c(-0.585,0.585),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept=0.301,lty=2,col="black",lwd=0.4) +  
  labs(x="log2(FC)",y="log2CPM")+guides(colour = guide_legend(override.aes = list(size = 4))) +
  geom_text_repel(data = top_genes_3, label = (top_genes_3$Gene), box.padding = unit(.100, "lines"),hjust= 0.30, size = 3, force = 1, nudge_y = 1, color ="black") 

ma3


### KO_aPD-L1 vs aPD-L1 ###
top_genes_4 <- bind_rows(
  FC4 %>% 
    filter(Expression == 'Up-regulated' ) %>% 
    arrange(desc(logFC)) %>%  head(top), FC4 %>% 
    filter(Expression == 'Down-regulated' ) %>% 
    arrange(logFC) %>%  head(top))

top_genes_4

top_genes_4 <- top_genes_4 %>%
  rownames_to_column(var = "Gene")

ma4 <- ggplot(FC4, aes(x=logFC, y=logCPM, colour=Expression)) + 
  geom_point(alpha=0.8) + 
  ylim(-3,20) +
  theme_bw() + 
  theme(
    legend.position = "right",            
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("#2166ac", "grey","#b2182b")) +
  labs(colour="Differential Expression") + 
  ggtitle(paste0("MA Plot (KO_aPD-L1 vs aPD-L1)")) +
  geom_vline(xintercept=c(-0.585,0.585),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept=0.301,lty=2,col="black",lwd=0.4) +  
  labs(x="log2(FC)",y="log2CPM")+guides(colour = guide_legend(override.aes = list(size = 4))) +
  geom_text_repel(data = top_genes_4, label = (top_genes_4$Gene), box.padding = unit(.100, "lines"),hjust= 0.30, size = 3, force = 1, nudge_y = 1, color ="black") 

ma4

# # Session.information ---看当前环境的没啥b用不用管---------------------------------------------
# sessionInfo()