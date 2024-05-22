# This script will construct a Volcano plot that describes the genes found differentially partitioned between mesophyll and bundle sheath cell types
# at either the 0h time point or 12 h time point. 

setwd("~/Desktop/Salk/Project - Light/R files")
load("L3_sorghum_T0_response_10X.RData")
load("L3_sorghum_T12_response_10X.RData")
T0_response_specific<-subset(T0_response, T0_response$p_val_adj<0.05)
T12_response_specific<-subset(T12_response, T12_response$p_val_adj<0.05)

# Volcano Plot
sorghum_lightdark_extensive_C4<-subset(sorghum_lightdark_extensive, sorghum_lightdark_extensive$category == 'C4')
sorghum_lightdark_extensive_calvin<-subset(sorghum_lightdark_extensive, sorghum_lightdark_extensive$category %in% c('calvin_cycle','photorespiration'))
sorghum_lightdark_extensive_light<-subset(sorghum_lightdark_extensive, sorghum_lightdark_extensive$category %in% c('light_reaction'))

# T0
T0_response_specific[,6]<-ifelse(rownames(T0_response_specific) %in% sorghum_lightdark_extensive_C4$gene_ID, 2, 0)
T0_response_specific[,7]<-ifelse(rownames(T0_response_specific) %in% sorghum_lightdark_extensive_calvin$gene_ID, 3, 0)
T0_response_specific[,8]<-ifelse(rownames(T0_response_specific) %in% sorghum_lightdark_extensive_light$gene_ID, 4, 0)
T0_response_specific[,9]<-T0_response_specific[,6]+T0_response_specific[,7]+T0_response_specific[,8]
T0_response_specific$p_val_adj<--log(T0_response_specific$p_val_adj,10)
T0_response_specific[T0_response_specific == Inf] <- 300

T0_plot_sorghum<-ggplot(T0_response_specific, aes(x=avg_log2FC, y=p_val_adj, color = as.factor(T0_response_specific[,9]), size = as.factor(T0_response_specific[,9]))) + geom_point() + theme_bw() + scale_color_manual(values=c('grey80','red','purple3','yellow3')) + xlim(c(-3,3)) + geom_vline(xintercept = 0, color = "black", size=.5) + xlab("log2 FC") + ylab("-log10 p-value") + theme(legend.position = 'none') 
# T12
T12_response_specific[,6]<-ifelse(rownames(T12_response_specific) %in% sorghum_lightdark_extensive_C4$gene_ID, 2, 0)
T12_response_specific[,7]<-ifelse(rownames(T12_response_specific) %in% sorghum_lightdark_extensive_calvin$gene_ID, 3, 0)
T12_response_specific[,8]<-ifelse(rownames(T12_response_specific) %in% sorghum_lightdark_extensive_light$gene_ID, 4, 0)
T12_response_specific[,9]<-T12_response_specific[,6]+T12_response_specific[,7]+T12_response_specific[,8]

T12_response_specific$p_val_adj<--log(T12_response_specific$p_val_adj,10)
T12_response_specific[T12_response_specific == Inf] <- 300
T12_plot_sorghum<-ggplot(T12_response_specific, aes(x=avg_log2FC, y=p_val_adj, color = as.factor(T12_response_specific[,9]), size = as.factor(T12_response_specific[,9]))) + geom_point() + theme_bw() + scale_color_manual(values=c('grey80','red','purple3','yellow3')) + xlim(c(-3,3))  + geom_vline(xintercept = 0, color = "black", size=.5) + xlab("log2 FC") + ylab("-log10 p-value") + theme(legend.position = 'none')
summary(T12_response_specific$avg_log2FC>0)
grid.arrange(T0_plot_sorghum, T12_plot_sorghum,nrow=2)

summary(T12_response_specific$V9>0)
summary(T0_response_specific$V9>0)


# Scatterplot (Extended Data)
T0_T12<-merge(T0_response, T12_response, by='row.names')
rownames(T0_T12)<-T0_T12[,1]
T0_T12<-T0_T12[,-1]
T0_T12<-subset(T0_T12, rownames(T0_T12) %in% c(rownames(T0_response_specific),rownames(T12_response_specific)))
T0_T12[,11]<-ifelse(rownames(T0_T12) %in% sorghum_lightdark_extensive_C4$gene_ID, 2, 0)
T0_T12[,12]<-ifelse(rownames(T0_T12) %in% sorghum_lightdark_extensive_calvin$gene_ID, 3, 0)
T0_T12[,13]<-ifelse(rownames(T0_T12) %in% sorghum_lightdark_extensive_light$gene_ID, 4, 0)
T0_T12[,14]<-T0_T12[,11]+T0_T12[,12]+T0_T12[,13]
T0_T12$avg_log2FC.x<-T0_T12$avg_log2FC.x*-1
T0_T12$avg_log2FC.y<-T0_T12$avg_log2FC.y*-1

ggplot(T0_T12, aes(x=avg_log2FC.x, y=avg_log2FC.y, color = as.factor(T0_T12[,14]), size = as.factor(T0_T12[,14]))) + geom_point(size=2) + theme_bw() + scale_color_manual(values=c('grey80','red4','purple3','yellow3')) + xlab("T0 M/BS") + ylab("T12 M/BS") + theme(legend.position = 'none') + geom_abline(slope=1, intercept = 0) + geom_vline(xintercept = 0, color = "black", size=.5) + geom_hline(yintercept = 0, color = "black", size=.5) 
