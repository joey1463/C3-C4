# This script will construct a Volcano plot that describes the genes found differentially partitioned between mesophyll and bundle sheath cell types
# at either the 0h time point or 12 h time point. 

setwd("~/Desktop/Salk/Project - Light/R Files")
load("L3_rice_T0_response_10X.RData")
load("L3_rice_T12_response_10X.RData")
T0_response_specific<-subset(T0_response, T0_response$p_val_adj<0.05)
T12_response_specific<-subset(T12_response, T12_response$p_val_adj<0.05)

# Volcano Plot
rice_lightdark_extensive_calvin<-subset(rice_lightdark_extensive, rice_lightdark_extensive$category %in% c('calvin_cycle','photorespiration'))
rice_lightdark_extensive_light<-subset(rice_lightdark_extensive, rice_lightdark_extensive$category %in% c('light_reaction'))

# T0
T0_response_specific[,6]<-ifelse(rownames(T0_response_specific) %in% rice_lightdark_extensive_calvin$gene.ID, 2, 0)
T0_response_specific[,7]<-ifelse(rownames(T0_response_specific) %in% rice_lightdark_extensive_light$gene.ID, 3, 0)
T0_response_specific[,8]<-T0_response_specific[,6]+T0_response_specific[,7]

T0_response_specific$p_val_adj<--log(T0_response_specific$p_val_adj,10)
T0_plot_rice<-ggplot(T0_response_specific, aes(x=avg_log2FC, y=p_val_adj, color = as.factor(T0_response_specific[,8]), size = as.factor(T0_response_specific[,8]))) + geom_point() + theme_bw() + scale_color_manual(values=c('grey80','purple3','yellow3'))  + xlim(c(-3,3)) + geom_vline(xintercept = 0, color = "black", size=.5) + ylim(c(0,60)) + xlab("log2 FC") + ylab("-log10 p-value") + theme(legend.position = 'none')

# T12
T12_response_specific[,6]<-ifelse(rownames(T12_response_specific) %in% rice_lightdark_extensive_calvin$gene.ID, 2, 0)
T12_response_specific[,7]<-ifelse(rownames(T12_response_specific) %in% rice_lightdark_extensive_light$gene.ID, 3, 0)
T12_response_specific[,8]<-T12_response_specific[,6]+T12_response_specific[,7]
T12_response_specific$p_val_adj<--log(T12_response_specific$p_val_adj,10)
T12_plot_rice<-ggplot(T12_response_specific, aes(x=avg_log2FC, y=p_val_adj, color = as.factor(T12_response_specific[,8]), size = as.factor(T12_response_specific[,8]))) + geom_point() + theme_bw() + scale_color_manual(values=c('grey80','purple3','yellow3','hotpink3')) + xlim(c(-3,3)) + geom_vline(xintercept = 0, color = "black", size=.5) + ylim(c(0,60)) + xlab("log2 FC") + ylab("-log10 p-value") + theme(legend.position = 'none')

grid.arrange(T0_plot_rice, T12_plot_rice, nrow=2)

summary(T12_response_specific$V8>0)
summary(T0_response_specific$V8>0)


# Scatterplot (Extended Data)
T0_T12<-merge(T0_response, T12_response, by='row.names')
rownames(T0_T12)<-T0_T12[,1]
T0_T12<-T0_T12[,-1]
T0_T12<-subset(T0_T12, rownames(T0_T12) %in% c(rownames(T0_response_specific),rownames(T12_response_specific)))
T0_T12[,11]<-ifelse(rownames(T0_T12) %in% rice_lightdark_extensive_calvin$gene.ID, 2, 0)
T0_T12[,12]<-ifelse(rownames(T0_T12) %in% rice_lightdark_extensive_light$gene.ID, 3, 0)
T0_T12[,13]<-T0_T12[,11]+T0_T12[,12]
T0_T12$avg_log2FC.x<-T0_T12$avg_log2FC.x*-1
T0_T12$avg_log2FC.y<-T0_T12$avg_log2FC.y*-1

ggplot(T0_T12, aes(x=avg_log2FC.x, y=avg_log2FC.y, color = as.factor(T0_T12[,13]), size = as.factor(T0_T12[,13]))) + geom_point(size=2) + theme_bw() + scale_color_manual(values=c('grey80','purple3','yellow3')) + xlab("T0 M/BS") + ylab("T12 M/BS") + theme(legend.position = 'none') + geom_abline(slope=1, intercept = 0) + geom_vline(xintercept = 0, color = "black", size=.5) + geom_hline(yintercept = 0, color = "black", size=.5) + xlim(c(-2,2)) + ylim(c(-2,2))
