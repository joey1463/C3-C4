# This script can create cell-type specific gene expression profiles for candidate genes. 

chosen_gene<-"LOC-Os03g57220" # glycolate oxidase
chosen_gene<-"LOC-Os12g17600" # Rubisco subunit

to_plot_1<-time_average[chosen_gene,]
to_plot_2<-as.data.frame(t(to_plot_1))
to_plot_3<-to_plot_2[str_detect(rownames(to_plot_2), "mesophyll"),][1:8]
to_plot_4<-to_plot_2[str_detect(rownames(to_plot_2), "bundle"),][1:8]
output<-c(to_plot_3, to_plot_4)
to_plot_5<-as.data.frame(cbind(output,c(rep(c(0,0.5,1,12,2,4,6,8),2))))
to_plot_5[,3]<-c(rep('mesophyll',8),rep('bundle_sheath',8))  
colnames(to_plot_5)<-c('expression','time','tissue')
#to_plot_5$expression<-log(to_plot_5$expression, 2)
plot_5<-ggplot(data=to_plot_5, aes(x=time, y=expression, group=tissue)) 
ggplot(data=to_plot_5, aes(x=time, y=expression, group=tissue))  + geom_point(size=2,aes(color=tissue))   + theme_bw() + theme(legend.position = "none") + ylab("normalized counts") + geom_smooth(method = "loess", se = FALSE, aes(color=tissue), size=1) + xlab("time (h)") + scale_color_manual(values=c("blue3", "green4")) +  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) + theme(axis.text=element_text(size=12),axis.title=element_text(size=16))


# HY5s and PIFs (Supplementary) (starch analysis deleted)

chosen_gene<-"LOC-Os02g10860" # HY5
chosen_gene<-"LOC-Os12g41650" # PIF
chosen_gene<-"LOC-Os02g12790" # GNC
chosen_gene<-"LOC-Os06g24070" # GLK1 
chosen_gene<-"LOC-Os01g13740" # GLK2 

to_plot_1<-time_average[chosen_gene,]
to_plot_2<-as.data.frame(t(to_plot_1))
to_plot_2[,2]<-c(rep(0,6),rep(0.5,6),rep(1,6),rep(12,6),rep(2,6),rep(4,6),rep(6,6),rep(8,6))
to_plot_2[,3]<-c('bundle_sheath','epidermis','guard','mesophyll','phloem','xylem')  
colnames(to_plot_2)<-c('expression','time','tissue')
ggplot(data=to_plot_2, aes(x=time, y=expression, group=tissue)) + geom_line(aes(color=tissue)) + geom_point(size=0) + theme_bw() + theme(axis.text=element_text(size=10),axis.title=element_text(size=15)) + xlab("time (h)") + ylab("normalized expression")   +  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))+ theme(legend.position="none")
