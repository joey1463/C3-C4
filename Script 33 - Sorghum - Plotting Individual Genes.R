# This script can create cell-type specific gene expression profiles for candidate genes. 

chosen_gene<-"Sobic.002G036000" # glycolate oxidase 
chosen_gene<-"Sobic.005G042000" # Rubisco 

to_plot_1<-time_average[chosen_gene,]
to_plot_2<-as.data.frame(t(to_plot_1[order(colnames(to_plot_1))]))
to_plot_2[,2]<-c(rep(0,6),rep(0.5,6),rep(1,6),rep(2,6),rep(4,6),rep(6,6),rep(8,6),rep(12,6))
to_plot_2[,3]<-c('bundle_sheath','epidermis','guard','mesophyll','phloem','xylem')    
colnames(to_plot_2)<-c('expression','time','tissue')
#to_plot_2$expression<-log(to_plot_2$expression,2)
to_plot_2<-subset(to_plot_2, to_plot_2$tissue %in% c('bundle_sheath','mesophyll'))
ggplot(data=to_plot_2, aes(x=time, y=expression, group=tissue))  + geom_point(size=2,aes(color=tissue))   + theme_bw() + theme(legend.position = "none") + ylab("normalized counts") + geom_smooth(method = "loess", se = FALSE, aes(color=tissue), size=1) + xlab("time (h)") + scale_color_manual(values=c("blue3", "green4")) +  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) + theme(axis.text=element_text(size=10),axis.title=element_text(size=15))

# Supplementary Figure
chosen_gene<-"Sobic.010G183600" #HY5
chosen_gene<-"Sobic.003G136900" #PIF
chosen_gene<-"Sobic.004G094100" # GNC
chosen_gene<-"Sobic.010G096300" # GLK1
chosen_gene<-"Sobic.003G002600" # GLK2

to_plot_1<-time_average[chosen_gene,]
to_plot_2<-as.data.frame(t(to_plot_1[order(colnames(to_plot_1))]))
to_plot_2[,2]<-c(rep(0,6),rep(0.5,6),rep(1,6),rep(12,6),rep(2,6),rep(4,6),rep(6,6),rep(8,6))
to_plot_2[,3]<-c('bundle_sheath','epidermis','guard','mesophyll','phloem','xylem')    
colnames(to_plot_2)<-c('expression','time','tissue')
ggplot(data=to_plot_2, aes(x=time, y=expression, group=tissue)) + geom_line(aes(color=tissue)) + geom_point(size=0) + theme_bw() + theme(axis.text=element_text(size=10),axis.title=element_text(size=15)) +  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))  + xlab("time (h)") + ylab("normalized expression")  
