library(ggplot2)
library(tidyr)
library(dplyr)
out <- "docs/revision/gsea_unique_ids_v1"
x <- readRDS(file.path(out,"summary.rds"))
stopifnot(nrow(x)==32L,all(x$n_trials==1000L))
x$cell_line <- factor(x$cell_line,levels=c("TK6","MCF-7","HepaRG","Liver spheroids"))
long <- pivot_longer(x,c(ora_rate,gsea_rate,archived_rate),names_to="method",values_to="rate")
labels <- c(ora_rate="Reactome ORA",gsea_rate="Unique-ID GSEA",archived_rate="Archived duplicate-ID GSEA")
p <- ggplot(long,aes(n,rate,color=method,shape=method,linetype=method))+
  geom_hline(yintercept=.05,linetype="dashed",color="grey45",linewidth=.4)+
  geom_line(linewidth=.7)+geom_point(size=2)+facet_wrap(~cell_line,nrow=1)+
  scale_color_manual(values=c(ora_rate="#0072B2",gsea_rate="#D55E00",archived_rate="#999999"),labels=labels)+
  scale_shape_manual(values=c(ora_rate=16,gsea_rate=17,archived_rate=1),labels=labels)+
  scale_linetype_manual(values=c(ora_rate="solid",gsea_rate="longdash",archived_rate="dotted"),labels=labels)+
  scale_x_continuous(breaks=3:10)+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2))+
  labs(x="Sample size per group",y="Experiment-level false-positive rate",color=NULL,shape=NULL,linetype=NULL,
    caption="1,000 frozen mock trials per cell model and sample size; within-trial size-restricted BH p < 0.05.\nORA retains its insufficient-input gate; different method inputs prevent algorithm-only attribution.")+
  theme_bw(base_size=11)+theme(legend.position="top",panel.grid.minor=element_blank(),
    strip.background=element_rect(fill="grey95"),plot.caption=element_text(hjust=0,size=9))
ggsave(file.path(out,"reactome_ora_gsea_comparison.pdf"),p,width=11,height=4.1,device=cairo_pdf)
ggsave(file.path(out,"reactome_ora_gsea_comparison.png"),p,width=11,height=4.1,dpi=400)
c <- pivot_longer(x,c(both,ora_only,gsea_only,neither),names_to="outcome",values_to="count")
c$outcome <- factor(c$outcome,levels=c("both","ora_only","gsea_only","neither"),labels=c("Both","ORA only","GSEA only","Neither"))
p2 <- ggplot(c,aes(factor(n),count/1000,fill=outcome))+geom_col(width=.85)+facet_wrap(~cell_line,nrow=1)+
 scale_fill_manual(values=c("Both"="#009E73","ORA only"="#0072B2","GSEA only"="#E69F00","Neither"="#D9D9D9"))+
 labs(x="Sample size per group",y="Proportion of paired mock trials",fill=NULL)+
 theme_bw(base_size=11)+theme(legend.position="top",panel.grid.minor=element_blank(),strip.background=element_rect(fill="grey95"))
ggsave(file.path(out,"reactome_concordance.pdf"),p2,width=11,height=3.8,device=cairo_pdf)
ggsave(file.path(out,"reactome_concordance.png"),p2,width=11,height=3.8,dpi=400)
capture.output(sessionInfo(),file=file.path(out,"figure_sessionInfo.txt"))
