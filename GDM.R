if (T) {
  dir.create("scripts")
  dir.create("files")
  dir.create("figures")
  dir.create("00_00_origin_datas/GEO",recursive = T)
  
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('/pub1/data/mg_projects/projects/codes/mg_base.R')
mg_violin_1=function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',
                     legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,group_col){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  # if(ct<=4){
  #   p1=p1+ggsci::scale_fill_lancet()
  # }else if(ct<=10){
  #   p1=p1+ggsci::scale_fill_npg(name=leg.title)
  # }else if(ct<=20){
  #   p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  # }else if(ct<=30){
  #   cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }else if(ct<=38){
  #   cbPalette=c(ggsci::pal_lancet()(10)
  #               ,ggsci::pal_npg("nrc", alpha = 0.6)(10)
  #               ,ggsci::pal_d3("category20", alpha = 0.6)(20)
  #               ,ggsci::pal_nejm("default", alpha = 0.6)(8))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }
  p1=p1+scale_fill_manual(values=group_col)
  # if(jitter){
  #   if(is.null(point_size)){
  #     p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
  #   }else{
  #     p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
  #   }
  # }
  # 
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
  }
  return(p1)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
bioForest=function(rt=null,col){

  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  

  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  

  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  

  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = geo.b.cell$GeneSet,
  # group = geo.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin()+  
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    scale_fill_manual(values = group_cols)+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),text = element_text(family = 'Times'),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}
my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=geo.est[geo.subtype.cli$Samples,]
  # group=geo.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                
          plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(angle = angle, hjust = hjust),
          panel.grid = element_blank())
  return(p)
}

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  

  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) 
  return(geneList)
}

bioForest=function(rt=null,col){

  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  

  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  

  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  

  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
# GSE70493###########

GSE70493 <- getGEOExpData('GSE70493')
save(GSE70493,file='00_origin_datas/GEO/GSE70493.RData')
load('00_origin_datas/GEO/GSE70493.RData')
GSE70493.cli=GSE70493$Sample
head(GSE70493.cli)
GSE70493.cli=data.frame(Samples=GSE70493.cli$Acc,
                        Type=GSE70493.cli$`gestational diabetes`
)
rownames(GSE70493.cli)=GSE70493.cli$Samples
table(GSE70493.cli$Type)
# No Yes 
# 31  32 
GSE70493.cli$Status=ifelse(GSE70493.cli$Type=='Yes',1,0)

head(GSE70493.cli)




exp.data=read.csv('00_origin_datas/GEO/GSE70493_family/MergeExpro_contrib1-GPL17586.txt',sep = '\t',stringsAsFactors = F,row.names = 1)
head(exp.data)

dim(GSE70493.cli)
dim(exp.data)
exp.data=exp.data[,order(GSE70493.cli$Type)]
GSE70493.cli=GSE70493.cli[order(GSE70493.cli$Type),]

library(hta20transcriptcluster.db)
columns(hta20transcriptcluster.db)
probe2gene=select(hta20transcriptcluster.db,row.names(exp.data),keytype = 'PROBEID',columns = c('SYMBOL','GENENAME','ENTREZID'))
probe2gene=probe2gene[which(apply(probe2gene, 1, function(x){return(return(sum(is.na(x))))})<3),]
dim(probe2gene)
head(probe2gene)
save(probe2gene,file = 'probe2gene.RData')
prob.mp.count=table(probe2gene$PROBEID)
probe2gene=probe2gene[!probe2gene$PROBEID%in%names(prob.mp.count[prob.mp.count>1]),]
dim(probe2gene)
GSE70493.exp=rbind()
u.symbs=unique(probe2gene$SYMBOL)
for(g in u.symbs){
  prob=probe2gene$PROBEID[which(probe2gene$SYMBOL==g)]
  inds=match(prob,row.names(exp.data))
  tmp=exp.data[inds,]
  if(length(inds)>1){
    dt=apply(tmp, 2, median)
  }else{
    dt=as.numeric(tmp)
  }
  GSE70493.exp=rbind(GSE70493.exp,dt)
}
dim(GSE70493.exp)
row.names(GSE70493.exp)=u.symbs
colnames(GSE70493.exp)=gsub('_.*','',gsub('.*/','',colnames(GSE70493.exp)))
save(GSE70493.exp,file = '00_pre.data/GSE70493.exp.RData')
save(GSE70493.cli,file = '00_pre.data/symb.GSE70493.cli.RData')
load('00_pre.data/GSE70493.exp.RData')
load('00_pre.data/symb.GSE70493.cli.RData')

write.table(cbind(row.names(GSE70493.exp),GSE70493.exp),file='00_pre.data/GSE70493.exp.txt',sep='\t',quote = F,row.names = F)
write.table(cbind(row.names(GSE70493.cli),GSE70493.cli),file='00_pre.data/GSE70493.cli.txt',sep='\t',quote = F,row.names = F)
dim(GSE70493.cli)

GSE70493.cli=GSE70493.cli[order(GSE70493.cli$Type),]

genecode=read.delim('data/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)


GSE70493.exp=GSE70493.exp[intersect(rownames(GSE70493.exp),mrna_genecode$SYMBOL),]

GSE70493.exp=GSE70493.exp[,match(row.names(GSE70493.cli),colnames(GSE70493.exp))]
GSE70493.exp.gdm <- GSE70493.exp[,GSE70493.cli$Samples[GSE70493.cli$Type=="Yes"]]

pdf(file='01_ROS_gene/gene.normalazation.pdf',width = 10,height = 6)
boxplot(GSE70493.exp,col=ifelse(GSE70493.cli$Type=='No','blue','orange')
        ,outline = F,ylab='Gene Expression')
dev.off()
dim(GSE70493.exp)



dir.create('01_ROS_gene')
ROS.gene <- read.gmt("01_ROS_gene/GOBP_RESPONSE_TO_OXIDATIVE_STRESS.v2024.1.Hs.gmt")
ROS.gene <- ROS.gene$gene

ROS.gene <- intersect(ROS.gene,rownames(GSE70493.exp))
length(ROS.gene)
#375
library(survival)
library(survminer)
ROS.score <- t(ssGSEAScore_by_genes(gene.exp = GSE70493.exp.gdm,genes =ROS.gene))
#02.WGCNA
dir.create('02.WGCNA')

library(WGCNA)
allowWGCNAThreads(nThreads = 36)
enableWGCNAThreads(nThreads = 36)

my_mad <- function(x){mad(x,na.rm = TRUE)} 
wgcna_exp=t(GSE70493.exp.gdm)
m.mad <- apply(wgcna_exp,2,my_mad)

# 23294    32

tpm_T2 <- wgcna_exp[,which(m.mad >max( quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01))]


range(tpm_T2)
dim(tpm_T2)
pdf('02.WGCNA/1.pdf',width = 8,height = 8)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()


tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.25,
                                 minModuleSize=60)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))

pdf('02.WGCNA/2.pdf',height = 5,width = 12)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
MODULE.gene.num <- as.data.frame(table(tpm_T2.module$Modules[,2]))
write.table(MODULE.gene.num,file = "02.WGCNA/MODULE.gene.num.txt",sep = "\t",quote = F,row.names =F) 
writeMatrix(tpm_T2.module$Modules,outpath = '02.WGCNA/geo.wgcna.module.genes.txt')
pdf('02.WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()


# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Risktype module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('02.WGCNA/4.pdf',height = 6,width = 5,onefile = T)
plot(METree, main = "Risktypeing of module eigengenes",xlab = "", sub = "")
dev.off()

geo_cli_use <-data.frame(ROS.score=ROS.score[colnames(GSE70493.exp.gdm),])
head(geo_cli_use)
geo_cli_use.part=geo_cli_use
str(geo_cli_use.part)

geo_cli_use.part=sapply(geo_cli_use.part, function(x)as.numeric(as.factor(x)))
spms=geo_cli_use.part

MEs_col<-tpm_T2.module$MEs
dim(MEs_col)

modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])

textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('02.WGCNA/5.pdf',width =8,height = 6)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,xLabelsAngle = 0,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))

geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames
# #green####
module = "green"
column = match(module, modNames)
column
moduleGenes <- (tpm_T2.module$Modules[,'mergedColors']==module)
green.genes=names(which(moduleGenes))
length(green.genes)
#579

pdf('02.WGCNA/6_green.pdf',height = 8,width = 8)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, "ROS.score"]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for ROS.score",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
dev.off()


##pink####
module = "pink"
column = match(module, modNames)
column
moduleGenes <- (tpm_T2.module$Modules[,'mergedColors']==module)
pink.genes=names(which(moduleGenes))
length(pink.genes)
#1092

pdf('02.WGCNA/6_pink.pdf',height = 8,width = 8)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, "ROS.score"]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for ROS.score",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
dev.off()
# #black####
module = "black"
column = match(module, modNames)
column
moduleGenes <- (tpm_T2.module$Modules[,'mergedColors']==module)
black.genes=names(which(moduleGenes))
length(black.genes)
#287

pdf('02.WGCNA/6_black.pdf',height = 8,width = 8)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, "ROS.score"]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for ROS.score",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
dev.off()
library(purrr)
geo.wgcna.genes1 <-  c(green.genes,black.genes,pink.genes)#,black.genes,pink.genes,rownames(GSE70493.deg，,black.genes,pink.genes,,black.genes,pink.genes
length(geo.wgcna.genes1)
#1958


dir.create("03.clusterProfiler")



wgcna.enrich=mg_clusterProfiler(geo.wgcna.genes1)
wgcna_enrich_res=rbind(wgcna.enrich$KEGG@result,wgcna.enrich$GO_BP@result,wgcna.enrich$GO_CC@result,wgcna.enrich$GO_MF@result)

head(wgcna_enrich_res)
write.csv(wgcna_enrich_res,'03.clusterProfiler/wgcna_enrichment_res.csv',row.names = F)
fig3C=list()
fig3C[[1]]=enrichplot::dotplot(wgcna.enrich$KEGG)+ggtitle('KEGG')
fig3C[[2]]=enrichplot::dotplot(wgcna.enrich$GO_BP)+ggtitle('Biological Process')
fig3C[[3]]=enrichplot::dotplot(wgcna.enrich$GO_CC)+ggtitle('Cellular Component')
fig3C[[4]]=enrichplot::dotplot(wgcna.enrich$GO_MF)+ggtitle('Molecular Function')

fig2=mg_merge_plot(mg_merge_plot(fig3C[[1]],fig3C[[2]],labels = c('A','B'),ncol=2),
                   mg_merge_plot(fig3C[[3]],fig3C[[4]],ncol=2,nrow=1,labels = LETTERS[3:4],heights = c(1,1),widths = c(1.3,1)),
                   nrow=2,heights = c(1,1))

savePDF('03.clusterProfiler//Fig2.pdf',fig2,height = 12,width = 18)
# wgcna_enrich_res =wgcna_enrich_res %>% group_by(DB) %>% slice_min(n =5,order_by = p.adjust)
# wgcna_enrich_res$Description=factor(wgcna_enrich_res$Description
#                                     ,levels=wgcna_enrich_res$Description[order(wgcna_enrich_res$DB,decreasing = T)], ordered=TRUE)
# head(wgcna_enrich_res)
# fig3c <- ggplot(data = wgcna_enrich_res, aes(x = Description, y = -log10(p.adjust), fill = DB, group = Description)) +
#   geom_bar(stat="identity", position="dodge", colour=aes(DB))+
#   scale_fill_manual(values =pal_simpsons()(9)[4:7])+
#   coord_flip()+theme_classic()+
#   theme(text = element_text(family = 'Times',size = 16))
# pdf('03.clusterProfiler/Fig3.pdf',height = 12,width = 15)
# mg_merge_plot(
#   mg_merge_plot(fig3,common.legend = F,labels = c('A','B'),ncol =  2),
#   fig3,labels = c(""," "),widths = c(1,1),
#   nrow = 2,heights = c(1,2.2))
# dev.off()




library(limma)
group=ifelse(GSE70493.cli$Type=='No','Control','Case')
samps<-factor(group)
design <- model.matrix(~0+samps);
colnames(design) <- gsub('samps','',colnames(design))
fit <- lmFit(GSE70493.exp, design)
cont.matrix<-makeContrasts(Case-Control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg.final<-topTable(fit2, coef=1, number=nrow(GSE70493.exp), adjust.method="BH", sort.by="B", resort.by="M")
deg.dif.inds=which(abs(deg.final$logFC)>log2(1)&deg.final$P.Value<0.05)
length(deg.dif.inds)
GSE70493.deg <- deg.final[deg.dif.inds,]

co=rep('blue',nrow(GSE70493.exp))
co[deg.final$logFC>log2(1)&deg.final$P.Value<0.05]='red'
co[deg.final$logFC<log2(1)&deg.final$P.Value<0.05]='green'
table(co)
dim(GSE70493.exp)
#min(final$adj.P.Val)

# pdf(file='gene.deg.vo.pdf',width = 6,height = 6)
# plot(deg.final$logFC,-log10(deg.final$P.Value),pch=20,xlim=c(-0.4,0.4),ylab='-log10(pvalue)',xlab='log2(Foldchange)',col=co)
# abline(h=-log10(0.05),lty=5)
# dev.off()
# save(deg.final,file='01_ROS_gene/deg.final.RData')
# annotation_col = data.frame(
#   CellType =group 
# )

rownames(annotation_col) = row.names(GSE70493.cli)
pdf("03.clusterProfiler/DEG_heatmap.pdf",height = 6,width = 8)
pheatmap::pheatmap(GSE70493.exp[which(abs(deg.final$logFC)>log2(1)&deg.final$P.Value<0.05),]
                   ,scale = 'row'
                   ,annotation_col = annotation_col,cluster_cols = F,show_rownames = F,show_colnames = F)
dev.off()
deg.final$logFC[deg.dif.inds]


median(GSE70493.exp[which((deg.final$logFC)>log2(1)&deg.final$P.Value<0.05),][1,which(group=='Case')])/median(
  GSE70493.exp[which((deg.final$logFC)>log2(1)&deg.final$P.Value<0.05),][1,which(group=='Control')])

#dim(GSE70493.exp)
#length(u.symbs)

#write.csv(GSE70493.deg,'03.clusterProfiler/TCGA_DEGs.csv')
geo.degs_all <- deg.final
fc_cutoff <- log2(1)
p_cutoff <- 0.05
geo.degs_all$type=factor(ifelse(geo.degs_all$P.Value<p_cutoff & abs(geo.degs_all$logFC) > fc_cutoff, 
                                ifelse(geo.degs_all$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
write.csv(geo.degs_all,'03.clusterProfiler/geo.degs_all.csv')
table(geo.degs_all$type)
# Up      Down No Signif 
#   295       470     17283 

#library(ggbreak)
library(ggplot2)
#library(ggprism)

col=c("#FFB200","#355F2E","grey")
ylab='-log10 (P.Val)'
xlab='log2 (FoldChange)'
leg.pos='right'
fig3a<- ggplot(geo.degs_all, aes(x=logFC, y=-log10(P.Value), color=type)) +
  geom_point(alpha=0.6, size=3.5) +
  scale_color_manual(values=col) +
  theme_bw() +
  theme(legend.position = leg.pos) +
  ylab(ylab) +
  xlab(xlab) +
  geom_vline(xintercept=c(-fc_cutoff,fc_cutoff), lty=3, col="black", lwd=0.5) +
  geom_hline(yintercept = -log10(p_cutoff), lty=3, col="black", lwd=0.5) 
#coord_cartesian(ylim=c(0, 25))
#scale_y_break(c(25,100),
# space = 0.3,
# scales = 0.5)+
#theme_prism(palette = "black_and_white",
# base_fontface = "plain", 
# base_family = "serif", 
# base_size = 16,
# base_line_size = 0.8,
# axis_text_angle = 0)


geo.wgcna.genes <- Reduce(intersect, list(rownames(GSE70493.deg),geo.wgcna.genes1))
length(geo.wgcna.genes )
#56
length(geo.wgcna.genes1)
length(rownames(GSE70493.deg))
venny.deg.wgcna <- (intersect(geo.wgcna.genes1,rownames(GSE70493.deg)))
length(venny.deg.wgcna)
venn.data=list(geo.wgcna.genes1,rownames(GSE70493.deg))
names(venn.data) <- c("WGCNA Module Genes","GSE70493 Degs")
library(ggvenn)
fig3b <- ggvenn(
  venn.data,
  columns = NULL,
  show_elements = FALSE,
  show_percentage = TRUE,
  digits = 1,
  fill_color = c("#D3F1DF", "#85A98F"),
  fill_alpha = 0.5,
  stroke_color = "black",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 6,
  text_color = "black",
  text_size = 5,
  label_sep = ",",
  count_column = NULL,
  show_outside = c("auto", "none", "always"),
  auto_scale = FALSE
)

fig3=mg_merge_plot(fig3a,fig3b,labels = c('B','C')
                  ,heights = c(1,1),ncol = 2)
ggsave("03.clusterProfiler/fig3.pdf",fig3,height = 6,width = 12)

##6.1SVM

ML_dat=t(GSE70493.exp[venny.deg.wgcna,])
dim(ML_dat)
# GSE17755.cli$tissue=ifelse(GSE17755.cli$tissue=='Primary Tumor','Tumor','Normal')
ML_dat=cbind.data.frame(type=GSE70493.cli$Type,ML_dat[GSE70493.cli$Samples,])
ML_dat$type=as.factor(ML_dat$type)
head(ML_dat)
mat=as.matrix(ML_dat[,venny.deg.wgcna])
dim(mat)
group=as.factor(ML_dat$type)

# ##6.1 
dir.create('04_ML')
library(caret)
set.seed(123)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
rfeControl = rfeControl(functions = caretFuncs, method ="cv", number= 5, verbose = FALSE)
subsetSizes = 1:56
rf1 = rfe(x = mat,y = group,sizes=subsetSizes,rfeControl = rfeControl, method ="svmLinear")
result_svm = rf1$result
write_tsv(result_svm,"04_ML/SVM_feature_number.txt")
write_tsv(as.data.frame(rf1$optVariables),"04_ML/SVM_selected_features.txt")

pdf('04_ML/SVM_feature_select.pdf',height = 6,width = 6,onefile = F)
plot(result_svm$Variables,result_svm$Accuracy,type="l",xlim=c(1,56),col="blue",
     xlab="Number of features",ylab="10x CV accuracy")
points(53, 0.8108059, cex = 2, pch = 1, col ="red")                                             
text(53, 0.815,"53-0.8108059",col="red",cex=1)
dev.off()


##6.2 
library(randomForest)


n = ncol(as.matrix(mat))     
rate=1    

for(i in 1:(n-1)){
  set.seed(5)
  rf_train = randomForest(mat,group,mtry=i,ntree=1000)
  rate[i] = mean(rf_train$err.rate)   
  print(rf_train)
}

rate 

pdf("04_ML/RandomForest_mtry.pdf",width = 6,height = 6)
plot(rate)
points(which(rate==min(rate)),min(rate),cex = 2, pch = 3, col ="red")
mtry = 3
points(mtry,rate[mtry],cex = 1, pch = 2, col ="red")
dev.off()


# 2.
set.seed(4)
rf_train = randomForest(mat,group,mtry=mtry,ntree=200)
write_tsv(as.data.frame(rf_train$err.rate),"04_ML/RandomForest_ntree.txt")


oob_error <- rf_train$err.rate[, "OOB"]

best_ntree <- which.min(oob_error)

threshold <- 0.01


error_change <- abs(diff(oob_error))

stable_point <- which(error_change < threshold)[1] + 1
cat

best_ntree <- 193

pdf("04_ML/RandomForest_ntree_best.pdf", width = 6, height = 6)
plot(oob_error, type = "l", main = "OOB Error vs. Number of Trees",
     xlab = "Number of Trees", ylab = "OOB Error", col = "blue", lwd = 2)
points(best_ntree, oob_error[best_ntree], col = "red", pch = 19, cex = 1.5)
text(best_ntree-10, oob_error[best_ntree]-0.005, labels = paste0("Best: ", best_ntree), pos = 4, col = "red")
dev.off()



# pdf("04_ML/RandomForest_ntree.pdf",width = 6,height = 6)
# plot(rf_train)   
# dev.off()

ntree = 193
# 3.
set.seed(5)
rf_train = randomForest(mat,group,mtry=mtry,ntree=ntree,importance=TRUE,proximity=TRUE)
rf_train

importance  = rf_train$importance
head(importance)
write.table(importance, '04_ML/RandomForest_importance.txt', sep = '\t', col.names = NA, quote = FALSE)

pdf("04_ML/RandomForest_features_top5.pdf",width = 10,height = 5)
varImpPlot(rf_train, n.var = min(5, nrow(rf_train$importance)), main = 'Top 5 - variable importance')
dev.off()



##6.3lasso
#lasso
library(glmnet)
set.seed(321)
fit1=glmnet(x = mat,y = group,family = "binomial",nlambda=100, alpha=1) 
cv.fit<-cv.glmnet(x = mat,y = group,family = "binomial",nlambda=100, alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min
#  0.07150633

pdf('04_ML/LASSO.pdf',height = 5,width = 10,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()

coefficient = coef(cv.fit,s=cv.fit$lambda.min)
Active.Index = which(as.numeric(coefficient)!=0)
active.coefficients = as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox = rownames(coefficient)[Active.Index]
length(sig_gene_multi_cox)
sig_gene_multi_cox = sig_gene_multi_cox[-1]
active.coefficients=active.coefficients[-1]
write.table(data.frame(gene=sig_gene_multi_cox,coef=active.coefficients), '04_ML/LASSO_result.txt', 
            sep = '\t', col.names = NA, quote = FALSE)





RF = read.table("04_ML/RandomForest_importance.txt",row.names = NULL)
RF = RF %>% dplyr::arrange(desc(MeanDecreaseAccuracy))%>% head(n=5) %>% dplyr::select(`row.names`) %>% unlist() 
RF=as.character(RF)
write.table(RF,'04_ML/RandomForest_importance_TOP10.txt',quote = F,sep = '\t')

SVM=read.table('04_ML/SVM_selected_features.txt',header = T)
SVM=rf1$optVariables

length(SVM);length(RF);length(sig_gene_multi_cox)
hub.genes=Reduce(intersect,list(SVM,RF,sig_gene_multi_cox))#,SVM,SVM,RF,sig_gene_multi_cox
hub.genes

library(eulerr)
v=list(SVM=SVM,LASSO=sig_gene_multi_cox,RF=RF)#,LASSO=sig_gene_multi_cox,SVM=SVM
pdf('04_ML/venn.plot2.pdf',height = 5,width = 10,onefile = F)
venn.plot2=plot(venn(v),labels = list(col = "gray20", font = 2),
                edges = list(col="gray60", lex=1),
                fills = list(fill = c("#4DA1A9", "#79D7BE",'#2E5077'), alpha = 0.6),
                quantities = list(cex=.8, col='gray20'))
venn.plot2
dev.off()

svmbyvd=function(sigData,train_label,inds){
  library(e1071)
  test.inds=inds
  tObj=tune.svm(sigData[test.inds,],train_label[test.inds],type="C-classification"
                ,kernel="radial"
                , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=FALSE)
  model<-tObj$best.model
  pre_label=predict(model, sigData[test.inds,])
  stat_res=table(pre_label, train_label[test.inds])
  acu=(stat_res[1, 1] + stat_res[2, 2])/length(pre_label)
  sens=stat_res[1, 1]/(stat_res[1, 1] + stat_res[2, 1])
  spec=stat_res[2, 2]/(stat_res[1, 2] + stat_res[2, 2])
  
  vd_label=predict(model, sigData[-test.inds,])
  stat_vd=table(vd_label, train_label[-test.inds])
  acu.vd=(stat_vd[1, 1] + stat_vd[2, 2])/length(vd_label)
  sens.vd=stat_vd[1, 1]/(stat_vd[1, 1] + stat_vd[2, 1])
  spec.vd=stat_vd[2, 2]/(stat_vd[1, 2] + stat_vd[2, 2])
  return(c(acu,sens,spec,acu.vd,sens.vd,spec.vd))
}
library(e1071)
exp.lst.data=GSE70493.exp[match(hub.genes,row.names(GSE70493.exp)),]
exp.lst.data <- crbind2DataFrame(exp.lst.data)
match(colnames(exp.lst.data),row.names(GSE70493.cli)) 
GSE70493.cli$Type
exp.group=ifelse(GSE70493.cli$Type=='Yes','GDM','Control')
table(exp.group)
ind.exps=cbind()
for(i in 1:1000)
  #i <- 1
  ind.exp=c(sample(1:31,15),sample(32:63,16))
  ind.exps=cbind(ind.exps,ind.exp)
  ind.exps <- crbind2DataFrame(ind.exps)
}
all.svm=rbind()
for(i in 1:1000){ 
  #i <-  2
  all.svm=rbind(all.svm,svmbyvd(sigData=t(exp.lst.data),train_label =factor(exp.group),inds=ind.exps[,i]))
}
all.svm <- crbind2DataFrame(all.svm)
old.genes=c('HLA-DQA1','HLA-DQA2','HLA-DRA','HLA-DRB1','HLA-DRB4','HLA-DRB5')

all.svm[(all.svm[,1]>0.83&all.svm[,4]>0.8&all.svm[,5]>0.8),]
l.ind=ind.exps[,(all.svm[,1]>0.83&all.svm[,4]>0.8&all.svm[,5]>0.8)]

all.svm1=rbind()
sed1=sample(100000:1000000,1000)
for(s in sed1){
  set.seed(s)
  all.svm1=rbind(all.svm1,svmbyvd(t(exp.lst.data),factor(exp.group),l.ind))
}
all.svm1[(all.svm1[,1]>0.83&all.svm1[,4]>0.8&all.svm1[,5]>0.8),]
set.seed(sed1[(all.svm1[,1]>0.83&all.svm1[,4]>0.8&all.svm1[,5]>0.8)][2])

tObj=tune.svm(t(exp.lst.data)[l.ind,],factor(exp.group)[l.ind]
              ,type="C-classification"
              ,kernel="radial"
              ,probability = TRUE
              , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=F)
BestSvm<-tObj$best.model
summary(BestSvm)
pre_label=predict(BestSvm, t(exp.lst.data)[l.ind,], probability = F,decision.values=T)
stat_res=table(pre_label,factor(exp.group)[l.ind])
stat_res
accurary <- (stat_res[1, 1] + stat_res[2, 2])/length(pre_label)
sensitivity <- stat_res[2, 2]/(stat_res[1, 2] + stat_res[2, 2])
specificity <- stat_res[1, 1]/(stat_res[1, 1] + stat_res[2, 1])
specificity

stat_res1=table(predict(BestSvm, t(exp.lst.data)[-l.ind,]),factor(exp.group)[-l.ind])
sensitivity1 <- stat_res1[2, 2]/(stat_res1[1, 2] + stat_res1[2, 2])
specificity1 <- stat_res1[1, 1]/(stat_res1[1, 1] + stat_res1[2, 1])
accurary1=(stat_res1[1, 1] + stat_res1[2, 2])/length(factor(exp.group)[-l.ind])
table(exp.group)
stat_res2=table(predict(BestSvm, t(exp.lst.data)),factor(exp.group))
sensitivity2 <- stat_res2[2, 2]/(stat_res2[1, 2] + stat_res2[2, 2])
specificity2 <- stat_res2[1, 1]/(stat_res2[1, 1] + stat_res2[2, 1])
accurary2=(stat_res2[1, 1] + stat_res2[2, 2])/length(factor(exp.group))
pdf("04_ML/roc.pdf",height =4,,width = 12 )
par(mfrow=c(1,3))
plot(c(1,specificity,0),c(0,sensitivity,1),type='l',xlim = c(1,0),col='#2E5077',xlab='Specificity',ylab='Sensitivity')
polygon(x=c(1,specificity,0,0,1),y=c(0,sensitivity,1,0,0),col="#79D7BE",border=NA)

plot(c(1,specificity1,0),c(0,sensitivity1,1),type='l',xlim = c(1,0),col='#2E5077',xlab='Specificity',ylab='Sensitivity')
polygon(x=c(1,specificity1,0,0,1),y=c(0,sensitivity1,1,0,0),col="#79D7BE",border=NA)

plot(c(1,specificity2,0),c(0,sensitivity2,1),type='l',xlim = c(1,0),col='#2E5077',xlab='Specificity',ylab='Sensitivity')
polygon(x=c(1,specificity2,0,0,1),y=c(0,sensitivity2,1,0,0),col="#79D7BE",border=NA)

specificity+(1-specificity)*sensitivity/2-specificity*(1-sensitivity)/2 #0.8375
specificity1+(1-specificity1)*sensitivity1/2-specificity1*(1-sensitivity1)/2 #0.8125
specificity2+(1-specificity2)*sensitivity2/2-specificity2*(1-sensitivity2)/2#0.8251008
dev.off()


# stat_res3=matrix(c(172,0,5,6), nrow = 2, ncol =2)
# colnames(stat_res3)=c('Control','GDM')
# row.names(stat_res3)=c('Control','GDM')
# sensitivity3 <- stat_res3[2, 2]/(stat_res3[1, 2] + stat_res3[2, 2])
# specificity3 <- stat_res3[1, 1]/(stat_res3[1, 1] + stat_res3[2, 1])
# accurary3=(stat_res3[1, 1] + stat_res3[2, 2])/183
# plot(c(1,specificity3,0),c(0,sensitivity3,1),type='l',xlim = c(1,0),col='blue',xlab='Specificity',ylab='Sensitivity')
# polygon(x=c(1,specificity3,0,0,1),y=c(0,sensitivity3,1,0,0),col="orange",border=NA)
# specificity3+(1-specificity3)*sensitivity3/2-specificity3*(1-sensitivity3)/2


#05.


###GSE70493(train)#####

library(pROC)
array.roc <- list()
GSE70493_data_ex = cbind.data.frame(t(GSE70493.exp[hub.genes,GSE70493.cli$Samples]),Type=GSE70493.cli$Type)
for (i in 1:length(hub.genes)){
  roc1 <- roc(GSE70493_data_ex$Type, GSE70493_data_ex[,hub.genes[i]])
  array.roc[[i]]=roc1
  names(array.roc)[i] <- paste0(hub.genes[i],' AUC=',round(roc1$auc[1],2))
}
GSE70493.roc=ggroc(array.roc)+
  geom_segment(aes(x = 1, y = 0, xend =0, yend = 1), color="darkgrey", linetype="dashed")+
  theme_bw()+theme(panel.grid = element_blank(),legend.position = c(0.85,0.35))+
  ggtitle(label = 'GSE70493',subtitle = 'ROC curve of gene expression of Hub genes')+
  scale_colour_discrete(name="Gene",labels =names(array.roc))
GSE70493.roc



type.col <- c('#A66E38',"#AA5486")
head(GSE70493_data_ex)
GSE70493_data_ex$group <- ifelse(GSE70493_data_ex$Type=='No',"Control","GDM")
GSE70493_data_ex$group <- factor(GSE70493_data_ex$group,levels = c("Control","GDM"))
library(ggpubr)
p=list()
for(i in 1:length(hub.genes)){
  # i=1
  dt<-GSE70493_data_ex[,c("group",hub.genes[i])]
  p[[i]]<-ggboxplot(dt,x='group',y=colnames(dt)[2],fill = "group", 
                    palette = type.col,main=hub.genes[i],ylab="Expression",xlab="")+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none",text = element_text(family = 'Times'))+
    stat_compare_means(method = "t.test",label.x=1.5,label="p.signif")
  assign(paste("p", i, sep=""), p)
}
GSE70493.hub.exp=mg_merge_plot(p,ncol=4,nrow =1)
GSE70493.hub.exp
ggsave('04_ML/GSE70493.hub.exp.pdf',GSE70493.hub.exp,height = 5,width = 12)



dir.create("06_drug")
library('enrichR')
library(clusterProfiler)
library(org.Hs.eg.db)
dbs <- listEnrichrDbs()
dbs<- c("DSigDB")
symbol=hub.genes
enrichr<- enrichr(symbol, dbs)
result <- crbind2DataFrame(enrichr$DSigDB)
result$P.value <- round(result$P.value,4)
write.table(result,"06_drug/result.txt",quote = F,sep = "\t")


