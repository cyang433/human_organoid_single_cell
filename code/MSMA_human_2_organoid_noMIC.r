.libPaths(c("/home/che82/R/x86_64-pc-linux-gnu-library/4.1","/usr/local/lib/R/site-library","/usr/lib/R/site-library","/usr/lib/R/library" ))

#scripts for multi-stages manifold-alignment(MSMA)
library(rdist)
library(SCORPIUS)
library(Rmagic)
library(pheatmap)
library(Linnorm)
#######on alignment##########
set.seed(1000)
library(reshape2)
Sys.setenv(RETICULATE_PYTHON = "/home/che82/miniconda3/bin/python")
library(reticulate)
#set pathon path#########hcf
path_to_python<-'/home/che82/miniconda3/bin/python'
use_python(path_to_python)
library(ManiNetCluster)
library(plyr)
library(RColorBrewer)
library(stringr)

library(ComplexHeatmap)
library(circlize)
source('func.r')



##########load and preprocess the data[human(Nowaski,2017) VS organoid(kanton.2019)]#########
load('start3.RData') #has organoid cell-types tranfered accoring to human cell-types with Seurat TransferData #'transfer.ctp.r'
human.dxg = unique(unlist(read.table('nowaski.2017.genes')))
#org.dxg = unique(unlist(read.table('kanton.2019.genes')))




############select genes of interest for further analysis
## co-expression genes
sel.genes = intersect(intersect(row.names(form.data1),row.names(form.data2)),unique(all.rec$gene))
##marker genes only
mk.genes = intersect(row.names(form.data1),row.names(form.data2))
##all overlapped genes
all.genes = intersect(row.names(tmp.data1),row.names(tmp.data2))

###############select expression based genes,log transform and reorder
sel.data1 = form.data1[row.names(form.data1) %in% sel.genes,]; sel.data1 = sel.data1[!duplicated(row.names(sel.data1)),]; sel.meta1=form.meta1  #on human data
sel.data2 = form.data2[row.names(form.data2) %in% sel.genes,]; sel.data2 = sel.data2[!duplicated(row.names(sel.data2)),]; sel.meta2=form.meta2  #on organoid data
#log transform
exp1 = log(sel.data1+1)
exp2 = log(sel.data2+1)
#reorder
exp1 = exp1[order(row.names(exp1)),order(sel.meta1$time)];sel.meta1=sel.meta1[order(sel.meta1$time),]
exp2 = exp2[order(row.names(exp2)),order(sel.meta2$time)];sel.meta2=sel.meta2[order(sel.meta2$time),]




#############START aligning
####STEP 1: for single cells, form pseudo cells; for bulk, no need to form pseudo cells
tgex = t(exp1)
tmeta = sel.meta1
ps.mat = c()
ps.tag = c()
ps.belong = c()
ps.gex = c()


for(t in unique(tmeta$time)) {
        print(t)
        tmp.gex = tgex[tmeta$time==t,]
        tmp.meta = tmeta[tmeta$time==t,]
        #ps.res = get_pseudo_cells2(tmp.gex,tag=t,ctp=tmp.meta$ctp,keep.perc=(1*min(table(tmeta$time)))/dim(tmp.gex)[1],ctp.collapse.thr=.5)
        ps.res = get_pseudo_cells3(tmp.gex,tag=t,ctp=tmp.meta$ctp,keep.perc = 70/dim(tmp.gex)[1],ctp.collapse.thr=.5)
        ps.mat = rbind(ps.mat,ps.res[[1]])
        ps.tag = rbind(ps.tag,ps.res[[2]])
	ps.belong = c(ps.belong,ps.res[[3]])
	ps.gex = rbind(ps.gex,ps.res[[4]])

}
ps.mat1=ps.mat
ps.time1=ps.tag[,1]
ps.ctp1 = ps.tag[,2]
ps.belong1 = ps.belong
ps.gex1 = ps.gex
ps.sc.time1 = tmeta$time
ps.sc.ctp1 = data.frame('orig_ctp'=tmeta$WGCNAcluster, 'ctp'=tmeta$ctp)
#######form pseudo cells 2#####
tgex = t(exp2)
tmeta = sel.meta2
ps.mat = c()
ps.tag = c()
ps.belong = c()
ps.gex = c()

for(t in unique(tmeta$time)) {
        print(t)
        tmp.gex = tgex[tmeta$time==t,]
        tmp.meta = tmeta[tmeta$time==t,]
        #ps.res = get_pseudo_cells2(tmp.gex,tag=t,ctp=tmp.meta$ctp,keep.perc=(1*min(table(tmeta$time)))/dim(tmp.gex)[1],ctp.collapse.thr=.5)
        ps.res = get_pseudo_cells3(tmp.gex,tag=t,ctp=tmp.meta$ctp,keep.perc = 71/dim(tmp.gex)[1],ctp.collapse.thr=.5)
        ps.mat = rbind(ps.mat,ps.res[[1]])
        ps.tag = rbind(ps.tag,ps.res[[2]])
	ps.belong = c(ps.belong,ps.res[[3]])
	ps.gex = rbind(ps.gex,ps.res[[4]])
}
ps.mat2=ps.mat
ps.time2=ps.tag[,1]
ps.ctp2 = ps.tag[,2]
ps.belong2 = ps.belong
ps.gex2 = ps.gex
ps.sc.time2 = tmeta$time
ps.sc.ctp2 = data.frame('orig_ctp' = tmeta$NewCellType, 'ctp'=tmeta$ctp)


#########step 2: MSMA---just replace the cell-type,no need to re align
#algn_res = runMSMA(ps.mat1,ps.mat2,method='cor')
algn_res = runMSMA_cor(ps.mat1,ps.mat2)


#######downstream analysis---not in the main
df2=algn_res[[3]]
df2$time = c(ps.time1,ps.time2)
df2$ctp = c(ps.ctp1,ps.ctp2)
df2$name = c(row.names(ps.mat1),row.names(ps.mat2))
########calculate pairwise distances between cells after 2nd MA#####
pair_dist = apply(df2[df2$data=='sample1',c(3:5)],1,function(x) {
        d = apply(df2[df2$data=='sample2',c(3:5)],1,function(y) eu.dist(x,y))
})
row.names(pair_dist)=row.names(ps.mat2)
colnames(pair_dist)=row.names(ps.mat1)
sim_dist = 1/(1+pair_dist)

#############FOR cells+cells similarity bi-clustering#######
sim_mat = t(sim_dist)
write.table(sim_mat,file='dist.csv',row.names=T,col.names=T,sep=',',quote=F)


#######run python bi-cluster##############!!!!!!!!!!!!!!!!!!!!!!!!
##python biclust.py #!!!!!!!!!!!
####then,python bicluster###########!!!!!!!!!!!!!!!!!!
#cls1 = as.vector(unlist(read.table('cls.cells1.csv',sep=',')))
#cls2 = as.vector(unlist(read.table('cls.cells2.csv',sep=',')))
#df2$cls = c(cls1,cls2)
#
####kmeans?
##kk = kmeans(df2[,c(3:5)],7)
##df2$cls = kk$cluster
#
#
##df2$avg_time = 0#reorder the clusters based on the average sampling time of each cluster
##cls_time = c()
##for (c in unique(df2$cls)) {
##	tmp = df2[df2$cls==c,]
##	m = mean(tmp$time)
##	df2[df2$cls==c,]$avg_time = m
##	cls_time = rbind(cls_time,c(c,m))
##}
##cls_time = data.frame(cls_time[order(cls_time[,2]),])
##names(cls_time)=c('cls','avg_time')
##cls_time$clsr = c(1:dim(cls_time)[1])
##df2$clsr=mapvalues(df2$cls,from=cls_time$cls,to=cls_time$clsr)
##
##df2_1 = df2[df2$data=='sample1',]; 
##df2_2 = df2[df2$data=='sample2',]
##
##cls_sim_mat = sim_mat[order(df2_1$clsr),order(df2_2$clsr)]
##df3 = rbind(df2_1[order(df2_1$clsr),],df2_2[order(df2_2$clsr),])
##
##library(RColorBrewer)
##cols0 = rainbow(9);names(cols0)=paste('cluster',c(1:9)) #for clusters
##cols1 = brewer.pal(n=7,name='Blues') ;names(cols1)=paste('time',c(1:7)) #for time points-human
##cols2 = brewer.pal(n=7,name='YlOrBr') ;names(cols2)=paste('time',c(1:7)) #for time points-organoid
##annot01 = rowAnnotation(time_human=paste('time',df3[df3$data=='sample1',]$time),cls_human=paste('cluster',df3[df3$data=='sample1',]$clsr),col=list(cls_human=cols0,time_human=cols1))
##annot02 = columnAnnotation(time_org=paste('time',df3[df3$data=='sample2',]$time),cls_org=paste('cluster',df3[df3$data=='sample2',]$clsr),col=list(cls_org=cols0,time_org=cols2))
###pdf('algn2.dist.pdf',height=10,width=10)
##Heatmap(cls_sim_mat,,show_row_names=F,show_column_names=F,cluster_rows=F,cluster_columns=F,left_annotation=annot01,top_annotation=annot02,row_names_gp=gpar(fontsize=2),column_names_gp=gpar(fontsize=2))
###dev.off()
##
##
################summarize cluster into a table for downstream analysis#######
#####pcell to scell#
##ps.df1 = df2[df2$data=='human',]
##row.names(ps.df1) = ps.df1$name
##cs.df1 = ps.df1[ps.belong1,]           #correspond to ps.gex1
##row.names(cs.df1) = row.names(ps.gex1)
##cs.df1$ctp = ps.sc.ctp1$ctp
##cs.df1$orig_ctp = ps.sc.ctp1$orig_ctp
##
##ps.df2 = df2[df2$data=='organoid',]
##row.names(ps.df2) = ps.df2$name
##cs.df2 = ps.df2[ps.belong2,]           #correspond to ps.gex1
##row.names(cs.df2) = row.names(ps.gex2)
##cs.df2$ctp = ps.sc.ctp2$ctp
##cs.df2$orig_ctp = ps.sc.ctp2$orig_ctp
##
##
#########association of cell type on each cluster[hypergeometric test]####
##cal.hyper<-function(ctps_list,cls_list,unq.cls = c(1:9)) {
##	ctps = sort(unique(ctps_list))
##	cls = sort(unique(cls_list))   
##	pvs = matrix(NA,nrow=length(unq.cls),ncol=length(ctps))
##	row.names(pvs)=unq.cls;colnames(pvs) = ctps
##	for (cl in unq.cls) {
##		for (ctp in ctps) {
##			#print (paste(cl,ctps))
##			x = length(intersect(which(ctps_list==ctp),which(cls_list==cl)))
##			m = length(which(ctps_list==ctp))
##			n = length(ctps_list)-m
##			k = length(which(cls_list==cl))
##			if (k==0) {
##				pvs[which(unq.cls ==cl),which(ctps==ctp)]=NA
##			}
##			else {
##				p = phyper(x,m,n,k,lower.tail=F,log.p=F)
##				pvs[which(unq.cls ==cl),which(ctps==ctp)] = p
##			}
##		}
##	}
##	return(pvs)
##}
##
##
##
########association of new cell types on each cluster---grouped celltypes to higher orders#######
##high.ctp<-function(ctp) {
##	nctp = ctp
##	if (ctp %in% c('oRG','RG-div1','RG-div2','RG-early','tRG','vRG')) {
##		nctp = 'RG'
##	} else if (ctp %in% c('nEN-early1','nEN-early2','nEN-late')) {
##		nctp = 'nEN'
##	} else if (ctp %in% c('EN-V1-1','EN-V1-2','EN-V1-3')) {
##		nctp = 'EN-V1'
##	} else {
##		nctp = str_replace(nctp,'\\d','')
##	} 
##	return(nctp)
##}
##cs.df1$ctp = sapply(cs.df1$orig_ctp,high.ctp)
##cs.df2$ctp = sapply(cs.df2$orig_ctp,high.ctp)
##
##hp3 = cal.hyper(cs.df1$ctp,cs.df1$clsr)
##hp4 = cal.hyper(cs.df2$ctp,cs.df2$clsr)
##
##hp3 = hp3[,order(colnames(hp3))]
##hp4 = hp4[,order(colnames(hp4))]
##hp3 = hp3[,nchar(colnames(hp3))>0]
##hp4 = hp4[,nchar(colnames(hp4))>0]
##
##
##cols = brewer.pal(9,'Set1')
##names(cols) = row.names(hp1)
##annot.cls = rowAnnotation(cls = row.names(hp1),col=list(cls=cols))
##ht3=Heatmap(-1*log10(hp3+1e-10),cluster_rows=F,cluster_columns=F,name='Human',col=colorRamp2(c(2:10),c(brewer.pal(9, "Blues"))))
##ht4=Heatmap(-1*log10(hp4+1e-10),cluster_rows=F,cluster_columns=F,name='Organoid',right_annotation=annot.cls,col=colorRamp2(c(2:10),c(brewer.pal(9, "Blues"))))
##draw(ht3+ht4, ht_gap = unit(1, "cm"))
##
##
##
##
#####single source 3D---organoid#######color the dots based on cluster[re-ordered cluster:clsr] information
###pdf('human_vs_humanOrg_organoid.pdf')
##time.cols = brewer.pal(9,'Set1')
##res = data.frame(df2[df2$data=='organoid',])
##library(plot3D)
##s3d<-scatter3D(x=res[,3],y=res[,4],z=res[,5],col='grey',bg=time.cols[as.numeric(mapvalues(res$clsr,names(table(res$clsr)),c(1:5)))],lwd=1,pch=24,
##colkey=F,theta = 30, phi = 0,cex=1.5,
##xlim=c(min(df2$Val0),max(df2$Val0)),
##ylim=c(min(df2$Val1),max(df2$Val1)),
##zlim=c(min(df2$Val2),max(df2$Val2)))
##legend("top", legend = levels(as.factor(res$data)), pch = c(16, 17),inset = -0.1, xpd = TRUE, horiz = TRUE)
##legend("bottom", legend = levels(as.factor(res$clsr)), col = time.cols[c(1:9)],pch=16,inset = -0.1, xpd = TRUE, horiz = TRUE)
###dev.off()
##
##
##
##
##
