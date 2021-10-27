#this is the official analysis scripts for single human VS organoid after alignment!!!!!!!!!!!!!!!!!!!!

#######run python bi-cluster##############!!!!!!!!!!!!!!!!!!!!!!!!
#python biclust.py
####then,python bicluster###########!!!!!!!!!!!!!!!!!!
#bicluster to 5 clusters??

cls1 = as.vector(unlist(read.table('cls.cells1.csv',sep=',')))
cls2 = as.vector(unlist(read.table('cls.cells2.csv',sep=',')))
df2$cls = c(cls1,cls2)

##kmeans
#kk = kmeans(df2[,c(3:5)],7)
#df2$cls = kk$cluster


###########################re-order clusters based on time
#df2$avg_time = 0#reorder the clusters based on the average sampling time of each cluster
#cls_time = c()
#for (c in unique(df2$cls)) {
#	tmp = df2[df2$cls==c,]
#	m = mean(tmp$time)
#	df2[df2$cls==c,]$avg_time = m
#	cls_time = rbind(cls_time,c(c,m))
#}
#cls_time = data.frame(cls_time[order(cls_time[,2]),])
#names(cls_time)=c('cls','avg_time')
#cls_time$clsr = c(1:dim(cls_time)[1])
#df2$clsr=mapvalues(df2$cls,from=cls_time$cls,to=cls_time$clsr)

###############re-order clusters for visualization
clsr2 = rep(0,dim(df2)[1])
clsr2[df2$cls==1]=1
clsr2[df2$cls==3]=2
clsr2[df2$cls==4]=3
clsr2[df2$cls==5]=4
clsr2[df2$cls==2]=5
df2$clsr = clsr2


#########################scatter plot#############
###single source 3D---organoid#######color the dots based on cluster[re-ordered cluster:clsr] information
pdf('human_vs_humanOrg_human2.pdf')
time.cols = brewer.pal(5,'Set1')
res = data.frame(df2[df2$data=='sample2',])
library(plot3D)
s3d<-scatter3D(x=res[,3],y=res[,4],z=res[,5],col='grey',bg=time.cols[as.numeric(mapvalues(res$clsr,names(table(res$clsr)),c(1:5)))],lwd=1,pch=21,
colkey=F,theta = 30, phi = 0,cex=1.5,
xlim=c(min(df2$Val0),max(df2$Val0)),
ylim=c(min(df2$Val1),max(df2$Val1)),
zlim=c(min(df2$Val2),max(df2$Val2)))
legend("top", legend = levels(as.factor(res$data)), pch = c(16, 17),inset = -0.1, xpd = TRUE, horiz = TRUE)
legend("bottom", legend = levels(as.factor(res$clsr)), col = time.cols[c(1:5)],pch=16,inset = -0.1, xpd = TRUE, horiz = TRUE)
dev.off()



##############summarize cluster into a table for downstream analysis#######
###pcell to scell#
ps.df1 = df2[df2$data=='human',]
row.names(ps.df1) = ps.df1$name
cs.df1 = ps.df1[ps.belong1,]           #correspond to ps.gex1
row.names(cs.df1) = row.names(ps.gex1)
cs.df1$ctp = ps.sc.ctp1$ctp
cs.df1$orig_ctp = ps.sc.ctp1$orig_ctp

ps.df2 = df2[df2$data=='organoid',]
row.names(ps.df2) = ps.df2$name
cs.df2 = ps.df2[ps.belong2,]           #correspond to ps.gex1
row.names(cs.df2) = row.names(ps.gex2)
cs.df2$ctp = ps.sc.ctp2$ctp
cs.df2$orig_ctp = ps.sc.ctp2$orig_ctp


#######association of cell type on each cluster[hypergeometric test]####
cal.hyper<-function(ctps_list,cls_list,unq.cls = c(1:5)) {
        ctps = sort(unique(ctps_list))
        cls = sort(unique(cls_list))
        pvs = matrix(NA,nrow=length(unq.cls),ncol=length(ctps))
        row.names(pvs)=unq.cls;colnames(pvs) = ctps
        for (cl in unq.cls) {
                for (ctp in ctps) {
                        #print (paste(cl,ctps))
                        x = length(intersect(which(ctps_list==ctp),which(cls_list==cl)))
                        m = length(which(ctps_list==ctp))
                        n = length(ctps_list)-m
                        k = length(which(cls_list==cl))
                        if (k==0) {
                                pvs[which(unq.cls ==cl),which(ctps==ctp)]=NA
                        }
                        else {
                                p = phyper(x,m,n,k,lower.tail=F,log.p=F)
                                pvs[which(unq.cls ==cl),which(ctps==ctp)] = p
                        }
                }
        }
        return(pvs)
}

hp1 = cal.hyper(cs.df1$orig_ctp,cs.df1$clsr)
hp2 = cal.hyper(cs.df2$orig_ctp,cs.df2$clsr)

hp1 = hp1[,order(colnames(hp1))]
hp2 = hp2[,order(colnames(hp2))]

pdf('human_vs_humanOrg.ctp.all.pdf',width=15)
cols = brewer.pal(5,'Set1')
names(cols) = row.names(hp1)
annot.cls = rowAnnotation(cls = row.names(hp1),col=list(cls=cols))
ht1=Heatmap(-1*log10(hp1+1e-10),cluster_rows=F,cluster_columns=F,name='Human',col=colorRamp2(c(2:7),c(brewer.pal(9, "Greys"))[2:7]),column_names_rot=90)
ht2=Heatmap(-1*log10(hp2+1e-10),cluster_rows=F,cluster_columns=F,name='Organoid',right_annotation=annot.cls,col=colorRamp2(c(2:7),c(brewer.pal(9, "Greys"))[2:7]),column_names_rot=90)
draw(ht1+ht2, ht_gap = unit(1, "cm"))
dev.off()



######association of new cell types on each cluster---grouped celltypes to higher orders#######
high.ctp<-function(ctp) {
	nctp = ctp
	if (ctp %in% c('oRG','RG-div1','RG-div2','RG-early','tRG','vRG')) {
		nctp = 'RG'
	} else if (ctp %in% c('nEN-early1','nEN-early2','nEN-late')) {
		nctp = 'nEN'
	} else if (ctp %in% c('EN-V1-1','EN-V1-2','EN-V1-3')) {
		nctp = 'EN-V1'
	} else {
		nctp = str_replace(nctp,'\\d','')
	} 
	return(nctp)
}
cs.df1$ctp = sapply(cs.df1$orig_ctp,high.ctp)
cs.df2$ctp = sapply(cs.df2$orig_ctp,high.ctp)

hp3 = cal.hyper(cs.df1$ctp,cs.df1$clsr)
hp4 = cal.hyper(cs.df2$ctp,cs.df2$clsr)

hp3 = hp3[,order(colnames(hp3))]
hp4 = hp4[,order(colnames(hp4))]
hp3 = hp3[,nchar(colnames(hp3))>0]
hp4 = hp4[,nchar(colnames(hp4))>0]

pdf('human_vs_humanOrg.ctp.high_noMIC.pdf',width=10)  #########whether to use MICROGLIA CASE
cols = brewer.pal(5,'Set1')
names(cols) = row.names(hp3)
annot.cls = rowAnnotation(cls = row.names(hp3),col=list(cls=cols))
ht3=Heatmap(-1*log10(hp3+1e-10),cluster_rows=F,cluster_columns=F,name='Human',col=colorRamp2(c(2:7),c(brewer.pal(9, "Greys"))[2:7]),column_names_rot=45)
ht4=Heatmap(-1*log10(hp4+1e-10),cluster_rows=F,cluster_columns=F,name='Organoid',right_annotation=annot.cls,col=colorRamp2(c(2:7),c(brewer.pal(9, "Greys"))[2:7]),column_names_rot=45)
draw(ht3+ht4, ht_gap = unit(1, "cm"))
dev.off()



#######cluster---time, association
#count on single cells within each bicluster, then on time p-value
#######on time distribution[hypergeometric test]####
cal.time.hyper<-function(times_list,cls_list,unq.cls = c(1:5)) {
        times = sort(unique(times_list))
        cls = sort(unique(cls_list))
        pvs = matrix(NA,nrow=length(unq.cls),ncol=length(times))
        row.names(pvs)=unq.cls;colnames(pvs) = times
        for (cl in unq.cls) {
                for (time in times) {
                        x = length(intersect(which(times_list==time),which(cls_list==cl)))
                        m = length(which(times_list==time))
                        n = length(times_list)-m
                        k = length(which(cls_list==cl))
                        if (k==0) {
                                pvs[which(unq.cls ==cl),which(times==time)]=NA
                        } else {
                                p = phyper(x,m,n,k,lower.tail=F,log.p=F)
                                pvs[which(unq.cls ==cl),which(times==time)] = p
                        }
                }
        }
        return(pvs)
}

time.hp1 = cal.time.hyper(cs.df1$time,cs.df1$clsr)
time.hp2 = cal.time.hyper(cs.df2$time,cs.df2$clsr)

p.adj1 = matrix(p.adjust(time.hp1,method='BH'),ncol=7);row.names(p.adj1)=row.names(time.hp1);colnames(p.adj1)=colnames(time.hp1)
p.adj2 = matrix(p.adjust(time.hp2,method='BH'),ncol=7);row.names(p.adj2)=row.names(time.hp2);colnames(p.adj2)=colnames(time.hp2)

#pheatmap(-1*log10(time.hp1+1e-100),cluster_rows=F,cluster_cols=F)
#pheatmap(-1*log10(time.hp2+1e-100),cluster_rows=F,cluster_cols=F)
pdf('cluster-time.correspond.human.pdf')
pheatmap(1*(p.adj1<0.05),cluster_rows=F,cluster_cols=F,color=c('white','grey'))
dev.off()

pdf('cluster-time.correspond.org.pdf')
pheatmap(1*(p.adj2<0.05),cluster_rows=F,cluster_cols=F,color=c('white','grey'))
dev.off()


##################################################################################
#######on genes#############
#############selected pcells sign back to single cells---FOR LATER, MAY NEED TO PERFORM ANALYSIS ON SINGLE CELLS && genes other than the selected.genes
############extract other genes gxp from the aligned cells############
mk.ps.gex1 = t(log(form.data1[mk.genes,row.names(ps.gex1)]+1))
mk.ps.gex2 = t(log(form.data2[mk.genes,row.names(ps.gex2)]+1))

mk.ps.mat1 = sapply(ps.df1$name,function(t) {
                tmp = mk.ps.gex1[ps.belong1==t,]
                if (is.null(dim(tmp))) {
                        return(tmp)
                } else {
                        return(colMeans(mk.ps.gex1[ps.belong1==t,]))
                }
})
mk.ps.mat2 = sapply(ps.df2$name,function(t) {
                tmp = mk.ps.gex2[ps.belong2==t,]
                if (is.null(dim(tmp))) {
                        return(tmp)
                } else {
                        return(colMeans(mk.ps.gex2[ps.belong2==t,]))
                }
})


#############SELECT on genes of interest##################
scale.ps.gex1 = scale(mk.ps.gex1)
scale.ps.gex2 = scale(mk.ps.gex2)

#genes on biclusters --- from gordon alignment||||||||filter 1
library(presto)
cls_res1 = wilcoxauc(t(scale.ps.gex1),cs.df1$clsr)
cls_res2 = wilcoxauc(t(scale.ps.gex2),cs.df2$clsr)
val_res11 = cls_res1[cls_res1$padj<=0.05 & cls_res1$logFC > .3,]
val_res12 = cls_res2[cls_res2$padj<=0.05 & cls_res2$logFC > .3,]


#another filter, based on time|||||||||||||||||||filter 2
time_res1 = wilcoxauc(t(scale.ps.gex1),cs.df1$time)
time_res2 = wilcoxauc(t(scale.ps.gex2),cs.df2$time)

val_res21 = time_res1[time_res1$padj<=0.05 & time_res1$logFC > .3,]
val_res22 = time_res2[time_res2$padj<=0.05 & time_res2$logFC > .3,]

#presto p-value based on cell.type|||||||||||||||||||DO NOT filter, just for visualization
ctp_res1 = wilcoxauc(t(scale.ps.gex1),cs.df1$ctp)
ctp_res2 = wilcoxauc(t(scale.ps.gex2),cs.df2$ctp)



##select correspondence based on p-value
corr1 = which(p.adj1<=0.05,arr.ind=T);colnames(corr1)=c('bicluster','time')
corr2 = which(p.adj2<=0.05,arr.ind=T);colnames(corr2)=c('bicluster','time')

library(ComplexHeatmap)
library(circlize)

###filter based on bicluster & time correspondence---get genes
sig.genes1 = c()
for (k in c(1:dim(corr1)[1])) {
        grp1 = corr1[k,1]
        grp2 = corr1[k,2]
        print(paste(k,grp1,grp2,' '))
        info1 = val_res11[val_res11$group==grp1,];row.names(info1) = info1$feature
        info2 = val_res21[val_res21$group==grp2,];row.names(info2) = info2$feature
        ovp_genes = intersect(info1$feature,info2$feature)
        tmp.info = data.frame('gene'=ovp_genes,'cluster'=rep(grp1,length(ovp_genes)),'time'=rep(grp2,length(ovp_genes)),'time.padj'=info1[ovp_genes,]$padj,'cluster.padj'=info2[ovp_genes,]$padj)
        sig.genes1 = rbind(sig.genes1,tmp.info)
}

sig.genes2 = c()
for (k in c(1:dim(corr2)[1])) {
        grp1 = corr2[k,1]
        grp2 = corr2[k,2]
        info1 = val_res12[val_res12$group==grp1,];row.names(info1) = info1$feature
        info2 = val_res22[val_res22$group==grp2,];row.names(info2) = info2$feature
        ovp_genes = intersect(info1$feature,info2$feature)
        tmp.info = data.frame('gene'=ovp_genes,'cluster'=rep(grp1,length(ovp_genes)),'time'=rep(grp2,length(ovp_genes)),'time.padj'=info1[ovp_genes,]$padj,'cluster.padj'=info2[ovp_genes,]$padj)
        sig.genes2 = rbind(sig.genes2,tmp.info)
}

write.table(sig.genes1,'human.time_cluster.genes.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(sig.genes2,'organoid.time_cluster.genes.tsv',col.names=T,row.names=F,sep='\t',quote=F)
####################################################
######draw the gene overlapping across seleted bicluster/time##############
int.cls1 = 3;int.time1 = 3 #human---selected cluster and time
int.cls2 = 3;int.time2 = 6 #organoid---select cluster and time

sel.genes1 = sig.genes1[sig.genes1$cluster %in% int.cls1 & sig.genes1$time %in% int.time1,]$gene
sel.genes2 = sig.genes2[sig.genes2$cluster %in% int.cls2 & sig.genes2$time %in% int.time2,]$gene

HO_con = intersect(sel.genes1,sel.genes2)
H_spec = sel.genes1[!(sel.genes1 %in% HO_con)]
O_spec = sel.genes2[!(sel.genes2 %in% HO_con)]
#find genes ctp enrichment p-value---save as a matrix
library(tidyr)
gen.mat<-function(glist,pst.res) {#given selected genes and presto res, generate a p-value matrix
        tmp.res = pst.res[pst.res$feature %in% glist & nchar(pst.res$group)>0,c('feature','group','padj')]
        tmp.res$feature=as.factor(tmp.res$feature)
        tmp.res$group=as.factor(tmp.res$group)
        mat = pivot_wider(tmp.res, names_from = group, values_from = padj)
        return(mat)
}

#HO_con_ctp = gen.mat(HO_con,ctp_res1)
#!!!use human brain as the gold standard
HO_con_ctp = data.frame(gen.mat(HO_con,ctp_res1));row.names(HO_con_ctp)=HO_con_ctp[,1];HO_con_ctp[,1]<-c()
H_spec_ctp = data.frame(gen.mat(H_spec,ctp_res1));row.names(H_spec_ctp)=H_spec_ctp[,1];H_spec_ctp[,1]<-c()
O_spec_ctp = data.frame(gen.mat(O_spec,ctp_res1));row.names(O_spec_ctp)=O_spec_ctp[,1];O_spec_ctp[,1]<-c()


####NOW, generate network for each gene list
library(igraph)

pvs = HO_con_ctp #change here

#way 1: use minimal p-value for each gene
#recs = as.matrix(data.frame('gene'=row.names(pvs),'ctp'=colnames(pvs)[apply(pvs,1,function(x) which(x==min(x))[1])]))  
#way 2: number of edges eq. to number of gene
library(dplyr)
library(igraph)
len = dim(pvs)[1]
thr = sort(unlist(pvs))[len]

idx = which(pvs<thr,arr.ind=T)
full.nodes = row.names(pvs)

recs = data.frame('gene'=row.names(idx),'ctp'=colnames(pvs)[idx[,2]])
full.nodes = data.frame(c(full.nodes,unique(recs$ctp)))
#recs[is.na(recs[,2]),2] = 'NA' #color other genes as grey

#g<-graph_from_data_frame(d=recs,vertices=full.nodes,directed=F)
#plot(g,vertex.size=5)
library("ggnetwork")

ctps = unique(c(colnames(H_spec_ctp),colnames(O_spec_ctp),colnames(HO_con_ctp)))

grs<-graph_from_data_frame(d=recs,vertices=full.nodes,directed=F)
V(grs)$size = 1
dgrs = data.frame('node'=names(V(grs)),'degree'=centralization.degree(grs)$res)
dgrs$color = dgrs$degree;dgrs$color[dgrs$degree==0]='grey';dgrs$color[dgrs$degree>0]='blue';dgrs$color[dgrs$node %in% ctps]='red'

ggdf = ggnetwork(grs, layout = layout_with_fr(grs), cell.jitter = 0)
#color manual
ggdf$group = sapply(ggdf$name,function(x) dgrs$color[dgrs$node==x])

pdf('HO_con.net.pdf',width=10)
ggplot(ggdf, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", curvature = 0.1, size = 0.5, alpha = 0.7)+
  geom_nodes(aes(x = x, y = y, size=size, color= factor(group)), alpha = 0.7) + scale_color_manual(values=c('blue','black','red'))+
  geom_nodelabel_repel(aes(label = name,color=factor(group)), size=6,alpha=1)+
  theme_blank()# + theme(legend.position = "none")
dev.off()



##check the distrbution of pvalues???????????MAY BE NEED TO DO this
#pvs = O_spec_ctp
#hist.data = hist(-1*log10(as.matrix(pvs)),breaks=100)
#hist.data$counts = log10(hist.data$counts)
#hist.data$counts[is.infinite(hist.data$counts)] = 0
#barplot(hist.data$counts)








###presto p-values of genes on each cluster,time and cell-type for heatmap
#generate the p-value matrix
#######following is for the 3 matrix visualization:for full version, check 'sig.genes.from.human.vs.organoid.algn.r' & 'tmp3.r'[time names]#######################
library(tidyr)
gen.mat<-function(glist,pst.res) {#given selected genes and presto res, generate a p-value matrix
        tmp.res = pst.res[pst.res$feature %in% glist & nchar(pst.res$group)>0,c('feature','group','pval')]
        tmp.res$feature=as.factor(tmp.res$feature)
        tmp.res$group=as.factor(tmp.res$group)
        mat = pivot_wider(tmp.res, names_from = group, values_from = pval)
        return(mat)
}

#visualize on human specific genes
mat11=data.frame(gen.mat(unique(sig.genes1$gene),cls_res1));row.names(mat11)=mat11[,1];mat11[,1]<-c()
mat12=data.frame(gen.mat(unique(sig.genes1$gene),time_res1));row.names(mat12)=mat12[,1];mat12[,1]<-c()
mat13=data.frame(gen.mat(unique(sig.genes1$gene),ctp_res1));row.names(mat13)=mat13[,1];mat13[,1]<-c()
log.mat11 = as.matrix(-1*log10(mat11+1e-50))
log.mat12 = as.matrix(-1*log10(mat12+1e-50))
log.mat13 = as.matrix(-1*log10(mat13+1e-50))
ht1=Heatmap(log.mat11,cluster_columns=F,show_row_names=F,col=colorRamp2(c(median(log.mat11),max(log.mat11)),c('#FFFFCC','red')))
ht2=Heatmap(log.mat12,cluster_columns=F,show_row_names=F,col=colorRamp2(c(median(log.mat12),max(log.mat12)),c('#FFFFCC','blue')))
ht3=Heatmap(log.mat13,cluster_columns=F,show_row_names=F,,col=colorRamp2(c(median(log.mat13),max(log.mat13)),c('#FFFFCC','black')))

ht1+ht2+ht3

#visualize on organoid genes
mat21=data.frame(gen.mat(unique(sig.genes2$gene),cls_res2));row.names(mat21)=mat21[,1];mat21[,1]<-c()
mat22=data.frame(gen.mat(unique(sig.genes2$gene),time_res2));row.names(mat22)=mat22[,1];mat22[,1]<-c()
mat23=data.frame(gen.mat(unique(sig.genes2$gene),ctp_res2));row.names(mat23)=mat23[,1];mat23[,1]<-c()
log.mat21 = as.matrix(-1*log10(mat21+1e-100))
log.mat22 = as.matrix(-1*log10(mat22+1e-100))
log.mat23 = as.matrix(-1*log10(mat23+1e-100))
ht4=Heatmap(log.mat21,cluster_columns=F,show_row_names=F,col=colorRamp2(c(median(log.mat21),max(log.mat21)),c('#FFFFCC','red')))
ht5=Heatmap(log.mat22,cluster_columns=F,show_row_names=F,col=colorRamp2(c(median(log.mat22),max(log.mat22)),c('#FFFFCC','blue')))
ht6=Heatmap(log.mat23,cluster_columns=F,show_row_names=F,,col=colorRamp2(c(median(log.mat23),max(log.mat23)),c('#FFFFCC','black')))

ht4+ht5+ht6
############################################################################















