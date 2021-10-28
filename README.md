# Single cell alignment
## Summary
We integrated recently published developmental scRNA-seq data of both human brains and organoids and aligned single cells across brains and organoids. Aligned and unaligned cells show certain developmental trajectories and form various cell clusters. We also annotated cell-type and developmental information for those clusters and identified differentially expressed genes in the clusters. Further gene set enrichment analyses reveal organoid-specific developmental functions and pathways. 
## Software Requirements
The analysis is based on R 4.0. Users should install the following packages prior to using the code:
```R
install.packages(c('rdist','SCORPIUS','Rmagic','pheatmap','Linnorm','reshape2','reticulate','ManiNetCluster','plyr','RColorBrewer','stringr','ComplexHeatmap','circlize'))
```
Besides,some functions for processing data are also needed for our project.
```R
source('../src/func.r')
```
## Demo
### Step 0: Data preprocessing
human(Nowaski,2017) VS organoid(kanton.2019)
```R
load('start3.RData')
human.dxg = unique(unlist(read.table('nowaski.2017.genes')))
```
We use the uploaded code to process the data,including relabel the time manually,selecting the gene of interest for analysis and 
```R
source('process_data.r')
```
Next we start to align.

### Step 1: Pseudo-cell formation

Unlike bulk analysis, single cell analysis must first construct pseudo cells, because single cells are highly random.
```R
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
```
### Step 2: Manifold alignment
```R
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
```
### Step 3: Use python to bi-cluster
```
python biclust.py
```
Read clustering results:
```R
cls1 = as.vector(unlist(read.table('cls.cells1.csv',sep=',')))
cls2 = as.vector(unlist(read.table('cls.cells2.csv',sep=',')))
df2$cls = c(cls1,cls2)
```
Re-order clusters for visualization:
```R
clsr2 = rep(0,dim(df2)[1])
clsr2[df2$cls==1]=1
clsr2[df2$cls==3]=2
clsr2[df2$cls==4]=3
clsr2[df2$cls==5]=4
clsr2[df2$cls==2]=5
df2$clsr = clsr2
```

### Step 4: Visualization
```R
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
```
