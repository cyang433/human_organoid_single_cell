###process human sc data
tmp.data1 = as.matrix(rdata1)
meta1$tag = meta1$Age_in_Weeks
tmp.form=get_form_data(tmp.data1,meta1,hvg=20000)
form.data1=tmp.form[[1]]
form.meta1=tmp.form[[2]]
form.data1 = form.data1[row.names(form.data1) %in% human.dxg,]


#relabel the time manually
sts=form.meta1$tag
for (i in c(1:length(sts))) {
        t = sts[i]
        if (t<8) {
                sts[i]=1
        } else if(t>=8 & t<10) {
                sts[i]=2
        } else if (t>=10 & t<13) {
                sts[i]=3
        } else if (t>=13 & t<16) {
                sts[i] = 4
        } else if (t>=16 & t<19) {
                sts[i] = 5
        } else if (t>=19 & t<24) {
                sts[i] = 6
        } else if (t>=24 & t<38) {
                sts[i] = 7
        } else {
                sts[i] = 8
        }
}


form.meta1$time = sts
form.meta1$ctp = unlist(sapply(form.meta1$WGCNAcluster, function(x) unlist(strsplit(x,'\\-'))[1]))
form.meta1$ctp = str_replace(form.meta1$ctp,'\\d','')
for (i in c(1:length(form.meta1$ctp))) {
        ctp = form.meta1[i,]$ctp

        if (is.na(ctp)) {
                form.meta1[i,]$ctp = 'NA'
        } else if (ctp %in% c('oRG','tRG','vRG','RG')) {
                form.meta1[i,]$ctp = 'RG'
        }else if (ctp == 'nEN') {
                form.meta1[i,]$ctp = 'EN'
        } else if (ctp == 'nIN') {
                form.meta1[i,]$ctp = 'IN'
        }
}


###process organoid sc data
time.code = data.frame('stage'=c('EB','iPSCs','Neuroectoderm','Neuroepithelium','Organoid-1M','Organoid-2M','Organoid-4M'),'time'=c(4,0,10,15,32,60,128))
times.info = meta2$Stage
for (i in 1:dim(time.code)[1]) {
        times.info[times.info==time.code[i,1]]=time.code[i,2]
}
meta2$tag = as.numeric(times.info)
#tmp.data2 = rdata2[,meta2$Line %in% c('409b2','H9')]
tmp.data2 = rdata2[,meta2$Line %in% c('H9')]
row.names(tmp.data2) = genes.info[,1]
#tmp.meta2 = meta2[meta2$Line %in% c('409b2','H9'),]
tmp.meta2 = meta2[meta2$Line %in% c('H9'),]   #!important to choose which
colnames(tmp.data2) = row.names(tmp.meta2)
tmp.form=get_form_data(tmp.data2,tmp.meta2,type='COUNTS',hvg=6000)
form.data2=tmp.form[[1]]
form.meta2=tmp.form[[2]]
#form.data2 = form.data2[row.names(form.data2) %in% org.dxg,]


##relabel the time manually
sts=form.meta2$tag
for (i in c(1:length(sts))) {
        t = sts[i]
        if (t==0) {
                sts[i]=1
        } else if (t==4) {
                sts[i]=2
        } else if (t==10) {
                sts[i]=3
        } else if (t==15) {
                sts[i] = 4
        } else if (t==32) {
                sts[i] = 5
        } else if (t==60) {
                sts[i] = 6
        } else if (t==128) {
                sts[i] = 7
        } else {
                sts[i] = 8
        }
}


form.meta2$time = sts
form.meta2$ctp = unlist(sapply(form.meta2$PredCellType, function(x) unlist(strsplit(x,'\\-'))[1]))
form.meta2$ctp = str_replace(form.meta2$ctp,'\\d','')
for (i in c(1:length(form.meta2$ctp))) {
        ctp = form.meta2[i,]$ctp
        if (is.na(ctp)) {
                form.meta2[i,]$ctp = 'NA'
        } else if (grepl('\\/',ctp)) {
                form.meta2[i,]$ctp = 'NA'
        }
}


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


