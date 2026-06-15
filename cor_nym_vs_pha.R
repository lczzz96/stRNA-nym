pha = readRDS("/Users/wg/pha/S1/filtered_feature_bc_matrix/pha.rds")
pha_cluster = read.table("/Users/wg/pha/S1/filtered_feature_bc_matrix/S1.map.40",header = T,row.names = 1)

pha = readRDS("/Users/wg/pha/S3/level_matrix/level_13/pha.rds")
pha_cluster = read.table("/Users/wg/pha/S3/S3.map.40",header = T,row.names = 1)

pha = NormalizeData(pha)
pha = FindVariableFeatures(pha, nfeatures = 1000, selection.method = "vst")
pha = ScaleData(pha)
pha$cluster =  0
pha@meta.data[rownames(pha_cluster),'cluster'] = pha_cluster$Cluster

pha_t = subset(pha,cluster %in% c('38','39'))
pha_t@meta.data[pha_t$cluster == '39','cluster'] = 'p';pha_t@meta.data[pha_t$cluster == '38','cluster'] = 's'
pha_t = subset(pha,cluster %in% c('1','36'))
pha_t@meta.data[pha_t$cluster == '1','cluster'] = 'p';pha_t@meta.data[pha_t$cluster=='36','cluster'] = 's'
pha_t = subset(pha,cluster %in% c('7','9','31','32','33','21','25'))
pha_t@meta.data[pha_t$cluster %in% c('7','9','31','32','33'),'cluster'] = 'p';pha_t@meta.data[pha_t$cluster %in% c('21','25'),'cluster'] = 's'
pha_t@active.ident = factor(pha_t$cluster)
pha_marker = FindAllMarkers(pha_t,min.pct = 0.25,only.pos = T)
pha_marker$gene = paste0(pha_marker$gene,'-mRNA1')
pha_marker = pha_marker[pha_marker$p_val_adj<=0.05,]

mcltabular = read.delim('../mcl_families.OGs.tabular',header = F)

marker = pha_marker
marker$group = 0
for (i in 1:nrow(pha_marker)) {
  sub = mcltabular[mcltabular$V2 == mcltabular[grep(pha_marker[i,7],mcltabular$V1),2],]
  marker[i,'group'] = sub$V2[1]
  sub = sub[grep("Nycol",sub$V1),]
  z = pha_marker[rep(i,nrow(sub)),];z$group = sub$V2
  z$gene = sub$V1
  marker = rbind(marker,z)
}
marker$gene[grep("Nycol",marker$gene)] = substr(marker$gene[grep("Nycol",marker$gene)],1,13)
marker$gene[grep("Nycol",marker$gene)] = paste0(marker$gene[grep("Nycol",marker$gene)],'v1.2')
comarker = ny_marker[ny_marker$gene %in% marker$gene,]
pmarker = marker[marker$group %in% marker[marker$gene %in% comarker$gene,'group'],]
pmarker = pmarker[grep('PAXXG',pmarker$gene),]

p1marker = pmarker
p2marker = pmarker
p3marker = pmarker

pmarker = rbind(p1marker,p2marker,p3marker)
pmarker = rbind(p2marker,p3marker)
pmarker = pmarker[!duplicated(pmarker[,c("cluster","gene")]),]


tmpdata = GetAssayData(pha_t,layer = 'counts')
tmpdata = as.data.frame(t(apply(tmpdata, 1, scale)))
colnames(tmpdata) = colnames(pha_t)
rownames(tmpdata) = paste0(rownames(tmpdata),'-mRNA1')
spring_data = as.data.frame(list(0,0,0,0),col.names=c('x','y','cluster','sample'));spring_data = spring_data[-1,]
for (i in colnames(tmpdata)) {
  tmp_sp = as.data.frame(list(0,0),col.names=c('x','y'))
  tmp_sp$y = mean(tmpdata[pmarker[pmarker$cluster =='p',7],i])
  tmp_sp$x = mean(tmpdata[pmarker[pmarker$cluster=='s',7],i])
  rownames(tmp_sp) = i
  tmp_sp$cluster = pha_t$cluster[i]
  tmp_sp$sample =pha_t$cluster[i]
  spring_data = rbind(spring_data,tmp_sp)
}
spring_data$id = rownames(spring_data)
spring_data$cluster = as.character(spring_data$cluster)

spring_data_1 = spring_data;spring_data_1$sample = paste0(spring_data_1$sample,'5')
spring_data_2 = spring_data;spring_data_2$sample = paste0(spring_data_2$sample,'6')
spring_data_3 = spring_data;spring_data_3$sample = paste0(spring_data_3$sample,'7')

spring_data = rbind(spring_data_2,spring_data_3)
spring_data_p = rbind(spring_data_2,spring_data_3)

ny_t = subset(ny,cluster %in% c('49','51','54','60','61','63','64','65','66'))
ny_t@meta.data[ny_t$cluster %in% c('49','51','54'),'cluster'] = 'p'
ny_t@meta.data[ny_t$cluster %in% c('60','61','63','64','65','66'),'cluster'] = 's'
nymeta = read.table("/Users/wg/clean_expression/stamen2petal_metadata.txt",header = T,row.names = 1)
ny$cell = rownames(ny@meta.data)
ny_t = subset(ny,cell %in% nymeta$cell)
nymeta = nymeta[rownames(ny_t@meta.data),]
tmpdata = GetAssayData(ny_t,layer = 'counts')
tmpdata = as.data.frame(t(apply(tmpdata, 1, scale)))
colnames(tmpdata) = colnames(ny_t)
spring_data = as.data.frame(list(0,0,0,0),col.names=c('x','y','cluster','sample'));spring_data = spring_data[-1,]
for (i in colnames(tmpdata)) {
  tmp_sp = as.data.frame(list(0,0),col.names=c('x','y'))
  tmp_sp$y = mean(tmpdata[ny_marker[ny_marker$cluster=='p',7],i])
  tmp_sp$x = mean(tmpdata[ny_marker[ny_marker$cluster=='s',7],i])
  rownames(tmp_sp) = i
  tmp_sp$cluster = nymeta[i,'seurat_clusters']
  tmp_sp$sample = nymeta[i,'sample']
  spring_data = rbind(spring_data,tmp_sp)
}
tmpdata = GetAssayData(pha_t,layer = 'counts')
tmpdata = as.data.frame(t(apply(tmpdata, 1, scale)))
colnames(tmpdata) = colnames(pha_t)
rownames(tmpdata) = paste0(rownames(tmpdata),'-mRNA1')
for (i in colnames(tmpdata)) {
  tmp_sp = as.data.frame(list(0,0),col.names=c('x','y'))
  tmp_sp$y = mean(tmpdata[pmarker[pmarker$cluster =='p',7],i])
  tmp_sp$x = mean(tmpdata[pmarker[pmarker$cluster=='s',7],i])
  rownames(tmp_sp) = i
  tmp_sp$cluster = pha_t$cluster[i]
  tmp_sp$sample =pha_t$cluster[i]
  spring_data = rbind(spring_data,tmp_sp)
}
spring_data$id = rownames(spring_data)
spring_data$cluster = as.character(spring_data$cluster)
spring_data = rbind(spring_data,spring_data_p)
spring_data = rbind(spring_data_p,spring_data)
ggplot(spring_data, aes(x=x, y=y,colour = cluster)) + 
  geom_point(shape=20,size=2) + scale_color_manual(values = palette_84)+
  theme(panel.background = element_rect(fill = "white"), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.key.size = unit(10,"pt"), legend.position="right")+xlim(-3,3)+ylim(-3,3)
ggplot(spring_data, aes(x=x, y=y,colour = cluster)) + 
  geom_point(shape=20,size=2) + scale_color_manual(values = palette_84)+
  theme(panel.background = element_rect(fill = "white"),panel.border = element_blank(),
        axis.title = element_blank(),
      legend.title = element_blank(),legend.position="None")+theme_bw()+xlim(-3,3)+ylim(-3,3)

draw_axis_line = function(length_x, length_y){
  axis_x_begin = -1*length_x
  axis_x_end = length_x
  
  axis_y_begin  = -1*length_y
  axis_y_end    = length_y
  
  # set zero point
  
  data = data.frame(x = 0, y = 0)
  p = ggplot(data = data) +
    
    # draw axis line
    geom_segment(y = 0, yend = 0, 
                 x = axis_x_begin, 
                 xend = axis_x_end,
                 size = 1) + 
    geom_segment(x = 0, xend = 0, 
                 y = axis_y_begin, 
                 yend = axis_y_end,
                 size = 1) +
    # labels
    theme_void() 
  return(p)
}
p = draw_axis_line(3,3)
p = p + 
  geom_point(data=spring_data, aes(x=x, y=y,colour = cluster),shape=20,size=2) + scale_color_manual(values = palette_84[c(13,1,2,9,3:7)])+
  theme(panel.background = element_rect(fill = "white"), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.key.size = unit(10,"pt"), legend.position="right")+xlim(-3,3)+ylim(-3,3)+
  labs(x = 'time',y='distal-proximal')
p
p = draw_axis_line(1,1)
p = p + 
  geom_point(data=spring_data[spring_data$sample%in% c('petal','stamen4'),], aes(x=x, y=y,colour = sample),shape=20,size=2) + scale_color_manual(values = palette_84)+
  theme(panel.background = element_rect(fill = "white"), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.key.size = unit(10,"pt"), legend.position="right")+xlim(-1,1)+ylim(-1,1)+
  labs(x = 'time',y='distal-proximal')
p
p = draw_axis_line(4,4)
p = p + 
  geom_point(data=spring_data, aes(x=x, y=y,colour = sample),shape=20,size=2) + scale_color_manual(values = palette_84[c(13,1,2,9,3:7)])+
  theme(panel.background = element_rect(fill = "white"), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),legend.key.size = unit(10,"pt"), legend.position="right")+xlim(-1,3)+ylim(-1,3)
p

N.colorata: stamen1,stamen2,stamen3,stamen4,petal
c("wheat2", "violetred2","plum3", "deepskyblue2", "olivedrab2")

P.aphrodite: stamen-bud6, stamen-bud7, petal-bud6, petal-bud7
c("lightsalmon","tomato","plum","mediumpurple")
pale = c("plum","mediumpurple", "olivedrab2","lightsalmon","tomato","wheat2", "violetred2","plum3", "deepskyblue2")
panel.grid.major = element_line(),
panel.grid.minor = element_line(size=0.5),
p = draw_axis_line(3,3)
p = p + 
  geom_point(data=spring_data[spring_data$sample=='stamen4',], aes(x=x, y=y,colour = sample),shape=20,size=10) + scale_color_manual(values = pale[9])+
  theme(panel.background = element_rect(fill = "white"), panel.grid.major=element_line(colour = "lightgrey",size=0.5),panel.grid.minor = element_line(colour = "lightgrey",size=0.5),
        panel.border = element_rect(color = "black", size = 5, fill = NA),plot.margin = margin(t = 10, r = 10, b = 10, l = 10),legend.key.size = unit(10,"pt"), legend.position="none")+xlim(-.7,1.5)+ylim(-.7,1.5)
p
p = draw_axis_line(3,3)
p = p + 
  geom_point(data=spring_data, aes(x=x, y=y,colour = sample),shape=20,size=2.5) + scale_color_manual(values = pale)+
  theme(panel.background = element_rect(fill = "white"), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),legend.key.size = unit(10,"pt"),legend.position="none")+xlim(-.7,1.5)+ylim(-.7,1.5)
p

p = draw_axis_line(3,3)
p = p + 
  geom_point(data=spring_data, aes(x=x, y=y,colour = sample),shape=20,size=5) + scale_color_manual(values = pale)+
  theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),legend.key.size = unit(10,"pt"),legend.position="none")+xlim(-.7,1.5)+ylim(-.7,1.5)
p

p = ggplot(data=spring_data, aes(x=x, y=y,colour = sample))+
  geom_point(shape=20,size=5) + scale_color_manual(values = pale)+
  theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),legend.key.size = unit(10,"pt"),legend.position="none",
        axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+xlim(-.7,1.5)+ylim(-.7,1.5)
p

p = ggplot(data=spring_data[spring_data$sample=='stamen4',], aes(x=x, y=y,colour = sample))+
  geom_point(shape=20,size=7) + scale_color_manual(values = pale[9])+
  theme(panel.background = element_rect(fill = "white"), panel.grid.major=element_line(colour = "lightgrey",size=0.5),panel.grid.minor = element_line(colour = "lightgrey",size=0.5),
        panel.border = element_rect(color = "black", size = 3, fill = NA),plot.margin = margin(t = 10, r = 10, b = 10, l = 10),legend.key.size = unit(10,"pt"),
        legend.position="none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+xlim(-.7,1.5)+ylim(-.7,1.5)
p
saveRDS(spring_data,'../springdata.rds')

marker = pmarker[,c(7,8,6)]
for (i in 1:328) {
  sub = mcltabular[mcltabular$V2 == marker[i,2],]
  sub = sub[grep("Nycol",sub$V1),]
  z = marker[rep(i,nrow(sub)),];z$group = sub$V1
  marker = rbind(marker,z)
}
marker = marker[grep("Nycol",marker$group),]
marker = marker[!duplicated(marker),]
rownames(marker) = 1:nrow(marker)
colnames(marker) = c('P.aphrodite gene','N.colorata gene','label')
marker$label = gsub('s','stamen',marker$label)
marker$label = gsub('p','IT',marker$label)
marker$`N.colorata gene` = substr(marker$`N.colorata gene`,1,13)
marker$`N.colorata gene` = paste0(marker$`N.colorata gene`,'v1.2')
marker = marker[marker$`N.colorata gene` %in% ny_marker$gene,]
