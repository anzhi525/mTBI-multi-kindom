
level_s=c("Control","mTBI")
pl_list=

#Fig2 A-----------------------------
library(ggplot2)
library(dplyr)
library(multcompView)
library(ggpubr)
library(ggrepel)
library(ggsci)

    pl = c("Control  mTBI")
    p_list = c("Shannon", "Simpson", "InvSimpson")
	
    for (i in p_list){
      pl2 <- as.character(unlist(strsplit(pl, split = "  ")))
      metadata = read.table("metadata.txt", header=T, row.names=1, sep="\t", comment.char="",check.names=F,quot="", stringsAsFactors = F)
      metadata00 <- subset(metadata, Group %in% pl2)
      metadata00 <- metadata00[order(metadata00$Group,decreasing = TRUE),,drop=FALSE]
      otutab = read.table(paste0("alpha.species.",pl,".csv"), header=T, row.names=1, sep=",",quot="", comment.char="")
      otutab <- otutab[rownames(metadata00),]
	  groupID = "Group"
	  index = i
	  
  otutab = otutab[rownames(metadata00), ]
  sampFile = as.data.frame(metadata00[, groupID], row.names = row.names(metadata00))
  df = cbind(otutab[rownames(sampFile), index], sampFile)
  colnames(df) = c(index, "group")

  
  df$group<- factor(df$group,levels =pl2,ordered = TRUE)
  p = ggplot(df, aes(x = group, y = .data[[index]], color = group)) + 
    geom_boxplot(alpha = 1,size = 1.2, width = 0.5, fill = "transparent",outlier.shape = NA)+
	scale_color_aaas(a=0.55)+
    scale_fill_aaas(a=0.55)+
    labs(x = "Groups", y = paste(index),color = groupID) +
    theme(axis.text=element_text(size=12),
          axis.title.y=element_text(size=15,face="bold",colour = 'black'),
          axis.title.x=element_blank(),
          axis.text.x=element_text(family = "sans",colour = 'black',size=12),
          axis.line = element_line(colour = 'black', size = 0.75),
          legend.position = "top",
          text = element_text(family = "sans",size = 10,colour = 'black'),
          panel.background = element_rect(fill = 'white'),
          legend.key = element_rect(fill = 'white'))+
    geom_jitter(position = position_jitter(0.17), size = 1,alpha = 0.7,shape=16)
  p=p+
    geom_signif(comparisons = list(c(pl2)),
                map_signif_level = c("***"=0.001,"**"=0.01, "*"=0.05, " "=2),
                test =wilcox.test,
                step_increase = 0.1,
                tip_length=0,
                size=0.9, 
                col = "black",
                textsize = 7,
                vjust = 0.65,
                position = "identity",
    )+
    geom_signif(comparisons = list(c(pl2)),
                map_signif_level = F,
                test =wilcox.test,
                step_increase = 0.1,
                tip_length=0,
                size=0.9, 
                col = "black",
                textsize = 4,
                vjust = 1.65,
                position = "identity",
    )
        ggsave(filename =paste0(i," index in species level.2.pdf"), p, width =100, height = 130, units="mm")
        ggsave(filename =paste0(i," index in species level.2.png"), p, width=100, height=130, units="mm")

    }
  


#Fig2 B-----------------------------
  library(vegan)
  library(ggtext)
  library(vegan)
  library(ggsci)
  library(ggplot2)

  pl2 <- c("Control","mTBI")
  
  otutab = read.table(paste0("out_kraken2_archaea_species.txt"), header=T, row.names=1, sep="\t", comment.char="",quot="",check.names=F)
  metadata=read.table("metadata.txt", header=T, row.names=1, sep="\t", comment.char="",check.names=F, quot="",stringsAsFactors=F)

  metadata <- subset(metadata, Group %in% pl2)
  metadata <- metadata[order(metadata$Group,decreasing = TRUE),,drop=FALSE]
  otutab <- otutab[, rownames(metadata)]

  otutab <- otutab[which(rowSums(otutab) > 0),]
  Transp_otu <- t(otutab)

  bray_dis <- vegdist(Transp_otu, method = "bray")
  bray_dis1<- as.matrix(bray_dis)

  write.xlsx(bray_dis1, paste0("pcoa_bray_distance_of_species.xlsx"),append = FALSE,col.names = TRUE,row.names = TRUE)
  idx = rownames(metadata) %in% rownames(bray_dis1)
  metadata = metadata[idx, , drop = F]
  metadata1=as.data.frame(lapply(metadata,as.character),stringsAsFactors = F)
  rownames(metadata1)<-rownames(metadata)
  bray_dis1 = bray_dis1[rownames(metadata1), rownames(metadata1)]

  sampFile = as.data.frame(metadata[, Group], row.names = row.names(metadata1))
  pcoa = cmdscale(bray_dis1, k = 3, eig = T)
  points = as.data.frame(pcoa$points)
  eig = pcoa$eig
  points1 = cbind(points, sampFile)
  colnames(points1) = c("V1", "V2", "V3", Group)
  points1$Group<- factor(points1$Group,
                          levels = pl2,ordered = TRUE)

   p0 = ggplot(points1, aes(x = V1, y = V2, color=Group))
   
      p1<-p0+
      geom_point(alpha = 0.68, size = 3,shape=16)+
      scale_color_aaas(a=0.55)+
      scale_fill_aaas(a=0.55)+
      ggtitle(paste0("PCoA of species"))+
      labs(x = paste("PCo 1"," (", format(100 * eig[1]/sum(eig),digits = 4), "%)", sep = ""), 
      y = paste("PCo 0"," (", format(100 * eig[2]/sum(eig), digits = 4), "%)",sep = ""),color = Group)  +
      theme_bw()+
      theme(
        plot.title = element_text(size=20,face="bold",colour = 'black'),
        axis.text.y = element_text(family = "sans",size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,family = "sans",colour = 'black'),
        axis.title=element_text(size=12,face="bold",colour = 'black'),
        legend.text = element_text(family = "sans",size = 9,colour = 'black'),
        legend.position = 'right')+
      stat_ellipse(level = 0.68)

    p1 = p1 + geom_text_repel(label = paste(rownames(points1)),colour = "black", size = 3.5)

  otu.adonis <- adonis2(t(otutab) ~ Group, data = metadata, permutations = 9999,method="bray")

  p4 <- p1 + ggtext::geom_richtext(
    aes(vjust = 1.1,hjust = 1.05,
        x = Inf , y = Inf,
        label = paste('PERMANOVA: <p>\n  p-value = ',
                      otu.adonis$`Pr(>F)`[1], sep = "")),
    size =3, family = "sans",colour = 'black')

      ggsave(filename = paste0("PCoA of species.pdf"), p4, width=150, height=120, units="mm")
      ggsave(filename = paste0("PCoA of species.png"), p4, width=150, height=120, units="mm")
  


#Fig2 C-----------------------------  
    library(ggplot2)
    library(reshape2)
    library(tidyverse)
    library(plotly)
    library(RColorBrewer)

	pl2 <- c("Control","mTBI")
    metadata = read.table("metadata.txt", header=T, row.names=1, sep="\t", comment.char="",check.names=F, quot="",stringsAsFactors = F)
    metadata <- subset(metadata, Group %in% pl2)
    otutab = read.table(paste0("out_kraken2_archaea_species.txt"), header=T, row.names=1, sep="\t", comment.char="",quot="",check.names=F)
    otutab <- otutab[,rownames(metadata)]
	topN=20
    groupID="Group"

    idx = rownames(metadata) %in% colnames(otutab)
    metadata = metadata[idx, , drop = F]
    otutab = otutab[, rownames(metadata)]
    otutab <- otutab[which(rowSums(otutab) > 0),]
    sampFile = as.data.frame(metadata[, groupID], row.names = row.names(metadata))
    colnames(sampFile)[1] = "group"
    mean_sort = as.data.frame(otutab[(order(-rowSums(otutab))), ])
    idx = grepl("unassigned|unclassified|unknown", rownames(mean_sort),ignore.case = T)
    mean_sort = rbind(mean_sort[!idx, ], mean_sort[idx, ])
    other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
    bartxt = mean_sort[1:topN,]
    mean_sort = mean_sort[1:(topN - 1), ]
    mean_sort = rbind(mean_sort, other)
    rownames(mean_sort)[topN] = c("Other")
	
      mean_sort$Taxonomy = rownames(mean_sort)
      data_all_sample = as.data.frame(melt(mean_sort, id.vars = c("Taxonomy")))
      data_all_sample$Taxonomy = factor(data_all_sample$Taxonomy, levels = rownames(mean_sort))
      data_all_sample = merge(data_all_sample, sampFile, by.x = "variable", 
                              by.y = "row.names")
      data_all_sample$group<- factor(data_all_sample$group,levels =c("Control","mTBI"),ordered = TRUE)
      p = ggplot(data_all_sample, aes(x = variable, y = value, fill = Taxonomy)) + 
        geom_bar(stat = "identity", width = 0.8,size =0.05) + 
        facet_grid(~group, scales = "free_x") +
        xlab("") + 
        ylab("Relative abundance(%)") + 
        theme(panel.grid = element_blank(),
              axis.text.y = element_text(family = "sans",size = 10,colour = 'black'),
              axis.text.x = element_text(family = "sans",size = 7,colour = 'black',angle=90),
              axis.title=element_text(size=12,face="bold"),
              axis.line = element_line(colour = 'black', size = 0.75),
              legend.text = element_text(family = "sans",size = 9,face="italic"),
              legend.title = element_text(family = "sans",size = 15,face="bold"),
              legend.key.size=unit(5,'mm'),
              legend.position = 'right',
              panel.background = element_rect(fill = 'white'))+
        guides(fill=guide_legend(ncol=1))+
        labs(fill = paste("Top", topN, "of species", sep = " "))
        getPalette = colorRampPalette(brewer.pal(12, "Set3"))
        p=p+scale_fill_manual(values = getPalette(topN))        
    p
	ggsave(path = newdir7, filename =paste0("Top.abundant.species.sample.1.pdf"), p, width=300, height=150, units="mm")
    ggsave(path = newdir7, filename =paste0("Top.abundant.species.sample.1.png"), p, width=300, height=150, units="mm")


#Fig2 F-----------------------------  
  library(ggrepel)
  library(rio)
  library(openxlsx)
  library(ropls)

    pl2 <- c("Control","mTBI")
    metadata_all <- read.table( "metadata.txt", header=T, sep="\t", comment.char="",check.names=F, quot="",stringsAsFactors = F)
    metadata_00 <- subset(metadata_all, Group %in% pl2)
    rownames(metadata_00) <- metadata_00$ID
    metadata_00 |>dplyr::select(-ID) -> metadata_00

    data_00 <- read.table(paste0("out_kraken2_archaea_species.txt"), header=T, row.names=1, sep="\t",quot="", comment.char="",check.names=F)
    data_11 <- data_00[,colnames(data_00)%in%row.names(metadata_00)]
    otu_table <- t(data_11)

  experiment <- pl2[2]
  control <- pl2[1]
  print(pl2)
  options(scipen = 200)
  otu_table <- as.data.frame(otu_table[which(rowSums(otu_table) > 0),])
  oplsda = opls(otu_table, metadata_00$Group, predI = 1, orthoI = 1)#, orthoI = NA)
  vip <- getVipVn(oplsda)
  vip_data<-as.data.frame(vip)
  
  nr <- ncol(otu_table)                   
  cn <- colnames(otu_table)               
  nl <- ncol(otu_table)                    

  otu_table2 <-as.data.frame(otu_table)
  metadata_00 <-as.data.frame(metadata_00)
  rownames(otu_table2)
  len1 <- length(pl2)*3+6
  len1
  out <- matrix(NA, nr, len1) 
  all10data<- NA
  
  for (i in 1:nr) {                 
    tem_otu <- merge(as.data.frame(otu_table2[i]),metadata_00,by="row.names")
    colnames(tem_otu)[2]<-"speciesid"
    tem_otu$speciesid =as.numeric(tem_otu$speciesid)
    tem_otu$Group = factor(tem_otu$Group,
                           levels =pl2,ordered = TRUE)
    tem_otu = as.data.frame(tem_otu[(order(tem_otu$Group)), ])
    if(length(pl2)==2){
      wilcox_pvalue <- wilcox.test(speciesid~Group,data = tem_otu)$p.value
      wilcox_statistic <- coin::wilcox_test(speciesid~Group,data = tem_otu,distribution = "exact")@statistic@teststatistic
    }
    mean_abun <- tapply(tem_otu$speciesid, tem_otu$Group, mean)
    fold_change <- mean_abun[experiment]/mean_abun[control]
    log_foldchange <- log(fold_change ,2)
    if(length(pl2)==2){
      effect.z <- coin::wilcox_test(speciesid~Group,tem_otu)@statistic@standardizedlinearstatistic
      effect_size <- effect.z/sqrt(length(rownames(tem_otu)))
    }

    len0 <- length(pl2)*2-1
    conclude <- rep(0, len0)
    seq(from=1,by=2,to=len0) 
    conclude[seq(from=1,by=2,to=len0)] <- rownames(mean_abun)[order(mean_abun)]
    if (wilcox_pvalue == 1 | is.na(wilcox_pvalue)) {
      conclude[seq(from=2,by=2,to=len0)] <- "<=>"
    }else if(wilcox_pvalue >= 0 & wilcox_pvalue < 0.05 ){
      conclude[seq(from=2,by=2,to=len0)] <- "<"
    }else{
      conclude[seq(from=2,by=2,to=len0)] <- "<="
    }
    compare <- paste(conclude, collapse = "")

    res <- c(wilcox_pvalue,wilcox_statistic,  mean_abun,fold_change, log_foldchange,effect_size,compare)#t_pvalue,t_statistic,

    all10data <- rbind(all10data,res)
  }
  out<-all10data
  out1 <- data.frame(out)
  out1 = out1[-1,] 
  p.wilcox <- out1[,1]
  p.adj.wilcox <- p.adjust(p.wilcox, "BH")
  p.adj.wilcox <- as.matrix(p.adj.wilcox)
  p.t <- out1[,1]
  p.adj.t <- p.adjust(p.t, "BH")
  p.adj.t <- as.matrix(p.adj.t)

  out1 <- cbind(out1, p.adj.wilcox)
  out1 <- cbind(out1, p.adj.t)  
  out1 <- out1[1:i,]
  row.names(otu_table)
  rownames(out1) <- row.names(t(otu_table))
  vip_data
  out2 <- merge(out1, vip_data,by="row.names")
  rownames(out2) <- out2$Row.names
  out2 <- out2[,-1] 

  cn <- c("wilcox_Pvalue","wilcox_statistic",  paste("mean_abun", pl2, sep = "_"),"FC", "log2_FC","effect_size","enrichment","Qwalue_wilcox","Qwalue_t","VIP")
  colnames(out2) <- cn  
  
  cut_off_pvalue = 0.05
  cut_off_logFC = 1    
  out2$change = ifelse(out2$wilcox_Pvalue <= cut_off_pvalue & abs(as.numeric(out2$log2_FC)) >= cut_off_logFC, 
                       ifelse(out2$log2_FC> cut_off_logFC ,'Up','Down'),
                       'Stable')
  out3 <- merge(out2, t(otu_table),by="row.names") 
  rownames(out3) <- out3$Row.names  

  out3 = as.data.frame(out3[order(out3$wilcox_Pvalue), ])

  newdir26=paste0("species analyze/volcano/",pl)
  write.xlsx(out3,file = paste0(newdir26,"/volcano of species.xlsx"))
  
  fit <- try(subset(out3, out3$change == "Up" | out3$change == "Down" ))
  if("try-error" %in% class(fit))
  {
    next
  }else
  {
    out4 <- subset(out3, out3$change == "Up" | out3$change == "Down" )
  }
  
  if (!is.null(volcano_picture)){
    otu_table1 <- out3[,c("wilcox_Pvalue","log2_FC","VIP","change")]
    otu_table1$name<-rownames(otu_table1)
    cut_off_pvalue = 0.05
    cut_off_logFC = 0.5849625    
    cut_off_VIP = 1    

    otu_table1$label = ifelse(otu_table1$wilcox_Pvalue <= cut_off_pvalue & abs(as.numeric(otu_table1$log2_FC)) >= cut_off_logFC & otu_table1$VIP >cut_off_VIP, 
                              otu_table1$name,"")
    otu_table1$log2_FC <-as.numeric(otu_table1$log2_FC)
    otu_table1$wilcox_Pvalue <-as.numeric(otu_table1$wilcox_Pvalue)

    p=ggplot(
      otu_table1, 
      aes(x = log2_FC, 
          y = -log10(wilcox_Pvalue), 
          colour=change)) +
      geom_point(alpha=0.9, size=3.5) +
      scale_color_manual(values=c("Down" = "#3B4992FF", "Stable"="#d2dae2","Up"="#EE0000FF"))+
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) 
      labs(x="log2(Fold change)",
           y="-log10 (P-value)")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5), 
            axis.text = element_text(color = 'black',size = 10, family = 'sans', face = 'plain'),
            axis.title = element_text(color = 'black',size = 15, family = 'sans', face = 'plain'),
            
            legend.position="right", 
            
            legend.title = element_blank()
      )+
      geom_text_repel(data=otu_table1,aes(label = label), 
                      colour = "black", size = 3)
    p
    ggsave(filename =paste0("volcano of species.pdf"), p, width = 200, height = 150, units="mm")
    ggsave(filename =paste0("volcano of species.png"), p, width = 200, height = 150, units="mm")
    }





#Fig2 G-----------------------------  
library(ppcor)
library("stringi")
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

 pheno = read.table("diff_pheno_new.txt",header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,check.names = FALSE,quot="")
 micro <-   read.table(paste0("diff_species.txt"),header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,check.names = FALSE,quot="")
 micro<- micro[,colnames(micro)%in%rownames(pheno1)]
 micro1 <- as.data.frame(t(micro))
 all10data = as.data.frame(matrix(nrow=0,ncol=8))
 mic_list = colnames(micro1)
 for (mic in mic_list){
   met_list=colnames(pheno1)
   for (met in met_list){           
     spear_out <- pcor.test(micro1[,mic],pheno1[,met],method = "spearman")
     spear_out1 <- cbind(mic,met,spear_out)
     all10data <- rbind(all10data,spear_out1)
   }
 }

datamat <- all10data
datamat1 <- datamat[,c(1:3)]
datamat2 <- datamat[,c(1,2,4)]
plotmat <- datamat1 %>% pivot_wider(names_from = met, values_from = estimate)
plotmat <- as.data.frame(plotmat)
rownames(plotmat) <- plotmat$mic
stars <- datamat2 %>% pivot_wider(names_from = met, values_from = p.value)
stars <- as.data.frame(stars)
rownames(stars) <- stars$mic
plotmat <- as.data.frame(t(plotmat))
stars <- as.data.frame(t(stars))

if (!is.null(stars)){
  ssmt <- stars< 0.01
  stars[ssmt] <-'**'
  smt <- stars >0.01& stars <0.05
  stars[smt] <- '*'
  stars[!ssmt&!smt]<- ''
} else {
  stars <- F
}

maxr <- round(max(plotmat),4)*100
minr <- abs(round(min(plotmat),4)*100)
mycol <- c(colorRampPalette(c('#31688E','white'))(minr),
           colorRampPalette(c('white','white'))(3),
           colorRampPalette(c('white','#FB9A99'))(maxr))
annotation2 = as.data.frame(read.table(paste0("pheatmap.species.ann.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))
annotation2 <- annotation2[colnames(plotmat),]

maxr <- round(max(plotmat),4)
minr <- round(min(plotmat),4)

Enrichment_col<- c("#1B9E77", "#D95F02")
names(Enrichment_col) = c("mTBI" , "Control" )
getPalette_kindom = colorRampPalette(brewer.pal(3, "Set1"))
kindom_col <- getPalette_kindom(3)
names(kindom_col) <- unique(annotation2$Kingdom)
getPalette_phylum = colorRampPalette(brewer.pal(6, "Set3"))
Phylum_col <- getPalette_phylum(6)
names(Phylum_col) <- unique(annotation2$Phylum)

getPalette_family = colorRampPalette(brewer.pal(9, "Set2"))
Family_col <- getPalette_family(9)
names(Family_col) <- unique(annotation2$Family)

species_ann <- HeatmapAnnotation(
  Kingdom= annotation2$Kingdom,
  Phylum = annotation2$Phylum,
  Family = annotation2$Family,
  Enrichment =annotation2$Enrichment,
  col = list( 
    Kingdom=kindom_col,
    Phylum = Phylum_col , 
    Family = Family_col,
    Enrichment = Enrichment_col
  )
)
cell_fun <-function(j, i, x, y, width, height, fill) {
  grid.text(stars[i,j], 
            x = x, 
            y =  y-height*0.3,
            gp = gpar(fontsize = 15,col="black"))
}
p=Heatmap(as.matrix(plotmat),cluster_columns = T,cluster_rows = F,
          width = unit(12, "cm"),
          height = unit(2, "cm"),
          row_names_gp = gpar(fontsize = 11,col ="black"),
          column_names_gp = gpar(fontsize = 11,col ="black"),
          row_names_max_width = unit(110,"mm"),
          heatmap_legend_param = list(
            at = c(-0.4,-0.2, 0, 0.2,0.4)
          ),
          cell_fun = cell_fun,
          rect_gp = gpar(col = "grey", lwd = 1),
          name = "Correlation",
           top_annotation = species_ann,
          col = mycol)
p
pdf("heatmap.2.pdf",  width=10, height=6)
p
dev.off()

png("heatmap.2.png",  width=700, height=500,units = "px")
p
dev.off()
