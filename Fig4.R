#Fig4 A-----------------------------
library(qgraph)
library(RColorBrewer)
library(tidyverse)

pl2 <- c("Control","mTBI")

metadata = read.table("metadata.txt", header=T, row.names=1, sep="\t", comment.char="",check.names=F, quot="",stringsAsFactors = F)
metadata00 <- subset(metadata, Group %in% pl2)
metadata00$SampleID <- rownames(metadata00)

outdata = read.table(paste0("diff_species.txt"), header=T, row.names=1, sep="\t", quot="",comment.char="",check.names=F)
outdata00 <- outdata[,rownames(metadata00)]
outdata00 <- outdata00[which(rowSums(outdata00)>0.0005),]

otu_data_all= read.table(paste0("taxnomy.txt"), header=F, sep="\t", comment.char="#", quot="",stringsAsFactors = F)
otu_data_all1 <- otu_data_all
character_to_count <- "|"

filtered_rows <- otu_data_all1[str_count(otu_data_all1$V1, fixed(character_to_count)) == 7,]
filtered_rows <- as.data.frame(filtered_rows)
otu_data_all2 <- separate(data=filtered_rows, col=filtered_rows, into=c("domain","Kindom","Phylum","Class","Order","Family","Genus","Species"), sep = "\\|", remove = FALSE, convert = FALSE)
otu_data_all2$Species_name <-gsub("s__", "", otu_data_all2$Species)

otu_data_all2$Kindom[otu_data_all2$Kindom == ""] <- "k__other"
otu_data_all2$Phylum[otu_data_all2$Phylum == ""] <- "p__other"
otu_data_all2$Class[otu_data_all2$Class == ""] <- "c__other"
otu_data_all2$Order[otu_data_all2$Order == ""] <- "o__other"
otu_data_all2$Family[otu_data_all2$Family == ""] <- "f__other"
otu_data_all2$Genus[otu_data_all2$Genus == ""] <- "g__other"
otu_data_all2$Species[otu_data_all2$Species == ""] <- "s__other"

taxonomy_table <- otu_data_all2 |>  dplyr::select(Species_name,domain,Kindom,Phylum,Order,Class,Family,Genus,Species) #|>

taxonomy_table <- taxonomy_table %>%
  mutate(phylum  = gsub("^p__", "", Phylum),
         class   = gsub("^c__", "", Class),
         order   = gsub("^o__", "", Order),
         family  = gsub("^f__", "", Family),
         genus   = gsub("^g__", "", Genus),             
         species = gsub("^s__", "", Species),                
  )

taxonomy_table$Kingdom <- taxonomy_table$domain
taxonomy_table$Kingdom <-gsub("d__", "",taxonomy_table$Kingdom)
taxonomy_table$Kingdom <-gsub("Eukaryota", "Fungi",taxonomy_table$Kingdom)
taxonomy_table <- taxonomy_table |>  dplyr::select("Species",Kingdom) |> distinct() 

colnames(taxonomy_table) <- c("Species_name","Kingdom")
taxonomy_table$Species_name <- make.names(taxonomy_table$Species_name)

top <- read.table("feature.csv",header=T)
top.tax2 <- top$feature

annotation <- taxonomy_table
colnames(annotation) <- c("node","Phylum")

metadata11 = metadata00
outdata = outdata00

otu11 <- outdata[,rownames(metadata)]
metadata <- subset(metadata11, Group %in% pl2[i])

otu <- otu11[,colnames(otu11)%in%rownames(metadata)]

A=t(otu)
C=A/rowSums(A)
otu2=t(C)
otu3 <- as.data.frame(otu2[which(rowSums(otu2) >= 0), ])
otu4 <- otu[rownames(otu3),]

CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat) 
  data.frame(
    from = rownames(cormat)[col(cormat)[ut]],
    to = rownames(cormat)[row(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}
otu0 <- otu4
occor <- corAndPvalue(t(otu0), use='pairwise', method="spearman") 
cor_df <- CorrDF(occor$cor , occor$p) 
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.5),] 
cor_df <- cor_df[which(cor_df$p < 0.05),] 
suppressWarnings(write.table(cor_df, file=paste0(pl2[i],".cor.0.5_p0.05.txt"), append = F, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T))

 
e2.case <- cor_df[,1:3]
colnames(e2.case) <- c("genus","variable","value")
e2.case = e2.case[which(e2.case$value != 0), ]
edge.case = e2.case[which(e2.case$genus != e2.case$variable), ]
edge.case2 = edge.case
node.case = data.frame(unique(c(
  as.character(edge.case2$genus),
  as.character(edge.case2$variable)
)))
colnames(node.case) = "node"
rownames(node.case) = node.case$node
node.case$score <- 1
node.case2 = node.case
node.case2$weight = abs(node.case2$score)
node.case2$class = 1 
node.case2[top.tax2, "class"] = 2 
node.case2 <- node.case2[complete.cases(node.case2), ]
edge.case2$weight = 0.1
for (kk in 1:length(top.tax2))
{
  for (mm in 1:nrow(edge.case2))
  {
    if (as.character(edge.case2[mm,1]) == top.tax2[kk] || as.character(edge.case2[mm,2]) == top.tax2[kk])
    {
      edge.case2[mm,"weight"] = 1
    }else
    {
      next
    }
  }
}
edge.case2$class = 1
edge.case2[which(edge.case2$value < 0), "class"] = 2  
node.case2 <-   left_join(node.case2,annotation,by = "node") 
g1 <- graph.empty()
g1 <- graph_from_data_frame(edge.case2, vertices = node.case2)
nodeSize <- 1.5
nodeDize <- 1.2
edgeSize <- 0.3
edgeDize <- 0
arrowSize = 0
my.layout = layout.sphere

EColor <- c("#ea66a6", "#2a5caa")

getPalette = colorRampPalette(brewer.pal(4, "Set1"))
my_color <- getPalette(length(levels(as.factor(V(g1)$Phylum))))
VText <- c(0.2,1)
V(g1)$size <- nodeSize + nodeDize * 10 * as.numeric(as.vector(node.case2$weight*0.5))

V(g1)$color <- my_color[factor(node.case2$Phylum)] 
V(g1)$label.cex <- VText[node.case2$class]
V(g1)$frame.color <- "black"
E(g1)$width <-
  edgeSize + (edgeDize * abs(3 * as.numeric(as.vector(
    edge.case2$weight
  ))))
E(g1)$color <- EColor[edge.case2$class]
E(g1)$arrow.size <- arrowSize

g1_2 <- g1
V(g1_2)$label <- ifelse(V(g1_2)$name %in% top.tax2, as.character(V(g1_2)$name), "")

e <- get.edgelist(g1_2,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g1_2),
                                       area=8*(vcount(g1_2)^2.3),repulse.rad=(vcount(g1_2)^3.1))
plot(
  g1_2,
  layout = l,
  edge.width=1.5, 
  vertex.label.color=c("black"),
  vertex.frame.color = "black",
  vertex.shapes = "none"
)
legend("right", legend =levels(as.factor(V(g1)$Phylum)),
       col = my_color,
       horiz = FALSE,pch = 16,text.width = 1 / 50,cex = 2,
       inset = 0.05, xpd =  FALSE,bty = "n",title = "Kingdom")

p = myplot({
  plot(g1_2,layout=layout.sphere,
       vertex.color=my_color,
       vertex.label.color=c("black"),
       vertex.frame.color = "black",
       vertex.shapes = "none"
  )
  legend("right", legend =levels(as.factor(V(g1)$Phylum)),
         col = my_color,
         horiz = FALSE,pch = 16,text.width = 1 / 50,cex = 1.3,
         inset = 0.1, xpd =  FALSE,bty = "n",title = "Kingdom")
})
p
plotsave(p,file = paste0("Co-occurrence.pdf"),width =350, height = 250, units="mm",limitsize=FALSE)
plotsave(p,file = paste0("Co-occurrence.png"),width = 350, height = 250, units="mm",limitsize=FALSE)
