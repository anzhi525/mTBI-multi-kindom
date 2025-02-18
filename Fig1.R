#####################Fig1 B##############################################################
control <- read_excel("mTBI_phenotypes.xlsx")

data <- read.delim("PCSS.txt")
data$Group <- factor(data$Group,levels = c("<10","10-20",">20"),ordered = TRUE)

color <- c("#8dd3c7","#bebada","#80b1d3","#fccde5","#feb24c","#ef6548","#f4a582")

p1 <- ggplot(data,aes(x="",y=n,fill=Group))+
  geom_col()+
  #geom_bar(stat="identity",width=1)+
  geom_text(aes(label=prop),position=position_stack(vjust=0.5),size=7)+
  scale_fill_manual(values = color) + 
  # scale_fill_brewer(palette="Dark2")+
  coord_polar("y",start=0)+
  theme_void()+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 22)) + 
  labs(fill="Group",x=NULL,y=NULL)

p1
ggsave(filename =paste0("PCSS.pdf"), p1, width = 100, height = 80, units="mm")
ggsave(filename =paste0("PCSS.png"), p1, width = 100, height=80, units="mm")


#####################Fig1 C##############################################################

color <- c("#8dd3c7","#bebada","#80b1d3","#fccde5","#feb24c","#ef6548","#f4a582")

p2 <- counts%>%
  ggplot(aes(x="",y=sum,fill=kingdom))+
  geom_col()+
  coord_polar(theta="y")+
  geom_label_repel(aes(label=prop,y=sum),
                   nudge_x=0.6,
                   nudge_y=0.6,
                   size=5,
                   show.legend=F,
                   segment.color="grey50"
  )+
  guides(fill=guide_legend(title="Kingdom"))+
  scale_fill_manual(values=color)+
  theme_void()+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 22))

p2
ggsave(filename =paste0("counts.pdf"), p2, width = 100, height = 80, units="mm")
ggsave(filename =paste0("counts.png"), p2, width = 100, height=80, units="mm")

#####################Fig1 D###############################

data <- read.delim("Phylum.txt")
p2 <- ggplot (data,aes (x = kingdom, y = Phylum, fill = kingdom)) + 
  geom_bar (stat="identity")+
  geom_text(aes(label = Phylum), vjust = -0.5, size = 5, color = "black", fontface = "bold") +
  scale_fill_manual(values=c("#8dd3c7","#bebada","#80b1d3","#fccde5"))+
  xlab("") + 
  ylab("") + 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(family = "sans",size = 18,colour = 'black'),
        axis.text.x = element_text(family = "sans",size = 18,colour = 'black',angle=90),
        axis.title=element_text(size=12,face="bold"),
        axis.line = element_line(colour = 'black', size = 1),
        legend.text = element_text(family = "sans",size = 20,face="italic"),
        legend.title = element_text(family = "sans",size = 20,face="bold"),
        legend.key.size=unit(5,'mm'),
        legend.position = 'right',
        panel.background = element_rect(fill = 'white')) +
  labs(fill = "Phylum")
p2

ggsave(filename =paste0("Phylum_sta.pdf"), p2, width = 130, height = 150, units="mm")
ggsave(filename =paste0("Phylum_sta.png"), p2, width=130, height=150, units="mm")


