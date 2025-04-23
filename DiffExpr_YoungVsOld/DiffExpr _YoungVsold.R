rm(list = ls())
setwd("/data3/lishuhan/30GSE_health_SRR/v3AS/")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
###1.Transcription analysis-----------------------------------------------------------------

a<-read.csv("./DEG/df_0_17.csv",row.names = 1)
b<-read.csv("./DEG/df_60_100.csv",row.names = 1)
colnames(b)[1]<-"Sample_Name"
c<-rbind(a,b)
write.csv(c,"./DEG/metadata.csv")

c<-read.csv("./DEG/metadata.csv",row.names = 1)
table(c$dataset)[order(names(table(c$dataset)))]
meta1<-c[order(c$dataset),]
df1<-y_old_count11[,rownames(meta1)]
meta1 <- meta1 %>%
  mutate(type = ifelse(Age >= 60, "old", "young"))
meta1<-meta1[,c(5,6)]

###Remove batch--------------------------------------------------------------------------
design <- model.matrix(~0+meta1$type+meta1$sex,data = df1)
rB_Expr<-log2(df1+1)

rB_Expr <- removeBatchEffect(rB_Expr,batch = factor(meta1$dataset,meta1$sex),design = design)
write.table(rB_Expr,
            file = "./DEG/rB_Expr.tab",
            quote = FALSE,sep="\t",row.names = TRUE)
###DEGs---------------------------------------------------------------------------------------
boxplot(rB_Expr,outline=FALSE,notch=T,col=factor(meta1$type,meta1$sex),las=2)
exprSet<-normalizeBetweenArrays(rB_Expr)
boxplot(exprSet,outline=FALSE,notch=T,col=factor(meta1$type,meta1$sex),las=2)
range(exprSet)
#exprSet<-log2(exprSet+1)
colnames(design)<-c("old","young")
contr.matrix <- makeContrasts(
  oldvsyoung=old-young,
  levels = colnames(design))
contr.matrix

lfit <- lmFit(exprSet, design)
vfit <- contrasts.fit(lfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
TRD_fit<-topTable(efit,sort.by = "P",adjust="BH" , coef = 'oldvsyoung',number=Inf)

TRD_pvalue<-subset(TRD_fit, adj.P.Val < 0.01 )
TRD_gene <-subset(TRD_fit, adj.P.Val < 0.01 & abs(logFC) >= log2(4))
TRD_down <- subset(TRD_fit, adj.P.Val < 0.01 & logFC <= -log2(4))
TRD_up <- subset(TRD_fit, adj.P.Val < 0.01 & logFC >= log2(4))
TRD_fit$threshold <- as.factor(ifelse(TRD_fit$P.Value< 0.05 & abs(TRD_fit$logFC) >= log2(2),
                                      ifelse(TRD_fit$logFC > log2(2) ,'Up','Down'),'Not'))
a<-read.csv("./DEG/y_old_DEG.csv")
write.csv(TRD_fit,"./DEG/y_old_DEG.csv")

##plot
df<-read.csv("./DEG/y_old_DEG.csv",row.names = 1)
df$gene <- rownames(TRD_fit)

gene <- read.delim("./DEG/top20gene.txt"
                   ,sep="\t",header = T,quote = "",fill=TRUE)$gene
df$label=""
df$label[match(gene,df$gene)] <- gene 
#write.table(df,"result.txt",sep='\t',row.names=F,quote=F) 

df$color <- ifelse(df$threshold == "Not" & df$label == "", "#bcbcbc",  
                   ifelse(df$threshold == "Up" & df$label == "", "#ffab84",   
                          ifelse(df$threshold == "Down" & df$label == "", "#8abddc",                                   ifelse(df$threshold == "Up" & df$label != "", "#be0001", "#0051a6"))))  

df <- df %>% arrange(color)
df$color <- factor(df$color, levels = c("#bcbcbc", "#ffab84", "#8abddc", "#be0001", "#0051a6"))

df$logP<- -log10(df$P.Value)
p <- ggscatter(df,
               x="logFC",  
               y="logP",  
               color = "color",  
               palette = c("#bcbcbc","#ffab84","#8abddc","#be000e","#0051a6"), 
               label = df$label,   
               font.label = c(15,"plain","black"),  
               repel = T ) +   
  labs( 
    x=expression(paste(Log[2], 'Fold Change'),color="black", size=12,face = "bold"),      y=expression(paste(-Log[10], 'P-value'),color="black", size=12,face = "bold"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="#bcbcbc",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="#bcbcbc",lwd=0.6)+
  theme(legend.position="none",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(color="black", size=18,face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"),
        axis.title.x = element_text(color="black", size=12,face = "bold"),
        axis.title.y = element_text(color="black", size=12,face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2))+
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50))+
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 2))

print(p)
ggsave(plot=p,"./DEG/plot/v2_y_old_DEG_vol.pdf",
       width = 10, height = 8, dpi = 600, device='pdf')

##heatmap
DEG<-read.csv("./DEG/y_old_DEG.csv",row.names = 1)
DEG_up<-DEG[DEG$threshold=="Up",]
DEG_down<-DEG[DEG$threshold=="Down",]
up<-DEG_up[order(DEG_up$P.Value,decreasing = FALSE),][1:10,]
down<-DEG_down[order(DEG_down$P.Value,decreasing = FALSE),][1:10,]
all<-rbind(DEG_up,DEG_down)

data <- read.delim("./DEG/rB_Expr.tab", header = TRUE, sep = "\t")


colnames(meta1)[2]<-"Group"
colnames(meta1)[1]<-"Dataset"

ann_color = list(group = c(young="#7293cb", old="#ed7c68"))
ann_color <- list(
  type = c("young"="#7293cb", "old"="#ed7c68"),
  dataset = c("GSE102114" = "#f6a889", "GSE123658" = "#fddb83","GSE124326" = "#db8b8f",
              "GSE134080" = "#fad2b6","GSE169687" = "#f2afa4","GSE181228" = "#fcd39b",
              "GSE182038" = "#f4d7c1","GSE191238" = "#f5b6b5","GSE209591" = "#f9c7a7",
              "GSE222889" = "#b4d0f7","GSE79362" = "#81a4eb","GSE94438" = "#afced7")
)


meta2<-meta1[order(meta1$Group),]
data2<-data[all$gene,rownames(meta2)]
pdf("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/DEG/plot/v3_DEG_heatmap.pdf",width = 10,height = 9)
phea<-pheatmap(as.matrix(data2),
               annotation_col=meta2,
               border_color = NA,
               color = colorRampPalette(c("blue","white","red"))(30),
               cluster_cols = F,
               cluster_rows = F,
               show_rownames = F,
               show_colnames = F,#c("c1","c2","c3","z1","z2","z3"),  
               scale = "row", ## none, row, column
               fontsize = 12,
               fontsize_row = 12,
               fontsize_col = 10,
               annotation_legend=TRUE,
               annotation_colors = ann_color,
               border = FALSE,
               use_raster=F)
phea
dev.off()

genes <- data.frame(genes=c("E2F1","IGF1","IGF2","RET","IGFBP3",
                            "TFAP2A","HOXB7","PAPPA","TP73" ,
                            "C1QA","MT1E","KL","FLT1","HESX1","LEP","PPARGC1A","IGFBP2"))
location<-which(rownames(data2) %in% genes$genes)
gene<-genes$genes[genes$genes %in% rownames(data2)]

phea+rowAnnotation(link = anno_mark(at = location, 
                                    labels = gene, labels_gp = gpar(fontsize = 10)))

pdf("./DEG/plot/DEG_heatmap.pdf",width = 10,height = 9)

p1<-phea+rowAnnotation(link = anno_mark(at = location, 
                                        labels = gene, labels_gp = gpar(fontsize = 10)))
p1
dev.off()
ggsave("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/DEG/plot/DEG_heatmap.pdf",p1,width = 10, 
       height = 9, dpi = 300, device='pdf')
# phea+rowAnnotation(link = anno_mark(at = which(rownames(data[all$gene,]) %in% genes$genes), 
#                                     labels = genes$genes, labels_gp = gpar(fontsize = 10)))




ggsave("./DEG/plot/DEG_heatmap.pdf",phea,width = 10, 
       height = 9, dpi = 300, device='pdf')
ggsave("./DEG/plot/DEG_heatmap.pheat.png",phea,width = 10, height = 9, dpi = 300, device='png')

###enrichment analysis-------------------------------------------------------------------------
rm(list=ls())
DEG<-read.csv("./DEG/y_old_DEG.csv",row.names = 1)
DEG1<-subset(DEG,DEG$threshold %in% c("Up","Down"))
##GO
load(file = "./Genome/enrichment_annotation/human_enrichment.RData")
term2gene <- gobp[,c(1,5)]
term2name <- gobp[,c(1,2)]


DEG_go_bp <- enricher(gene = DEG1$gene,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                      TERM2GENE = term2gene,TERM2NAME = term2name)
DEG_go_bp<-data.frame(DEG_go_bp@result)

term2gene <- gomf[,c(1,5)]
term2name <- gomf[,c(1,2)]

DEG_go_mf <- enricher(gene = DEG1$gene,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                      TERM2GENE = term2gene,TERM2NAME = term2name)
DEG_go_mf<-data.frame(DEG_go_mf@result)


library(viridis)
DEG_go_bp$ONTOLOGY<-"BP"
DEG_go_mf$ONTOLOGY<-"MF"

allgg<-rbind(DEG_go_bp,DEG_go_mf)

df_top5 <- allgg %>%
  group_by(ONTOLOGY) %>%
  arrange(pvalue) %>%
  slice_head(n = 10)

df_top5$ONTOLOGY <- factor(df_top5$ONTOLOGY, levels=c('BP', 'MF'))
df_top5$Description <- factor(df_top5$Description, levels = rev(df_top5$Description))

options(repr.plot.width=7, repr.plot.height=6.5)

mycol3 <- c('#6BA5CE', '#F5AA5F')
cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")

p <- ggplot(data = df_top5, aes(x = Count, y = Description, fill=ONTOLOGY)) +
  geom_bar(width = 0.5,stat = 'identity') +
  theme_classic() + 
  scale_x_continuous(expand = c(0,0.5)) +
  scale_fill_manual(values = alpha(mycol3, 0.66))

p <- p + theme(axis.text.y = element_blank()) + 
  geom_text(data = df_top5,
            aes(x = 0.1, y = Description, label = Description),
            size = 4.8,
            hjust = 0) 




p <- p + labs(title = 'Enriched top 5 BP and MF') +
  theme(
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.ticks.y = element_blank())


ggsave("./DEG/plot/GOtop10.pdf", p, w=7, h=6.5)

##KEGG
load(file = "./Genome/enrichment_annotation/human_enrichment.RData")
term2gene <- kegg[,c(1,5)]
term2name <- kegg[,c(1,2)]
DEG_kegg <- enricher(gene = DEG1$gene,
                     pvalueCutoff = 0.05,pAdjustMethod = "BH",
                     minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                     TERM2GENE = term2gene,TERM2NAME = term2name)
DEG_kegg<-data.frame(DEG_kegg@result)

DEG_kegg$ONTOLOGY<-"KEGG"



df_top5 <- DEG_kegg %>%
  group_by(ONTOLOGY) %>%
  arrange(pvalue) %>%
  slice_head(n = 20)

df_top5$ONTOLOGY <- factor(df_top5$ONTOLOGY, levels=c('BP', 'MF'))
df_top5$Description <- factor(df_top5$Description, levels = rev(df_top5$Description))

options(repr.plot.width=7, repr.plot.height=6.5)

mycol3 <- c('#6BA5CE')
cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")

p <- ggplot(data = df_top5, aes(x = Count, y = Description, fill=ONTOLOGY)) +
  geom_bar(width = 0.5,stat = 'identity') +
  theme_classic() + 
  scale_x_continuous(expand = c(0,0.5)) +
  scale_fill_manual(values = alpha(mycol3, 0.66))

p <- p + theme(axis.text.y = element_blank()) + 
  geom_text(data = df_top5,
            aes(x = 0.1, y = Description, label = Description),
            size = 4.8,
            hjust = 0) 




p <- p + labs(title = 'Enriched top 20 KEGG') +
  theme(
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.ticks.y = element_blank())


ggsave("./DEG/plot/KRGGtop20.pdf", p, w=7, h=6.5)

###GSEA------------------------------------------------------------------------------------------------
KEGG_kk_entrez <- gseKEGG(geneList     = geneList,
                          organism     = "hsa", 
                          pvalueCutoff = 0.25)   

KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                             OrgDb="org.Hs.eg.db",
                             keyType='ENTREZID')             

GO_kk_entrez <- gseGO(geneList     = geneList,
                      ont          = "ALL",  
                      OrgDb        = "org.Hs.eg.db",
                      keyType      = "ENTREZID",
                      pvalueCutoff = 0.25)   

GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb=org.Hs.eg.db,
                           keyType='ENTREZID')
save(KEGG_kk_entrez, GO_kk_entrez, file = "./GSEA/GSEA_result.RData")
load("./GSEA/GSEA_result.RData")

##KEGG
kk_gse <- GO_kk
kk_gse_entrez <- GO_kk_entrez


kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.25 & abs(kk_gse$NES)>1]
write.csv(kk_gse_cut,"./GSEA_KEGG.csv")
write.csv(kk_gse_cut,"./GSEA_GO.csv")
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]


down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),10),]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$NES),decreasing = T),10),]

up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T)),]

library(gseaplot2)
up_gsea$Description
i=2
gseap1 <- enrichplot::gseaplot2(kk_gse,
                                up_gsea$ID[5],
                                title = up_gsea$Description[5],
                                color = "red",
                                base_size = 20,
                                rel_heights = c(1.5, 0.5, 1),
                                subplots = 1:3,   
                                ES_geom = "line", 
                                pvalue_table = T) 
ggsave(gseap1, filename = './GSEA/gobp.pdf', width =10, height =8)


gseap2 <- gseaplot2(kk_gse,
                    up_gsea$ID,
                    title = "UP_GSEA_all",
                    color = "red",
                    base_size = 20,
                    rel_heights = c(1.5, 0.5, 1),
                    subplots = 1:3, 
                    ES_geom = "line",
                    pvalue_table = T) 
ggsave(gseap2, filename = "GSEA_up_all.pdf",width =12,height =12)

##GOBP
rm(list=ls())
load("./GSEA/GSEA_result.RData")
all<-GO_kk_entrez@result
bp<-all[all$ONTOLOGY=="BP",]

bp$log10Pvalue <- -log10(bp$pvalue)
up <- bp[bp$NES > 0, ]
down <- bp[bp$NES < 0, ]

up_sorted <- up[order(-up$log10Pvalue), ]
down_sorted <- down[order(-down$log10Pvalue), ]
write.csv(up_sorted,"./GSEA/up_go.csv")
write.csv(down_sorted,"./GSEA/down_go.csv")

top10up<-up_sorted[1:10,]
top10down<-down_sorted[1:10,]
top10up$type<-"Activated"
top10down$type<-"Suppressed"

data<-rbind(top10up,top10down)

data <- data %>%
  arrange(log10Pvalue) %>%
  mutate(Description = factor(Description, levels = unique(Description)))


p <- ggplot(data, aes(x = Description, y = log10Pvalue, size = log10Pvalue, color = NES)) +
  geom_point(alpha = 1) +
  scale_color_gradient(low = "blue", high = "red") +  
  scale_size_continuous(range = c(1, 7)) +  
  facet_grid(. ~ type, scales = "free", space = "fixed") +  
  coord_flip() + 
  theme_bw() +  
  theme(    
    plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
    strip.text.y = element_text(size = 20, color = "black", face = "bold"),  
    legend.position = "right",
    legend.title = element_text(size = 15, color = "black", face = "bold"),
    legend.text = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12,  color = "black", face = "bold"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, color = "black", face = "bold")
  ) +  
  labs(x = "Description", y = "-log10(pvalue)", size = "-log10(pvalue)", color = "NES")  

p
ggsave("./GSEA/GObptop10.pdf", p, w=10, h=6)

##Intersection of differential genes and aging genes------------------------------------------
rm(list=ls())
age_gene<-read.csv("./genage_human.csv",row.names = 1)
DEG<-read.csv("./DEG/y_old_DEG.csv",row.names = 1)
DEG <- DEG %>% dplyr::filter(threshold %in% c("Up", "Down"))
DEG_up<-DEG %>% dplyr::filter(threshold %in% c("Up"))
DEG_down<-DEG %>% dplyr::filter(threshold %in% c("Down"))


venn.diagram(
  x = list(
    age_gene$symbol,
    DEG$gene,
    DEG_up$gene,
    DEG_down$gene
  ),
  scaled = F, 
  alpha = 0.8, 
  lwd = 1,
  lty = 1,
  col = c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF"),
  label.col = 'black', 
  cex = 2, 
  fontface = "bold", 
  category.names = c("age_gene$symbol", "DEG$gene", "DEG_up$gene", "DEG_down$gene"), 
  output = TRUE,
  fill = c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF"), 
  filename = './plot/age_gene.png', 
  imagetype = "png", 
  resolution = 400,
  compression = "lzw" 
)

df <- Reduce(intersect, list(SE_as$gene, RI_as$gene, MXE_as$gene,A3SS_as$gene,A5SS_as$gene))

venn.diagram(
  x = list(
    age_gene$symbol,
    DEG$gene
  ),
  scaled = F,
  alpha = 0.8, 
  lwd = 1,
  lty = 1,
  col = c('#F1CC74','#AED2E2'), 
  cex = 2, 
  fontface = "bold", 
  category.names = c("age_gene$symbol", "DEG$gene"),  output = TRUE,
  fill = c('#F1CC74','#AED2E2'), 
  filename = './DEG/plot/age_deg.png',   imagetype = "png", 
  resolution = 400,   compression = "lzw" )
df <- Reduce(intersect, list(age_gene$symbol, DEG$gene))


##up
venn.diagram(
  x = list(
    age_gene$symbol,
    DEG_up$gene
  ),
  scaled = F, 
  alpha = 0.8, 
  lwd = 1,
  lty = 1,
  col = c('#F1CC74','#AED2E2'),   cex = 2,   fontface = "bold",   category.names = c("age_gene$symbol", "DEG_up$gene"),  output = TRUE,
  fill = c('#F1CC74','#AED2E2'),   filename = './plot/age_updeg.png', # 文件保存
  imagetype = "png", 
  resolution = 400,  compression = "lzw" 
)
df <- Reduce(intersect, list(age_gene$symbol, DEG_up$gene))

##down
venn.diagram(
  x = list(
    age_gene$symbol,
    DEG_down$gene
  ),
  scaled = F,   alpha = 0.8,
  lwd = 1,
  lty = 1,
  col = c('#F1CC74','#AED2E2'),   cex = 2, 
  fontface = "bold",
  category.names = c("age_gene$symbol", "DEG_down$gene"),   output = TRUE,
  fill = c('#F1CC74','#AED2E2'), 
  filename = '/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/DEG/plot/age_downdeg.png',   imagetype = "png", 
  resolution = 400,   compression = "lzw" )
df <- Reduce(intersect, list(age_gene$symbol, DEG_down$gene))

###10次随机采样----------------------------------------------
c<-read.csv("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/DEG/metadata.csv",row.names = 1)
table(c$dataset)[order(names(table(c$dataset)))]
meta1<-c[order(c$dataset),]

meta1 <- meta1 %>%
  dplyr::mutate(type = ifelse(Age >= 60, "old", "young"))
meta1<-meta1[,c(5,6)]
write.csv(meta1,"/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/file/metadata.csv")
write.csv(y_old_count11,"/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/file/count.csv")
df1<-y_old_count11[,rownames(meta1)]

###去批次
design <- model.matrix(~0+meta1$type,data = df1)
rB_Expr<-log2(df1+1)

rB_Expr <- removeBatchEffect(rB_Expr,batch = factor(meta1$dataset),design = design)
write.table(rB_Expr,
            file = "/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/file/time1.tab",
            quote = FALSE,sep="\t",row.names = TRUE)


###做差异
boxplot(rB_Expr,outline=FALSE,notch=T,col=factor(meta1$type),las=2)
exprSet<-normalizeBetweenArrays(rB_Expr)
boxplot(exprSet,outline=FALSE,notch=T,col=factor(meta1$type),las=2)
range(exprSet)
#exprSet<-log2(exprSet+1)
colnames(design)<-c("old","young")
contr.matrix <- makeContrasts(
  oldvsyoung=old-young, #谁在前谁比谁
  levels = colnames(design))


lfit <- lmFit(exprSet, design)
vfit <- contrasts.fit(lfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
TRD_fit<-topTable(efit,sort.by = "P",adjust="BH" , coef = 'oldvsyoung',number=Inf)


TRD_fit$threshold <- as.factor(ifelse(TRD_fit$P.Value< 0.05 & abs(TRD_fit$logFC) >= log2(2),
                                      ifelse(TRD_fit$logFC > log2(2) ,'Up','Down'),'Not'))

write.csv(TRD_fit,"/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/file/time1_DEG.csv")

##循环的代码
library(limma)

# 读入表达矩阵和meta信息
df1 <- read.table("df1.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta1 <- read.table("meta1.txt", header = TRUE, row.names = 1, sep = "\t")

# 保留老年样本固定的 SRR 编号
old_samples <- rownames(meta1[meta1$type == "old", ])

# 输出路径
output_dir <- "/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/file"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2025) # 固定随机种子以便复现

for (i in 1:10) {
  # 随机抽样 young 中的样本
  young_pool <- rownames(meta1[meta1$type == "young", ])
  young_samples <- sample(young_pool, 80)
  
  # 构建新的 meta 表
  selected_samples <- c(old_samples, young_samples)
  meta_sub <- meta1[selected_samples, , drop = FALSE]
  
  # 取出对应的表达矩阵，确保顺序一致
  df_sub <- df1[, rownames(meta_sub)]
  
  # 去批次效应
  design <- model.matrix(~0 + meta_sub$type)
  colnames(design) <- c("old", "young")
  
  rB_Expr <- log2(df_sub + 1)
  rB_Expr <- removeBatchEffect(rB_Expr, batch = factor(meta_sub$dataset), design = design)
  
  # 写入去批次后的表达矩阵
  write.table(rB_Expr,
              file = sprintf("%s/time%d.tab", output_dir, i),
              quote = FALSE, sep = "\t", row.names = TRUE)
  
  # 画箱图（节省空间，不输出图像，可注释掉）
  # pdf(sprintf("%s/time%d_boxplot.pdf", output_dir, i))
  # boxplot(rB_Expr, outline = FALSE, notch = TRUE, col = factor(meta_sub$type), las = 2)
  # dev.off()
  
  # 差异表达分析
  exprSet <- normalizeBetweenArrays(rB_Expr)
  # boxplot(exprSet, outline = FALSE, notch = TRUE, col = factor(meta_sub$type), las = 2)
  
  contr.matrix <- makeContrasts(oldvsyoung = old - young, levels = colnames(design))
  lfit <- lmFit(exprSet, design)
  vfit <- contrasts.fit(lfit, contrasts = contr.matrix)
  efit <- eBayes(vfit)
  TRD_fit <- topTable(efit, sort.by = "P", adjust = "BH", coef = 'oldvsyoung', number = Inf)
  
  TRD_fit$threshold <- as.factor(
    ifelse(TRD_fit$P.Value < 0.05 & abs(TRD_fit$logFC) >= log2(2),
           ifelse(TRD_fit$logFC > log2(2), 'Up', 'Down'),
           'Not'))
  
  # 写入 DEG 表
  write.csv(TRD_fit, file = sprintf("%s/time%d_DEG.csv", output_dir, i))
  
  cat(sprintf("Round %d complete...\n", i))
}



##10次随机采样的折线图或条形图--------------------------
# 1. 所有文件路径与组名
prefix <- "/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/file"
unrandom_file <- "/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/DEG/y_old_DEG.csv"

group_names <- c("unrandom", paste0("random", 1:10))
file_paths <- c(unrandom_file, file.path(prefix, paste0("time", 1:10, "_DEG.csv")))

# 2. 读取每个文件，统计 threshold 中 Up 和 Down 的数量
get_counts <- function(file, group) {
  df <- read.csv(file)
  df %>%
    filter(threshold %in% c("Up", "Down")) %>%
    count(threshold) %>%
    mutate(group = group)
}

# 合并所有结果
count_list <- mapply(get_counts, file_paths, group_names, SIMPLIFY = FALSE)
count_df <- do.call(rbind, count_list)

# 整理成标准格式
count_df <- count_df %>%
  rename(Direction = threshold, Count = n) %>%
  mutate(group = factor(group, levels = group_names))  # 保证顺序

# 设置统一主题
my_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),           # 去掉网格线
    axis.line = element_line(color = "black", size = 0.6),  # 横纵坐标线
    axis.text = element_text(color = "black"),              # 黑色坐标文字
    axis.title = element_text(color = "black"),             # 黑色坐标标题
    plot.title = element_text(color = "black", hjust = 0.5) # 黑色图标题，居中
  )

# 设置 y 轴范围和间距（0 到 1550，间距 310）
y_breaks <- seq(0, 1550, length.out = 5)

# 条形图
p1 <- ggplot(count_df, aes(x = group, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  my_theme +
  labs(x = "", y = "Gene Count", title = "Up/Down Gene Count per Comparison") +
  scale_fill_manual(values = c("Up" = "#e41a1c", "Down" = "#377eb8")) +
  scale_y_continuous(
    breaks = y_breaks,
    limits = c(0, 1550),
    labels = scales::number_format(accuracy = 1)  # 去掉小数点
  )

# 折线图
p2 <- ggplot(count_df, aes(x = group, y = Count, color = Direction, group = Direction)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  my_theme +
  labs(x = "", y = "Gene Count", title = "Up/Down Gene Count per Comparison") +
  scale_color_manual(values = c("Up" = "#e41a1c", "Down" = "#377eb8")) +
  scale_y_continuous(
    breaks = y_breaks,
    limits = c(0, 1550),
    labels = scales::number_format(accuracy = 1)  # 去掉小数点
  )

# 5. 保存图像
ggsave("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/plot/bar.pdf", p1, width = 10, height = 6)
ggsave("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/plot/line.pdf", p2, width = 10, height = 6)

# 6. 可选：保存表格
write.csv(count_df, "up_down_gene_counts.csv", row.names = FALSE)





##10次随机取样upset图-------------------------------------
install.packages("ComplexUpset")
remotes::install_github("krassowski/complex-upset")
library(patchwork) 
library(ComplexUpset)
# ✅ 设置文件路径和组名
prefix <- "/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/file"
unrandom_file <- "/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/DEG/y_old_DEG.csv"
group_names <- c("unrandom", paste0("random", 1:10))
file_paths <- c(unrandom_file, file.path(prefix, paste0("time", 1:10, "_DEG.csv")))

# ✅ 提取 Up / Down 基因列表
up_genes <- list()
down_genes <- list()

for (i in seq_along(file_paths)) {
  df <- read.csv(file_paths[i])
  up_genes[[group_names[i]]] <- df$X[df$threshold == "Up"]
  down_genes[[group_names[i]]] <- df$X[df$threshold == "Down"]
}

sample_list_up <- up_genes
sample_list_down <- down_genes
pdf("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/plot/up_upset.pdf", width = 12, height = 6)
UpSetR::upset(
  fromList(sample_list_up),
  nsets = length(sample_list_up),  
  nintersects = 20,                
  order.by = "degree",              
  decreasing = TRUE,               
  keep.order = FALSE,
  mb.ratio = c(0.6, 0.4),
  text.scale = 2,
  sets.bar.color = c(
    "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
    "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99"
  )[1:length(sample_list_up)],
  main.bar.color = "#6a6a6a"
)

dev.off()


pdf("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/plot/down_upset.pdf", width = 12, height = 6)
UpSetR::upset(
  fromList(sample_list_down),
  nsets = length(sample_list_down),  
  nintersects = 20,                
  order.by = "degree",              
  decreasing = TRUE,               
  keep.order = FALSE,
  mb.ratio = c(0.6, 0.4),
  text.scale = 2,
  sets.bar.color = c(
    "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
    "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99"
  )[1:length(sample_list_up)],
  main.bar.color = "#6a6a6a"
)

dev.off()

###10次随机采样的均值及差异基因占比----------------------------
rm(list=ls())
my_DEG<-read.csv("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/DEG/y_old_DEG.csv")
my_DEG_up<-subset(my_DEG,my_DEG$threshold=="Up")
my_DEG_down<-subset(my_DEG,my_DEG$threshold=="Down")

dir_path <- "/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/file/"
files <- paste0("time", 1:10, "_DEG.csv")
# 初始化合集
up_union <- c()
down_union <- c()

# 每组的上调/下调基因列表
up_list <- list()
down_list <- list()

# 用于统计每组基因数
up_counts <- c()
down_counts <- c()

# 读取 10 个文件并生成合集
for (i in 1:10) {
  df <- read.csv(file.path(dir_path, files[i]))
  up_genes <- unique(df$X[which(df$threshold == "Up")])
  down_genes <- unique(df$X[which(df$threshold == "Down")])
  
  up_list[[i]] <- up_genes
  down_list[[i]] <- down_genes
  
  up_union <- union(up_union, up_genes)
  down_union <- union(down_union, down_genes)
  
  up_counts <- c(up_counts, length(up_genes))
  down_counts <- c(down_counts, length(down_genes))
}

# 输出平均数
cat("10组中上调基因数目平均值:", mean(up_counts), "\n")#1017 
cat("10组中下调基因数目平均值:", mean(down_counts), "\n")#1166

# 读取对照数据
my_DEG <- read.csv(file.path(dir_path, "y_old_DEG.csv"))
my_DEG_up <- unique(my_DEG$X[which(my_DEG$threshold == "Up")])
my_DEG_down <- unique(my_DEG$X[which(my_DEG$threshold == "Down")])

# 计算占比
up_percent <- c(length(intersect(my_DEG_up, up_union)) / length(up_union))#0.56
down_percent <- c(length(intersect(my_DEG_down, down_union)) / length(down_union))#0.56

for (i in 1:10) {
  up_percent <- c(up_percent, length(intersect(up_list[[i]], up_union)) / length(up_union))
  down_percent <- c(down_percent, length(intersect(down_list[[i]], down_union)) / length(down_union))
}

# 准备绘图数据
comparison <- c("All sample", paste0("Subset", 1:10))
plot_df <- data.frame(
  Comparison = rep(comparison, times = 2),
  Percent = c(up_percent, down_percent),
  Direction = rep(c("Up", "Down"), each = 11)
)


plot_df$Label <- paste0(round(plot_df$Percent * 100, 1), "%")

# 绘图
ggplot(plot_df, aes(x = Comparison, y = Percent, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Label), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3.5) +  # vjust 控制垂直位置
  labs(title = "Proportion of Up/Down Genes in Union Set",
       y = "Proportion", x = "Comparison Group") +
  theme_minimal()
