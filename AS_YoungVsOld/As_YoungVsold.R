rm(list = ls())
setwd("/data3/lishuhan/30GSE_health_SRR/v3AS/")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages('maser')
library(maser)
####2.Analysis of alternative splicing results.R---------------------------------------------------
path <-("./output_rmats/")
Aging <- maser(path, c("old", "young"), ftype = "JC")
Aging_filt <- filterByCoverage(Aging, avg_reads  = 5)
Aging_top <- topEvents(Aging_filt, fdr = 0.05, deltaPSI = 0.1)
splicingDistribution(Aging_top)

####AS events volcano plot------------------------------------------------------------------------
load("./Aging.Rdata")

Aging_filt@RI_stats$logFDR<-(-log(Aging_filt@RI_stats$FDR))
Aging_filt@RI_stats$geneSymbol<-Aging_filt@RI_events$geneSymbol
data <- Aging_filt@RI_stats
names(data)

data$status <- ifelse(data$FDR<0.05&data$IncLevelDifference>0.1,'up',
                      ifelse(data$FDR<0.05&data$IncLevelDifference<(-0.1),'down','no significance'))

df <- data[is.finite(data$logFDR), ]
table(df$status)

df$status <- factor(df$status,levels =  c('no significance','up','down'))
names(df) <- c("ID","PValue","FDR","deltaPSI" ,"logFDR","geneSymbol" ,"status" )
write.csv(df,"/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/file/v2_RI_vol.csv")
df_up_label<-subset(df,deltaPSI>0.32 & logFDR>30)
df_down_label<-subset(df,deltaPSI< -0.1 & logFDR>20)
label_df<-rbind(df_up_label,df_down_label)
label_df<-label_df[-c(2,6),]

df$label=""
df$label[match(label_df$ID,df$ID)]<-df$geneSymbol[match(label_df$ID,df$ID)]


df$color <- ifelse(df$status == "no significance" & df$label == "", "#bcbcbc",                     ifelse(df$status == "up" & df$label == "", "#ffab84",                         
                                                                                                          ifelse(df$status == "down" & df$label == "", "#8abddc",                                
                                                                                                                 ifelse(df$status == "up" & df$label != "", "#be0001", "#0051a6"))))

# 绘图
df <- df %>% arrange(color)
df$color <- factor(df$color, levels = c("#bcbcbc", "#ffab84", "#8abddc", "#be0001", "#0051a6"))


p <- ggscatter(df,
               x="deltaPSI", 
               y="logFDR",  
               color = "color",  
               palette = c("#bcbcbc","#ffab84","#8abddc","#be000e","#0051a6"),  
               label = df$label,  
               font.label = c(15,"plain","black"),  
               repel = T ) +     labs( 
                 title = "RI",
                 x=expression(paste('deltaPSI'),color="black", size=12,face = "bold"),       y=expression(paste(-Log, 'FDR'),color="black", size=12,face = "bold"))+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="#bcbcbc",lwd=0.6)+
  geom_hline(yintercept = -log(0.05),lty=4,col="#bcbcbc",lwd=0.6)+
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
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5))

p
ggsave(filename = './plot/AS_volcano.pdf',p, width = 8,height = 6,dpi = 600)

###Splicing events and parental gene number statistics----------------------------------------
RI<-read.csv("./RI_vol.csv")
SE<-read.csv("./SE_vol.csv")
MXE<-read.csv("./MXE_vol.csv")
A3SS<-read.csv("./A3SS_vol.csv")
A5SS<-read.csv("./A5SS_vol.csv")

SE$AS_type <- "SE"
RI$AS_type <- "RI"
MXE$AS_type <- "MXE"
A3SS$AS_type <- "A3SS"
A5SS$AS_type <- "A5SS"

df_combined <- bind_rows(SE, RI, MXE, A3SS, A5SS)
df_combined$differ <- ifelse(df_combined$status == "old", 1, 0)

df_combined<-subset(df_combined,df_combined$status %in% c("old","young"))
df<-melt(table(df_combined$status,df_combined$AS_type))

p <- ggplot(df, aes(x = Var1, y = value, fill = Var2)) +
  #geom_bar(position = "do") +
  facet_wrap(~ Var1, ncol = 2, scales = "free_x") +
  geom_bar(stat = 'identity', position = 'dodge',
           width = 0.8, color = 'black') +
  geom_text(aes(label = value), size = 4,
            position = position_dodge(width = 0.8),
            vjust = -0.3) +
  labs(x = NULL) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(colour = 'black', size = 12, face = "bold"), 
    axis.text.y = element_text(colour = 'black', size = 12, face = "bold"),  
    strip.text = element_text(face = "bold", size = 12),  
    legend.title = element_text(face = "bold", size = 12),  
    legend.text = element_text(face = "bold", size = 10),  
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5) 
  ) +
  scale_fill_manual(values = rep(c('#D55E00', '#0072B2', '#F0E442',
                                   '#009E73', '#9270ca', '#CC79A7'), 3)) +
  labs(fill = "Type", y = "Number", x = "AS Events")

p
##gene
data<-melt(table(df_combined$status,df_combined$AS_type))
data$ga<-paste0(data$Var1,data$Var2)
df_combined$ga<-paste0(df_combined$status,df_combined$AS_type,df_combined$geneSymbol)
df_combined1<-df_combined[!duplicated(df_combined$ga),]
data<-melt(table(df_combined1$status,df_combined1$AS_type))

p1<-ggplot(data,aes(x=Var1,y=value,fill=Var2))+
  #geom_bar(position = "do") +
  facet_wrap(~Var1, ncol = 2, scales = "free_x") + 
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  geom_text(aes(label=value),size=4,
            position = position_dodge(width = 0.8), 
            vjust=-0.3)+ 
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+
  theme(
    axis.text.x = element_text(colour = 'black', size = 12, face = "bold"),  
    axis.text.y = element_text(colour = 'black', size = 12, face = "bold"),  
    strip.text = element_text(face = "bold", size = 12),  
    legend.title = element_text(face = "bold", size = 12), 
    legend.text = element_text(face = "bold", size = 10),  
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)  
  ) +
  scale_fill_manual(values = rep(c('#D55E00','#0072B2','#F0E442',
                                   '#009E73','#9270ca','#CC79A7'),3))+
  labs(fill="Type",y="Number",x="AS Genes")


(p+p1)
ggsave(filename = './plot/as_gene_sum.pdf',(p+p1),width = 15,height = 6,dpi=600)


###The proportion of splicing events and the number of parental genes-----------------
df_combined <- bind_rows(SE, RI, MXE, A3SS, A5SS)
df_combined<-subset(df_combined,df_combined$status %in% c("old","young"))

as_type_counts <- as.data.frame(table(df_combined$AS_type))
colnames(as_type_counts) <- c("AS_type", "Count")

total_count <- nrow(df_combined)


as_type_counts$Proportion <- (as_type_counts$Count / total_count) * 100



custom_colors <- c("A3SS" = "#D55E00", "A5SS" = "#0072b2", "MXE" = "#F0E442", "RI" = "#009E73", "SE" = "#9270ca")

p<-ggplot(data = as_type_counts, aes(x = reorder(AS_type, Count), y = Count, fill = AS_type)) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = Count), hjust = -0.3, color = "black", size = 6, fontface = "bold") +
  geom_text(aes(label = paste0(round(Proportion, 2), "%")), hjust = 1.3, color = "black", size = 6, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  xlab("AS Type") + 
  ylab("Number") +
  scale_y_continuous(limits = c(0, 4300), breaks = seq(0, 4300, 860)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.ticks = element_line(size = 1.5),
    axis.line = element_line(size = 1.5),
    panel.grid = element_blank()
  )

ggsave(filename = '/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/plot/precent.pdf',
       p,
       width = 8,height = 6,dpi=600)

##gene

data<-melt(table(df_combined$status,df_combined$AS_type))
data$ga<-paste0(data$Var1,data$Var2)
df_combined$ga<-paste0(df_combined$status,df_combined$AS_type,df_combined$geneSymbol)
df_combined1<-df_combined[!duplicated(df_combined$ga),]data<-melt(table(df_combined1$status,df_combined1$AS_type))

summarized_df <- data %>%
  group_by(Var2) %>%
  summarise(total_value = sum(value))

as_type_counts <- as.data.frame(table(summarized_df$Var2))
colnames(as_type_counts) <- c("AS_type", "Count")


# 绘图

custom_colors <- c("A3SS" = "#D55E00", "A5SS" = "#0072b2", "MXE" = "#F0E442", "RI" = "#009E73", "SE" = "#9270ca")


total_sum <- sum(summarized_df$total_value)


summarized_df$Proportion <- (summarized_df$total_value / total_sum) * 100


p <- ggplot(data = summarized_df, aes(x = reorder(Var2, total_value), y = total_value, fill = Var2)) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = total_value), hjust = -0.3, color = "black", size = 6, fontface = "bold") +
  geom_text(aes(label = paste0(round(Proportion, 2), "%")), hjust = 1.3, color = "black", size = 6, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  xlab("AS Type") + 
  ylab("Number") +
  scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, 600)) +
  #scale_y_continuous(limits = c(0, max(summarized_df$total_value) * 1.1), breaks = seq(0, max(summarized_df$total_value) * 1.1, by = 500)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.ticks = element_line(size = 1.5),
    axis.line = element_line(size = 1.5),
    panel.grid = element_blank()
  )


ggsave(filename = './plot/gene_precent.pdf', p,width = 8,height = 6,dpi=600)

###upset----------------------------------------------------------------------------------
se_deduped <- data %>%
  dplyr::filter(AS_type == "SE") %>%
  distinct(gene_name, chr, strand.x, start, end, differ, AS_type, .keep_all = TRUE)


result <- data %>%
  dplyr::filter(AS_type != "SE") %>%
  bind_rows(se_deduped)

write.csv(result,"./AS.csv")
data<-read.csv("/AS.csv")
sample_list <- list(A3SS = alldata$geneSymbol[alldata$AS_type %in% "A3SS"], A5SS = alldata$geneSymbol[alldata$AS_type %in% "A5SS"], 
                    SE = alldata$geneSymbol[alldata$AS_type %in% "SE"],
                    RI = alldata$geneSymbol[alldata$AS_type %in% "RI"],MXE = alldata$geneSymbol[alldata$AS_type %in% "MXE"])
png("/data3/lishuhan/30GSE_health_SRR/v2AS/plot/y_old_upset.png", width = 15*600, height = 7*600, res = 600)
upset1<-upset(fromList(sample_list), 
              nsets = 5,     
              nintersects = 60, 
              order.by = c("freq","degree"),
              #group.by = c("sets"), cutoff = 7, 
              keep.order = F, 
              mb.ratio = c(0.6,0.4),   
              text.scale = 2, 
              #sets.bar.color = c("#ff6600","#a1cc04","#00b0f0","#e10fd7","#7030a0"),
              sets.bar.color = c("#6aa5df","#6db385","#70ab37","#a48b30","#8e77fb"),
              main.bar.color = c("#6a6a6a")
)
upset1
dev.off()

###sample information plot--------------------------------------------------------------------------
rm(list = ls())
a<-read.csv("./df_0_17.csv")
b<-read.csv("./df_17_59.csv")
c<-read.csv("./df_60_100.csv")
colnames(c)[2]<-"Sample_Name"
d<-rbind(a,b,c)
write.csv(d,"./SRR_info/health.csv")


freq <- table(d$dataset)

levels <- names(freq)[order(freq, decreasing = TRUE)]

d$dataset <- factor(d$dataset, levels = levels)
d$Age <- as.numeric(as.character(d$Age ))

color_vector <- c("Healthy" = "#30ba43")
g=ggplot(data = d,aes(x = dataset, y = Age , fill = status))+
  geom_point(aes(fill=status),size=4,shape=21,color ="grey20",
             position = position_jitter(width=0.2,height=0.1))+
  #geom_jitter(mapping = aes(x = dataset, y = Age, color = status), width = 0.2) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5)) +
  scale_fill_manual(values = color_vector) +  
  xlab("Dataset") +
  ylab("Age(years)")+
  # labs(colour = "Status")+
  theme_minimal() +
  theme(
    panel.background = element_blank(),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"), 
    axis.text.x = element_text(angle = 60, hjust = 1,face = "bold",color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1,face = "bold",color = "black")
  )+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.line = element_line(size = 1, colour = 'black'),
        axis.ticks = element_line(size = 1, colour = 'black'))+
  theme(
    #panel.border = element_rect(colour = "black", fill=NA, size=1.5),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    text = element_text(size=12, face="bold"),  
    plot.title = element_text(size=16, face="bold"),  
    legend.text = element_text(size=16, face="bold"), 
    axis.title = element_text(size=14, face="bold"), 
    axis.text = element_text(size=14, face="bold")  
  )
g

ggsave(plot=g,"./v2_info.pdf"
       ,width = 10,height=6,dpi=300)

p<-ggplot(a, aes(x = Age)) + 
  stat_density(geom = "area", alpha = 0.3, color = "#5c93bd",fill="#5c93bd") +  
  scale_x_continuous(limits = c(8, 17), breaks = seq(8, 17, by = 1)) +  
  labs(x = "Age (years)", y = "Density", title = "Density Plot of Age") +    theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    text = element_text(size=12, face="bold"), 
    plot.title = element_text(size=16, face="bold"),  
    legend.text = element_text(size=14, face="bold"),  
    axis.title = element_text(size=14, face="bold",colour = "black"),  
    axis.text = element_text(size=14, face="bold",colour = "black")    )
p
ggsave(plot=p,"./plot/v2_young_info.pdf",
       width = 10,height=6,dpi=300)

p2<-ggplot(c, aes(x = Age)) + 
  stat_density(geom = "area", alpha = 0.3, fill="#ff8670",color = "#ff8670") +    scale_x_continuous(limits = c(60, 87), breaks = seq(60, 87, by = 1)) +  
  labs(x = "Age (years)", y = "Density", title = "Density Plot of Age") + 
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    text = element_text(size=12, face="bold"),  
    plot.title = element_text(size=16, face="bold"),  
    legend.text = element_text(size=14, face="bold"),  
    axis.title = element_text(size=14, face="bold",colour = "black"),  
    axis.text = element_text(size=14, face="bold",colour = "black")  
  )
p2
ggsave(plot=p2,"./plot/v2_old_info.pdf",
       width = 10,height=6,dpi=300)
(p)/(p2)
ggsave(plot=(p)/(p2),
       "./plot/merge_yo_info.pdf",width = 10,height=12,dpi=300)


##Map of the distribution of splicing events on chromosomes---------------------------------
rm(list=ls())
chr_size = read.table("./DEG/hg38.chrom.sizes.txt",sep = '\t',header = F,row.names = 1)
chr_size$chr = rownames(chr_size)
colnames(chr_size)[1]='Total_length_bp'
df<-read.csv("./ALL_AS_chr.csv",row.names = 1)
library(ggrepel)
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
df$chr <- factor(df$chr, levels = chr_levels)
p <- ggplot()+
  geom_jitter(data = df,
              aes(x = chr, y = deltaPSI),
              size = 0.85,
              width =0.4)
p

df_positive <- df %>% dplyr::filter(deltaPSI > 0)
df_negative <- df %>% dplyr::filter(deltaPSI < 0)


max_deltaPSI_per_chr <- df_positive %>%
  group_by(chr) %>%
  summarise(max_deltaPSI = max(deltaPSI))


min_deltaPSI_per_chr <- df_negative %>%
  group_by(chr) %>%
  summarise(min_deltaPSI = min(deltaPSI))


chr_levels <- sort(unique(c(df_positive$chr, df_negative$chr)))
chr_to_num <- setNames(seq_along(chr_levels), chr_levels)

max_deltaPSI_per_chr <- max_deltaPSI_per_chr %>%
  mutate(x = chr_to_num[chr]) %>%
  arrange(x)

min_deltaPSI_per_chr <- min_deltaPSI_per_chr %>%
  mutate(x = chr_to_num[chr]) %>%
  arrange(x)


dfbar <- data.frame(
  x = max_deltaPSI_per_chr$x,
  y = max_deltaPSI_per_chr$max_deltaPSI
)
dfbar$x<-rownames(dfbar)
dfbar1 <- data.frame(
  x = min_deltaPSI_per_chr$x,
  y = min_deltaPSI_per_chr$min_deltaPSI
)
dfbar1$x<-rownames(dfbar1)

p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1


color_vector <- c("A3SS" = "#D55E00", "A5SS" = "#0072b2", "MXE" = "#F0E442", "RI" = "#009E73", "SE" = "#9270ca") 

p2 <- ggplot() +
  geom_col(data = dfbar,
           mapping = aes(x = x, y = y),
           fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = dfbar1,
           mapping = aes(x = x, y = y),
           fill = "#dcdcdc", alpha = 0.6) +
  geom_jitter(data = df,
              aes(x = factor(chr, levels = paste0("chr", c(1:22, "X", "Y"))), y = deltaPSI, color = AS_type),
              size = 0.85, width = 0.4) +
  scale_color_manual(values = color_vector) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1)
  )
p2

dfcol <- data.frame(
  x = as.numeric(chr_levels),
  y = 0,
  label = levels(chr_levels)
)


mycol <- colorRampPalette(c("#E64B35", "#4DBBD5", "#00A087", 
                            "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148"))(length(dfcol$label))

p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.18,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3

p6 <- p3+
  labs(x="chr",y="delPSI")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =5,
            color ="white")
p6
：
p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )+
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5))
p7
ggsave(filename = './plot/chr.pdf',p7, width = 15,height = 6,dpi = 600)

###PPI plot data construction----------------------------------------------------------------------
SE<-read.csv("./file/SE_vol.csv",row.names = 1)
RI<-read.csv("./file/RI_vol.csv",row.names = 1)
MXE<-read.csv("./file/MXE_vol.csv",row.names = 1)
A3SS<-read.csv("./file/A3SS_vol.csv",row.names = 1)
A5SS<-read.csv("./file/A5SS_vol.csv",row.names = 1)
df<-rbind(SE,RI,MXE,A3SS,A5SS)
library(dplyr)
result <- df %>%
  filter(status %in% c("old", "young")) %>%  
  arrange(desc(logFDR)) %>%  
  distinct(geneSymbol, .keep_all = TRUE) %>%  
  slice_head(n = 200)  
##GO
load(file = "./Genome/enrichment_annotation/human_enrichment.RData")
term2gene <- gobp[,c(1,5)]
term2name <- gobp[,c(1,2)]
DEG_go_bp <- enricher(gene = result$geneSymbol,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                      TERM2GENE = term2gene,TERM2NAME = term2name)
DEG_go_bp<-data.frame(DEG_go_bp@result)
write.csv(DEG_go_bp,"./DEG_go_bp.csv")
go_df<-DEG_go_bp[c(DEG_go_bp$pvalue<0.05),]
go_df1<-DEG_go_bp[c(DEG_go_bp$p.adjust<0.05),]

install.packages("httr2")
packageVersion('httr2')
install.packages("fanyi")
library(fanyi)
set_translate_option(appid = '20231223001918987', key = '2jzXfF7T1Eloro8WaPjb', source = 'baidu')
baidu_translate(go_df$Description, from = "en", to = "zh")
a<-baidu_translate(go_df$Description, from = "en", to = "zh")
protein<-go_df[c(1,3,4,7,34,35,33,23,20,143,73),]
aging<-go_df[c(9,12,84,31,34,46,49,81,98,197,83,137),]
immune<-go_df[c(5,6,7,9,12,11,21,22,29,63,41,43,45,113,115,157,153),]
protein$type<-"protein"
aging$type<-"aging"
immune$type<-"immune"
GO_df<-rbind(protein,aging,immune)
GO_df<-GO_df[!duplicated(GO_df$geneID),]

term2gene1<-term2gene[term2gene$GO %in% GO_df$ID,]
b<-c()
d<-c()
und<-c()
for (i in 1:length(GO_df$ID)) {
  name1<-term2gene1$SYMBOL[grep(GO_df$ID[i],term2gene1$GO)]
  
  und<-c(und,i)
  for (j in GO_df$ID[-i]) {
    a<-sum(name1 %in% term2gene1$SYMBOL[grep(j,term2gene1$GO)])
    if(a>0){
      b<-c(b,j)
      d<-c(d,GO_df$ID[i])
    }
  }
}
df<-data.frame(source=d,target=b)

df1<-merge(df,GO_df[,1:2],by.x="source",by.y="ID")
colnames(df1)<-c("source_ID","target_ID","source_desc")
df1<-merge(df1,GO_df[,1:2],by.x="target_ID",by.y="ID")
colnames(df1)[4]<-"target_desc"
count<-melt(table(df1$source_ID))
df1<-merge(df1,count,by.x="source_ID",by.y="Var1")
df1<-merge(df1,GO_df,by.x="source_ID",by.y="ID")
df2<-df1[,c(1,2,3,4,13)]
write.csv(df2,"./file/GO_ppi.csv")
write.table(df2,"./GO_ppi.txt",sep = "\t",row.names = F)


term2gene1<-term2gene[term2gene$GO %in% GO_df$ID,]
b<-c()
d<-c()
und<-c()
for (i in 1:length(GO_df$ID)) {
  name1<-unlist(strsplit(GO_df$geneID[i],"/"))
  und<-c(und,i)
  k<-c()
  for (j in GO_df$ID[-und]) {
    a<-sum(name1 %in% unlist(strsplit(GO_df$geneID[grep(j,GO_df$ID)],"/")))
    ####unlist(strsplit(GO_df$geneID[grep(j,GO_df$ID)],"/"))
    ####term2gene1$SYMBOL[grep(j,term2gene1$ID)]
    if(a>0){
      b<-c(b,j)
      d<-c(d,GO_df$ID[i])
    }
  }
}
df<-data.frame(source=d,target=b)
df1<-merge(df,GO_df[,1:2],by.x="source",by.y="ID")
colnames(df1)<-c("source_ID","target_ID","source_desc")
df1<-merge(df1,GO_df[,1:2],by.x="target_ID",by.y="ID")
colnames(df1)[4]<-"target_desc"



df1<-merge(df1,GO_df,by.x="source_ID",by.y="ID")
df2<-df1[,c(1,2,3,4,12,13)]


tar<-unique(df2$target_ID[!df2$target_ID %in% df2$source_ID])
num1=dim(df2)[1]+length(tar)
df3<-data.frame(source_ID=tar,
                target_ID="na",
                source_desc=GO_df[tar,"Description"],
                target_desc="na",
                Count=GO_df[tar,"Count"],
                type=GO_df[tar,"type"])
df2<-rbind(df2,df3)

write.table(df2,"./GO_ppi.txt",sep = "\t",row.names = F)

####Splicing factor---------------------------------------------------------------------------
DEG<-read.csv("./DEG/y_old_DEG.csv",row.names = 1)
DEG_df <- DEG %>%
  filter(threshold %in% c("Up", "Down"))
AS_factor<-read.csv("./file/AS_factor.csv")
DEG_AS_df<-merge(DEG_df,AS_factor,by.x="gene",by.y="gene_name")
write.csv(AS_factor,"./file/AS_factor1.csv")

SE_as<-read.csv("./SE_factor.csv")
RI_as<-read.csv("./RI_factor.csv")
MXE_as<-read.csv("./MXE_factor.csv")
A3SS_as<-read.csv("./A3SS_factor.csv")
A5SS_as<-read.csv("./A5SS_factor.csv")

venn.diagram(x=list(SE_as$gene,RI_as$gene,MXE_as$gene,A3SS_as$gene,A5SS_as$gene),
             
             scaled = F,              
             alpha= 0.8,              
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC"), 
             
             label.col ='black' , 
             abel.col=c('#FFFFCC','#CCFFFF',......)             
             cex = 2,              
             fontface = "bold", 
             
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC"),              
             category.names = c("SE_AS_factor", "RI_AS_factor","MXE_AS_factor","A3SS_AS_factor","A5SS_AS_factor") ,             
             cat.dist = c(0.2, 0.2, 0.2, 0.2, 0.2), 
             
             cat.pos = c(0, -10, 240, 120, 20), 
             
             cat.cex = 1,              
             cat.fontface = "bold",               
             cat.col=c('#f8c527','#29a3e7',"#de0067","#5e1b84", "#39ab43"),   
             
             cat.default.pos = "outer",               
             output=TRUE,
             
             filename='./plot/five.png',
             
             imagetype="png",  
             
             resolution = 400,  
             
             compression = "lzw"
             
)

grid.draw(data)

##Splcing factor correlation----------------------------------------------------------------------
rm(list=ls())
load("./rmaps/AS_correlation_list.RData")

filtered_results_list <- list()

for (gene_name in names(AS_results_list)) {
  df <- AS_results_list[[gene_name]]
  
  significant <- df[df$p_value < 0.05, ]
  
  significant <- significant[significant$correlation > 0, ]
  
  filtered_results_list[[gene_name]] <- significant
}

row_counts <- sapply(filtered_results_list, nrow)
sorted_row_counts <- sort(row_counts, decreasing = TRUE)
top10_filtered_results <- filtered_results_list[top10_names]


top10_max_correlation_genes <- list()


for (gene_name in top10_names) {
  df <- top10_filtered_results[[gene_name]]
  max_correlation_row <- df[which.max(df$correlation), ]
  top10_max_correlation_genes[[gene_name]] <- max_correlation_row$df2_geneName
}

unique_genes <- unique(unlist(top10_max_correlation_genes))


heatmap_matrix <- matrix(0, nrow = length(unique_genes), ncol = length(top10_names), 
                         dimnames = list(unique_genes, top10_names))


for (gene_name in top10_names) {
  df <- top10_filtered_results[[gene_name]]
  for (unique_gene in unique_genes) {
    if (unique_gene %in% df$df2_geneName) {
      correlation_value <- df$correlation[df$df2_geneName == unique_gene]
      heatmap_matrix[unique_gene, gene_name] <- correlation_value
    }
  }
}

rownames(heatmap_matrix)<-paste("A3SS", rownames(heatmap_matrix), sep = "_")

heatmap_matrix <- heatmap_matrix[, factor(colnames(heatmap_matrix), levels = top10_names)]

p <- pheatmap(heatmap_matrix,
              clustering_distance_rows = "correlation",
              clustering_distance_cols = "correlation",
              cluster_rows = TRUE,
              cluster_cols = F,
              fontsize = 12,          
              fontsize_row = 10,       
              fontsize_col = 10,       
              angle_col = 45, 
              color = colorRampPalette(c("#f7f7f7","#c0688c"))(50),               annotation_legend = TRUE
)

p
ggsave(plot=p,"./rmaps/AS.pdf",width = 4,height = 4)


###Calculation of correlation between splicing factors and gene expression values------
rm(list=ls())

df <- read.delim("./DEG/rB_Expr.tab"
                 , header = TRUE, sep = "\t")
df17<-read.csv("./df_0_17.csv",row.names = 1)
df60<-read.csv("./df_60_100.csv",row.names = 1)
colnames(df60)[1]<-"Sample_Name"
meta<-rbind(df60,df17)
df <- df[, rownames(meta)]

load("./Aging.Rdata")

Aging_filt@A3SS_stats$logFDR<-(-log(Aging_filt@A3SS_stats$FDR))
Aging_filt@A3SS_stats$geneSymbol<-Aging_filt@A3SS_events$geneSymbol
data <- Aging_filt@A3SS_stats
names(data)

data$status <- ifelse(data$FDR<0.05&data$IncLevelDifference>0.1,'up',
                      ifelse(data$FDR<0.05&data$IncLevelDifference<(-0.1),'down','no significance'))

A3SS_data <- data %>%
  dplyr::filter(status %in% c("up", "down"))


A3SS_psi<-as.data.frame(Aging_top@A3SS_PSI)
na_percentage <- rowSums(is.na(A3SS_psi)) / ncol(A3SS_psi)
df_filtered <- A3SS_psi[na_percentage <= 0.8, ]
col_means <- apply(df_filtered, 2, function(x) mean(x, na.rm = TRUE))
df_imputed <- apply(df_filtered, 2, function(x) ifelse(is.na(x), col_means, x))
A3SS_psi <- as.data.frame(df_imputed)



A5SS_psi<-as.data.frame(Aging_top@A5SS_PSI)
na_percentage <- rowSums(is.na(A5SS_psi)) / ncol(A5SS_psi)
df_filtered <- A5SS_psi[na_percentage <= 0.8, ]
col_means <- apply(df_filtered, 2, function(x) mean(x, na.rm = TRUE))
df_imputed <- apply(df_filtered, 2, function(x) ifelse(is.na(x), col_means, x))
A5SS_psi <- as.data.frame(df_imputed)

MXE_psi<-as.data.frame(Aging_top@MXE_PSI)
na_percentage <- rowSums(is.na(MXE_psi)) / ncol(MXE_psi)
df_filtered <- MXE_psi[na_percentage <= 0.8, ]
col_means <- apply(df_filtered, 2, function(x) mean(x, na.rm = TRUE))
df_imputed <- apply(df_filtered, 2, function(x) ifelse(is.na(x), col_means, x))
MXE_psi <- as.data.frame(df_imputed)

SE_psi<-as.data.frame(Aging_top@SE_PSI)
na_percentage <- rowSums(is.na(SE_psi)) / ncol(SE_psi)
df_filtered <- SE_psi[na_percentage <= 0.8, ]
col_means <- apply(df_filtered, 2, function(x) mean(x, na.rm = TRUE))
df_imputed <- apply(df_filtered, 2, function(x) ifelse(is.na(x), col_means, x))
SE_psi <- as.data.frame(df_imputed)

RI_psi<-as.data.frame(Aging_top@RI_PSI)
na_percentage <- rowSums(is.na(RI_psi)) / ncol(RI_psi)
df_filtered <- RI_psi[na_percentage <= 0.8, ]
col_means <- apply(df_filtered, 2, function(x) mean(x, na.rm = TRUE))
df_imputed <- apply(df_filtered, 2, function(x) ifelse(is.na(x), col_means, x))
RI_psi <- as.data.frame(df_imputed)

##Take A3SS for example
rm(list=ls())
A3SS_data$ID <- as.character(A3SS_data$ID)
id_to_gene <- setNames(A3SS_data$geneSymbol, A3SS_data$ID)
new_row_names <- paste(id_to_gene[rownames(A3SS_psi)], rownames(A3SS_psi), sep = "_")
rownames(A3SS_psi) <- new_row_names

result_list <- list()


for (row_name in rownames(A3SS_psi)) {
  
  gene_info <- strsplit(row_name, "_")[[1]]
  gene <- gene_info[1]
  id <- gene_info[2]
  psi_values <- as.numeric(A3SS_psi[row_name, ])
  
  
  if (gene %in% rownames(df)) {
    count_values <- as.numeric(df[gene, ])
    
    
    temp_df <- data.frame(
      Sample = colnames(A3SS_psi),
      PSI = psi_values,
      Count = count_values
    )
    
    
    result_list[[id]] <- temp_df
  }
}
A3SS_res<-result_list
save(A3SS_res,file="./A3SS_res.RData")

load("./A3SS_res.RData")

for (name in names(result_list)) {
  a <- result_list[[name]]
  a <- a[complete.cases(a$PSI, a$Count), ]
  p2 <- ggplot(a, aes(x = PSI, y = Count)) +
    geom_point(color = "navy", fill = "navy", shape = 21, size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", aes(alpha = 0.5), color = "navy", fill = "grey",
                formula = y ~ poly(x, 1, raw = FALSE), linetype = 2, alpha = 0.3) +
    stat_poly_eq(formula = y ~ poly(x, 1, raw = FALSE), 
                 aes(label = paste(after_stat(eq.label), after_stat(adj.rr.label), after_stat(p.value.label), sep = "~~~~")), 
                 parse = TRUE, label.x = "left") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = '#333c41', size = 12),
          legend.text = element_text(color = '#333c41', size = 12),
          legend.title = element_blank(),
          legend.position = c(0.9, 0.9)) +
    labs(x = "PSI_value", y = "count_value") +
    theme(axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"))
  
  
  ggsave(plot = p2, filename = paste0("./A3SS_plot/", name, ".pdf"), width = 4, height = 4)
}

###Interproscan---------------------------------------------------------------------------------------
rm(list=ls())

all_pro_locate <- read.delim('./AS/domain/uniprot-compressed_true_download_true_format_fasta_query__28_28taxon-2023.04.04-13.11.30.15.fasta.tsv',
                             header = T,sep='\t')
all_pro_locate_pfam <- all_pro_locate[all_pro_locate$Analysis=='Pfam',]
all_pro_locate_pfam2 <- separate(all_pro_locate_pfam,'ID',c('sp','protein_uniport','protein_uniport2'))
uniport <- as.data.frame(unique(all_pro_locate_pfam2$protein_uniport))

write.table(unique(all_pro_locate_pfam2$protein_uniport),
            file = './AS/domain/protein_uniport.txt',
            sep = '\n',
            quote = F,col.names = F,row.names = F)
change_id <- read.delim('./AS/domain/uniprot-compressed_true_download_true_format_tsv-2023.04.07-04.26.35.28.tsv')
names(change_id) <- c('protein_uniport','ensemble_uniport')
all_pro_locate_pfam4 <- merge(all_pro_locate_pfam2,change_id)
write.table(all_pro_locate_pfam4,
            file = './AS/domain/all_pro_locate_pfam4.txt',
            sep = '\t',
            quote = F,col.names = T,row.names = F)



rm(list=ls())
library(EnsDb.Hsapiens.v86)library(ensembldb)

all_pro_locate_pfam4 <- read.delim('./AS/domain/all_pro_locate_pfam4.txt')


###Take SE example
filterri_human=read.delim("./output_rmats/SE.MATS.JC.txt")
filterri_human=filterri_human[filterri_human$FDR<0.05&abs(filterri_human$IncLevelDifference)>0.1,]
filterri_human=na.omit(filterri_human)
names(filterri_human)
filterri_human$chr <- gsub('chr','',filterri_human$chr)
filterri_human2 <- filterri_human[,c(1:7)]


edb <- EnsDb.Hsapiens.v86

res_all <- NULL

for(i in 1:nrow(filterri_human[1:7,])){
  gnmRI <- GRanges(filterri_human2[i,4],
                   IRanges(start = filterri_human2[i,6],end = filterri_human2[i,7]-1))
  gnm_res <- genomeToProtein(gnmRI, edb)
  
  res <- as.data.frame(gnm_res);print(res)
  for (j in 1:dim(res)[1]) {
    res1<-res[j,,drop=F]
    if(!res$start[j]==-1){
      res1$col = i
      if(exists("res_all")==FALSE){
        res_all = res1
      }else{
        res_all <- rbind(res_all,res1)
      }}
  }
}

res_all2 <- res_all[res_all$cds_ok=='TRUE',]#323
res_all3 <- res_all2[,c(3,4,5,6,10,11)]
all_pro_locate_pfam4$ensemble_uniport <- gsub('\\.*','',all_pro_locate_pfam4$ensemble_uniport)
names(all_pro_locate_pfam4)
all_pro_locate_pfam4_2 <- all_pro_locate_pfam4[,c(3,9,10,15,16)]
all_pro_locate_pfam4_2$names <- substr(all_pro_locate_pfam4_2$ensemble_uniport,1,nchar(all_pro_locate_pfam4_2$ensemble_uniport)-1)

self_domain <- merge(res_all3,all_pro_locate_pfam4_2,by='names')
self_domain$domain_if <-  'other'

for(i in 1:nrow(self_domain)){
  self_domain[i,]$domain_if <- ifelse(self_domain[i,]$start>(self_domain[i,]$Start_location+self_domain[i,]$width+1)&(self_domain[i,]$start<self_domain[i,]$Stop_location),'TRUE','FALSE')
}
table(self_domain$domain_if)


self_domain_T <- self_domain[self_domain$domain_if=='TRUE',]
write.csv(self_domain_T,file = './interproscan/SE_domain.csv')


###Various intersection and scoring analysis
#intersection
SASP<- read.csv("./SASP_SenMayo.csv",sep = "\t")
SASP<-subset(SASP,SASP$Annotation=="SASP")
my_list <- list(DEG1$gene, SASP$Symbol)
my_list[[2]] <- unique(my_list[[2]])
names(my_list) <- c("DEGs", "SASP")


fit <- euler(my_list, shape = "circle")


plot(
  fit,
  
  fills = list(fill = c("#F98F34", "#0C4E9B"), alpha = 0.8),
  labels = list(col = "black", font = 2, fontfamily = "serif", cex = 1.5),
  edges = list(col = "black", lwd = 3, lty = 2),  
  quantities = list(type = "counts", cex = 1)  
) -> p2

p2

#scoring
rm(list=ls())
infa<-read.csv("./immune/infla_gene.csv")
rB_Expr<-read.table("./DEG/rB_Expr.tab",sep="\t")

#metadata
c<-read.csv("./DEG/metadata.csv",row.names = 1)
table(c$dataset)[order(names(table(c$dataset)))]
meta1<-c[order(c$dataset),]
meta1 <- meta1 %>%
  mutate(type = ifelse(Age >= 60, "old", "young"))
meta1<-meta1[,c(5,6)]

##GSVA
infa_genes <- list(infa$Gene) 
gsva_scores <- gsva(as.matrix(rB_Expr), infa_genes, method = "plage", kcdf = "Gaussian")
gsva_scores_vector <- as.vector(gsva_scores)

data <- data.frame(score = gsva_scores_vector, group = meta1$type)

p <- ggviolin(data, x="group", y="score", color = "group", 
              ylab="score",
              xlab="group",
              add.params = list(fill="white"),
              palette = c("#F98F34", "#0C4E9B"), 
              width=1, add = "boxplot")
#p=p+rotate_x_text(60)
p1 <- p + stat_compare_means(comparisons = list(c("old", "young")))
p1
ggsave("./SASP_score.pdf", plot = p1, width = 4, height = 3,dpi=600)


sasp_expression <- rB_Expr %>% 
  dplyr::filter(row.names(rB_Expr) %in% infa$Gene)


sasp_means <- colMeans(sasp_expression)


data_for_plot <- data.frame(
  Sample = names(sasp_means),
  GeneExp = sasp_means,
  Type = meta1$type  
)


p <- ggviolin(data_for_plot, x="Type", y="GeneExp", color = "Type", 
              ylab="Average expression of inflammatory factor sets",
              xlab="Type",
              add.params = list(fill="white"),
              palette = c("#0C4E9B", "#F98F34"), 
              width=1, add = "boxplot")

p1 <- p + stat_compare_means(comparisons = list(c("old", "young")))



ggsave("./infa_mean.pdf", p1, width = 6, height = 5,dpi=600)



###Balanced sampling---------------------------------------
rm(list=ls())
write.csv(df_combined,"/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/before.csv")
path <-("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/output_rmats1/")
Aging <- maser(path, c("old", "young"), ftype = "JC")
Aging

# A Maser object with 127386 splicing events.
# 
# Samples description: 
#   Label=young     n=80 replicates
# Label=old     n=83 replicates
# 
# Splicing events: 
#   A3SS.......... 8069 events
# A5SS.......... 5223 events
# SE.......... 92150 events
# RI.......... 4671 events
# MXE.......... 17273 events
Aging_filt <- filterByCoverage(Aging, avg_reads  = 5)
Aging_filt
# A Maser object with 76703 splicing events.
# 
# Samples description: 
#   Label=young     n=80 replicates
# Label=old     n=83 replicates
# 
# Splicing events: 
#   A3SS.......... 4395 events
# A5SS.......... 2695 events
# SE.......... 55430 events
# RI.......... 2960 events
# MXE.......... 11223 events
Aging_top <- topEvents(Aging_filt, fdr = 0.05, deltaPSI = 0.1)
Aging_top
# A Maser object with 4820 splicing events.
# 
# Samples description: 
#   Label=young     n=80 replicates
# Label=old     n=83 replicates
# 
# Splicing events: 
#   A3SS.......... 247 events
# A5SS.......... 202 events
# SE.......... 413 events
# RI.......... 1126 events
# MXE.......... 2832 events

##数目统计
# 定义事件类型
event_types <- c("RI", "SE", "MXE", "A3SS", "A5SS")

# 创建一个空列表用于存储各个结果数据框
all_dfs <- list()

# 循环处理每种事件类型
for (event in event_types) {
  # 获取对象，例如 Aging_filt@RI_stats
  stats <- slot(Aging_filt, paste0(event, "_stats"))
  events <- slot(Aging_filt, paste0(event, "_events"))
  
  # 计算-logFDR
  stats$logFDR <- -log(stats$FDR)
  # 添加gene symbol
  stats$geneSymbol <- events$geneSymbol
  
  # 添加状态列
  stats$status <- ifelse(stats$FDR < 0.05 & stats$IncLevelDifference > 0.1, 'old',
                         ifelse(stats$FDR < 0.05 & stats$IncLevelDifference < -0.1, 'young', 'Not significant'))
  
  # 过滤非有限值
  df <- stats[is.finite(stats$logFDR), ]
  # 设置因子顺序
  df$status <- factor(df$status, levels = c("Not significant", "old", "young"))
  # 重命名列
  names(df)[names(df) == "PValue"] <- "PValue"
  names(df)[names(df) == "FDR"] <- "FDR"
  names(df)[names(df) == "IncLevelDifference"] <- "deltaPSI"
  names(df)[names(df) == "logFDR"] <- "logFDR"
  
  # 确保列顺序一致（可根据你实际的列名进行调整）
  df <- df[, c("ID", "PValue", "FDR", "deltaPSI", "logFDR", "geneSymbol", "status")]
  
  # 保存单独的vol文件
  write.csv(df, paste0("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/", event, "_vol.csv"), row.names = FALSE)
  
  # 添加一个事件类型列，用于后续合并
  df$event_type <- event
  all_dfs[[event]] <- df
}

# 合并所有数据
merged_df <- do.call(rbind, all_dfs)

# 保存合并后的文件
write.csv(merged_df, "/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/all_AS_vol.csv", row.names = FALSE)


##事件和基因数目统计
# 添加 'differ' 列
merged_df$differ <- ifelse(merged_df$status == "old", 1, 0)

merged_df<-subset(merged_df,merged_df$status %in% c("old","young"))
df<-melt(table(merged_df$status,merged_df$event_type))
df<-subset(df,df$Var1 %in% c("old","young"))

p <- ggplot(df, aes(x = Var1, y = value, fill = Var2)) +
  facet_wrap(~ Var1, ncol = 2, scales = "free_x") +
  geom_bar(stat = 'identity', position = 'dodge',
           width = 0.8, color = 'black') +
  geom_text(aes(label = value), size = 4,
            position = position_dodge(width = 0.8),
            vjust = -0.3) +
  labs(x = NULL) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(colour = 'black', size = 12, face = "bold"),  # X轴文本
    axis.text.y = element_text(colour = 'black', size = 12, face = "bold"),  # Y轴文本
    strip.text = element_text(face = "bold", size = 12),  # 分面标题
    legend.title = element_text(face = "bold", size = 12),  # 图例标题
    legend.text = element_text(face = "bold", size = 10),  # 图例文本
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)  # 图形标题
  ) +
  scale_fill_manual(values = rep(c('#D55E00', '#0072B2', '#F0E442',
                                   '#009E73', '#9270ca', '#CC79A7'), 3)) +
  labs(fill = "Type", y = "Number", x = "AS Events")

p
##gene
data<-melt(table(merged_df$status,merged_df$event_type))
data$ga<-paste0(data$Var1,data$Var2)
merged_df$ga<-paste0(merged_df$status,merged_df$event_type,merged_df$geneSymbol)
df_combined1<-merged_df[!duplicated(merged_df$ga),]#这里就保证了是要符合status AS gene三种情况相同，就可以做出没有重复得
data<-melt(table(df_combined1$status,df_combined1$event_type))
data<-subset(data,data$Var1 %in% c("old","young"))

p1<-ggplot(data,aes(x=Var1,y=value,fill=Var2))+
  #geom_bar(position = "do") +
  facet_wrap(~Var1, ncol = 2, scales = "free_x") + 
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  geom_text(aes(label=value),size=4,
            position = position_dodge(width = 0.8), 
            vjust=-0.3)+ 
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+
  theme(
    axis.text.x = element_text(colour = 'black', size = 12, face = "bold"),  # X轴文本
    axis.text.y = element_text(colour = 'black', size = 12, face = "bold"),  # Y轴文本
    strip.text = element_text(face = "bold", size = 12),  # 分面标题
    legend.title = element_text(face = "bold", size = 12),  # 图例标题
    legend.text = element_text(face = "bold", size = 10),  # 图例文本
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)  # 图形标题
  ) +
  scale_fill_manual(values = rep(c('#D55E00','#0072B2','#F0E442',
                                   '#009E73','#9270ca','#CC79A7'),3))+
  labs(fill="Type",y="Number",x="AS Genes")


(p+p1)
ggsave(filename = '/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/plot/as_gene_sum.pdf',
       (p+p1),
       width = 15,height = 6,dpi=600)
###交集
before_df<-read.csv("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/before.csv",row.names = 1)
before_df<-subset(before_df,before_df$status %in% c("old","young"))

##old
old1<-subset(before_df,before_df$status %in% c("old"))
old2<-subset(merged_df,merged_df$status %in% c("old"))
my_list <- list(old1$ID, old2$ID)
my_list[[2]] <- unique(my_list[[2]])
my_list[[1]] <- unique(my_list[[1]])
names(my_list) <- c("unrandom", "random")

# 创建euler对象，计算两个list的交集
fit <- euler(my_list, shape = "circle")

# 绘图
plot(
  fit,
  #fills = list(fill = c("#4986b5", "#34bf49"), alpha = 0.6),
  fills = list(fill = c("#F98F34", "#0C4E9B"), alpha = 0.8),
  labels = list(col = "black", font = 2, fontfamily = "serif", cex = 1.5),
  edges = list(col = "black", lwd = 3, lty = 2),  # lty = 2 表示虚线
  quantities = list(type = "counts", cex = 1)  # 只显示 counts，不显示 percent
) -> p2

p2


ggsave("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/plot/up_veen.pdf", plot = p2, width = 4, height = 3,dpi=600)

##young
young1<-subset(before_df,before_df$status %in% c("young"))
young2<-subset(merged_df,merged_df$status %in% c("young"))
my_list <- list(young1$ID, young2$ID)
my_list[[2]] <- unique(my_list[[2]])
my_list[[1]] <- unique(my_list[[1]])
names(my_list) <- c("unrandom", "random")

# 创建euler对象，计算两个list的交集
fit <- euler(my_list, shape = "circle")

# 绘图
plot(
  fit,
  #fills = list(fill = c("#4986b5", "#34bf49"), alpha = 0.6),
  fills = list(fill = c("#F98F34", "#0C4E9B"), alpha = 0.8),
  labels = list(col = "black", font = 2, fontfamily = "serif", cex = 1.5),
  edges = list(col = "black", lwd = 3, lty = 2),  # lty = 2 表示虚线
  quantities = list(type = "counts", cex = 1)  # 只显示 counts，不显示 percent
) -> p2

p2





ggsave("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/plot/down_veen.pdf", plot = p2, width = 4, height = 3,dpi=600)

intersection_genes <- intersect(my_list[[1]], my_list[[2]])

##所有AS相关基因的交集
before_df<-read.csv("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/before.csv",row.names = 1)
merged_df<-read.csv("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/all_AS_vol.csv") 
before<-subset(before_df,before_df$status %in% c("old","young"))
merged<-subset(merged_df,merged_df$status %in% c("old","young"))

my_list <- list(before$geneSymbol, merged$geneSymbol)
my_list[[2]] <- unique(my_list[[2]])
my_list[[1]] <- unique(my_list[[1]])
names(my_list) <- c("original", "subset")

# 创建euler对象，计算两个list的交集
fit <- euler(my_list, shape = "circle")

# 绘图
plot(
  fit,
  #fills = list(fill = c("#4986b5", "#34bf49"), alpha = 0.6),
  fills = list(fill = c("#F98F34", "#0C4E9B"), alpha = 0.8),
  labels = list(col = "black", font = 2, fontfamily = "serif", cex = 1.5),
  edges = list(col = "black", lwd = 3, lty = 2),  # lty = 2 表示虚线
  quantities = list(type = "counts", cex = 1)  # 只显示 counts，不显示 percent
) -> p2

p2





ggsave("/data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/plot/AS_gene_veen.pdf", plot = p2, width = 4, height = 3,dpi=600)

