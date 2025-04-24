rm(list=ls())
setwd("./file/")
peptied<-read.csv("./DASEs_Peptides_anno.csv")
####Statistics of peptide segment numbers---------------------------
data = tibble(
  name = c("Brianna","Cameron","Daisy","Ethan","Fiona","Gabriel","Hannah","Isaac"),
  year2022 = c(234, 125, 234, 346, 450, 324, 189, 540),
  year2023 = c(256, 110, 321, 240, 600, 200, 210, 470)
)
mpg %>% 
  mutate(manufacturer = str_to_title(manufacturer)) %>% 
  group_by(manufacturer) %>% 
  summarise(n = n()) %>% 
  slice_max(n, n=10) %>% 
  arrange(n) %>% 
  mutate(id = 1:n())-> data
data %>% 
  ggplot(aes(x=id, y=n, fill = manufacturer)) +
  geom_bar(stat = 'identity', width = 0.7) +
  
  geom_text(aes(x=id, y=0, label = manufacturer),
            hjust=1.1, size=3.5, color="black") +
  scale_fill_manual(values = pal10, guide="none") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 50),
                     position = 'right') +
  scale_x_continuous(expand = c(0, 0), 
                     limits=c(-1, 10.5)) +
  coord_polar(theta = 'y') +
  labs(title = 'CircularBar Chart',
       subtitle = "Number of models per manufacturer, 1999-2008")+
  theme_minimal() +
  theme(plot.background = element_rect(fill='white', color='white'),
        axis.title = element_blank(),
        axis.text = element_blank())

###Peptide density-----------
different_df <- peptied %>%
  filter(CDS_same == "Different") %>%
  mutate(peptide_count = sapply(strsplit(dif_cut_seq, ","), length))  # 计算每行的肽段数目
different_df<-different_df[-11,]

max_peptides <- max(different_df$peptide_count, na.rm = TRUE)#横坐标最大5208

sorted_peptides <- sort(unique(different_df$peptide_count), decreasing = TRUE)
second_max_peptides <- sorted_peptides[2]  # 第二大的数值4509
third_max_peptides <- sorted_peptides[3]#三大4200

# 
p<-ggplot(different_df, aes(x = peptide_count, fill = AS_type, color = AS_type)) +
  geom_density(alpha = 0.5) +
  labs(x = "Number of Peptides",
       y = "Density") +
  scale_fill_manual(values = c("SE" = "#4986b5", "MXE" = "#65bd54", "A3SS" = "#ff6347", "A5SS" = "#9370db", "RI" = "#ffa500")) +
  scale_color_manual(values = c("SE" = "#4986b5", "MXE" = "#65bd54", "A3SS" = "#ff6347", "A5SS" = "#9370db", "RI" = "#ffa500")) +
  scale_x_continuous(limits = c(0, max_peptides)) +  # 设置 x 轴的范围
  theme_minimal() +
  theme(legend.title = element_blank())  +
  theme_minimal() +  # 使用更简洁的主题
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )+ # 设置图例标题
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # 移除灰色背景
    axis.text = element_text(face = "bold"), # 加粗坐标轴数字
    axis.title = element_text(face = "bold"),
    panel.border = element_blank(), # 移除黑色边框线
    axis.line = element_line(size = 1, colour = 'black'),
    axis.ticks = element_line(size = 1, colour = 'black'),
    legend.position = "right" # 根据需要调整图例位置
  ) +
  theme(
    axis.text = element_text(size = 18, colour = 'black'),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  ) +
  scale_x_continuous(limits = c(0, 5208), breaks = seq(0, 5208, by = 1302))
p+ylim(0,0.03)
ggsave(plot=p+ylim(0,0.03),"./plot/diff_pep.pdf",width = 10,height = 6)


p<-ggplot(different_df, aes(x = peptide_count, fill = AS_type, color = AS_type)) +
  geom_density(alpha = 0.5) +
  labs(x = "Number of Peptides",
       y = "Density") +
  scale_fill_manual(values = c("SE" = "#4986b5", "MXE" = "#65bd54", "A3SS" = "#ff6347", "A5SS" = "#9370db", "RI" = "#ffa500")) +
  scale_color_manual(values = c("SE" = "#4986b5", "MXE" = "#65bd54", "A3SS" = "#ff6347", "A5SS" = "#9370db", "RI" = "#ffa500")) +
  scale_x_continuous(limits = c(0, 1000)) +  # 设置 x 轴的范围
  theme_minimal() +
  theme(legend.title = element_blank())  +
  theme_minimal() +  # 使用更简洁的主题
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )+ # 设置图例标题
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # 移除灰色背景
    axis.text = element_text(face = "bold"), # 加粗坐标轴数字
    axis.title = element_text(face = "bold"),
    panel.border = element_blank(), # 移除黑色边框线
    axis.line = element_line(size = 1, colour = 'black'),
    axis.ticks = element_line(size = 1, colour = 'black'),
    legend.position = "right" # 根据需要调整图例位置
  ) +
  theme(
    axis.text = element_text(size = 18, colour = 'black'),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  )
ggsave(plot=p+ylim(0,0.04),"./plot/1000_diff_pep.pdf",width = 10,height = 6)


different_df <- different_df %>%
  mutate(peptide_count = ifelse(peptide_count > 1000, 1000, peptide_count))

# 计算最大值用于设置x轴范围
max_peptides <- max(different_df$peptide_count)

# 绘制密度图
p <- ggplot(different_df, aes(x = peptide_count, fill = AS_type, color = AS_type)) +
  geom_density(alpha = 0.5) +
  labs(x = "Number of Peptides",
       y = "Density") +
  scale_fill_manual(values = c("SE" = "#4986b5", "MXE" = "#65bd54", "A3SS" = "#ff6347", "A5SS" = "#9370db", "RI" = "#ffa500")) +
  scale_color_manual(values = c("SE" = "#4986b5", "MXE" = "#65bd54", "A3SS" = "#ff6347", "A5SS" = "#9370db", "RI" = "#ffa500")) +
  scale_x_continuous(limits = c(0, max_peptides)) +  # 设置 x 轴的范围
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # 移除灰色背景
    axis.text = element_text(face = "bold", size = 18), # 加粗坐标轴数字
    axis.title = element_text(face = "bold", size = 20),
    panel.border = element_blank(), # 移除黑色边框线
    axis.line = element_line(size = 1, colour = 'black'),
    axis.ticks = element_line(size = 1, colour = 'black'),
    legend.position = "right",
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  )
p
ggsave(plot=p+ylim(0,0.04),"./plot/v2_1000_diff_pep.pdf",width = 10,height = 6)



write.csv(different_df,"./file/diff_pep_df.csv")


###Re-count the data of the two libraries----------------------------
#################young
fasta <- readLines("./combined_peptides.fa")
transcript_ids <- fasta[grep("^>", fasta)]
transcript_ids <- gsub(">", "", transcript_ids) 
peptides <- fasta[grep(",", fasta)]
peptide_counts <- sapply(peptides, function(line) length(strsplit(line, ",")[[1]]))
total_peptide_count <- sum(peptide_counts)

###########################mid
fasta <- readLines("./combined_peptides.fa")
transcript_ids <- fasta[grep("^>", fasta)]
transcript_ids <- gsub(">", "", transcript_ids)  
peptides <- fasta[grep(",", fasta)]
peptide_counts <- sapply(peptides, function(line) length(strsplit(line, ",")[[1]]))
total_peptide_count <- sum(peptide_counts)

####减掉青年库剩下的------------------------------
A3SS<-read.csv("./A3SS3.csv",row.names = 1)
A3SS_df <- A3SS %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )
A5SS<-read.csv("./noyoung/A5SS3.csv",row.names = 1)
A5SS_df <- A5SS %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )
MXE<-read.csv("./noyoung/MXE3.csv",row.names = 1)
MXE_df <- MXE %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )

MXE1<-read.csv("./noyoung/MXE2.csv",row.names = 1)
MXE_df <- MXE1 %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )
RI<-read.csv("./noyoung/RI3.csv",row.names = 1)
RI_df <- RI %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )
RI1<-read.csv("./v2/noyoung/RI2.csv",row.names = 1)
RI_df <- RI1 %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )
SE<-read.csv("./noyoung/SE3.csv",row.names = 1)
SE_df <- SE %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )
SE1<-read.csv("./noyoung/SE2.csv",row.names = 1)
SE_df <- SE1 %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )
##########Pre-data organization of the middle-aged peptide library----------------------------------------
write.csv(A3SS,"./noyoung/csv/A3SS_df.csv")
write.csv(A5SS,"./noyoung/csv/A5SS_df.csv")
MXE2<-rbind(MXE,MXE1)
MXE2<-MXE2[which(MXE2$dif_cut_seq!=""),]
write.csv(MXE2,"./noyoung/csv/MXE_df.csv")
RI2<-rbind(RI,RI1)
RI2<-RI2[which(RI2$dif_cut_seq!=""),]
write.csv(RI2,"./noyoung/csv/RI_df.csv")
SE2<-rbind(SE,SE1)
SE2<-SE2[which(SE2$dif_cut_seq!=""),]
write.csv(RI2,"./noyoung/csv/SE_df.csv")
merge_df<-rbind(A3SS,A5SS,MXE2,RI2,SE2)
write.csv(merge_df,"./noyoung/csv/merge_df.csv")
##用merge_df跑uniq.r代码
##整理
A3SS_unmatched<-read.csv("./csv/A3SS3.csv")
A3SS_df <- A3SS_unmatched %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )

A5SS_unmatched<-read.csv("./csv/A5SS3.csv")
A5SS_df <- A5SS_unmatched %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )
MXE_matched<-read.csv("./csv/MXE2.csv")
MXE_unmatched<-read.csv("./v2/csv/MXE3.csv")

MXE_df <- MXE_unmatched %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )

RI_unmatched<-read.csv("./v2/csv/RI3.csv")
RI_df <- RI_unmatched %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )
SE_matched<-read.csv("./csv/SE2.csv")
SE_unmatched<-read.csv("./csv/SE3.csv")
SE_df <- SE_unmatched %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(dif_cut_seq, ","), length))
  )

all_noyoung_df<-rbind(A3SS_unmatched,A5SS_unmatched,MXE_matched,MXE_unmatched,RI_unmatched,SE_unmatched)
all_noyoung_df1<-all_noyoung_df[which(all_noyoung_df$dif_cut_seq!=""),]
write.csv(all_noyoung_df1,"./csv/noyoung_df.csv")

###The remaining peptide segments are reanalyzed-------------------------
rm(list=ls())
fasta<-read.table("./combined_peptides.fa",sep = "\t")
fasta1<-read.table("/./combined_peptides.fa",sep = "\t")
data<-read.csv("./noyoung_df.csv",row.names = 1)


a1<-fasta[seq(2,261610,2),]
a2<-unique(unlist(strsplit(a1,",")))
data$seq1<-"0"
for(i in unique(data$AS_type)){
  num1<-grep(i,data$AS_type)
  data1<-data[num1,]
  for(j in 1:length(num1)){
    a3=unlist(strsplit(data1$dif_cut_seq[j],","))
    data$seq1[num1[j]]<-paste0(a3[! a3 %in% a2],collapse = ",")
  }
}
##先减掉young库的
noy<-read.csv("./noyoung_pep.csv",row.names = 1)
noy1<-noy[which(noy$seq1!=""),]
summary_df <- noy1 %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(seq1, ","), length))
  )
write.csv(noy1,"./noy1.csv")
##减掉mid库所剩的
noym<-read.csv("./noym_pep.csv",row.names = 1)
noym1<-noym[which(noym$seq2!=""),]
summary_df <- noym1 %>%
  group_by(AS_type) %>%
  summarise(
    Tranid_Count = n_distinct(tranid),
    Peptide_Count = sum(sapply(strsplit(seq2, ","), length))
  )
write.csv(noym1,"./noym1.csv")

###Process the output file of netMHC------------------------------------
##cell member-----------------------------
net<-read.csv("./netMHC_out.csv",header = F)
b<-c()
d<-c()
list_net<-list()
for (i in grep("HLA-",net$V4)) {
  num<-which(net$V9[(i+2):(i+2291)]<=0.5)
  net$V9[i]<-length(num)
  b<-c(b,i,i+1,num+1+i)
  data<-net[c(i+1,num+1+i),]
  colnames(data)<-data[1,]
  list_net[net$V4[i]]<-list(data<-data[-1,])
}
net1<-net[b,]
write.csv(net1,"./CM_pep.csv")
##肽段统计
library(dplyr)
library(purrr)
#统计所有和HLA分子的唯一的肽段
unique_peptides <- list_net %>% 
  map(~ .x$Peptide) %>%     
  unlist() %>%              
  unique()                
print(unique_peptides)
length(unique_peptides)#421
save(list_net,file = "./CM_list.RData")

###将肽段mapping回去
result_df$seq2_list <- strsplit(as.character(result_df$seq2), split = ",", fixed = TRUE)


# 定义查找函数
find_peptide_positions <- function(peptide, dataframe) {
  # 使用sapply检查每行的seq2_list是否包含特定的肽段
  positions <- sapply(dataframe$seq2_list, function(peptide_list) {
    peptide %in% peptide_list
  })
  return(which(positions))
}

empty_lists <- replicate(length(unique_peptides), list(), simplify = FALSE)

# 创建数据框
results <- data.frame(Peptide = unique_peptides, Positions = I(empty_lists))

# 填充results数据框中的Positions列
for (i in seq_along(unique_peptides)) {
  results$Positions[[i]] <- find_peptide_positions(unique_peptides[i], result_df)
}
write.csv(results,"./MHC/CM_pep_pos.csv")

expanded_results <- results %>%
  unnest(c(Positions))

grouped_results <- expanded_results %>%
  group_by(Positions) %>%
  summarise(Peptides = toString(Peptide),
            NumPeptides = n()) %>%
  arrange(Positions)

result_df <- result_df %>%
  mutate(RowName = row_number())  

grouped_results$Positions <- as.numeric(grouped_results$Positions)
final_results <- left_join(grouped_results, result_df, by = c("Positions" = "RowName"))
final_results <- final_results[,c(1,2,3,4,7,8)]
write.csv(final_results,"./MHC/CM_pep_info.csv")

##画图展示
summary_df <- final_results %>%
  group_by(AS_type) %>%
  summarise(TotalNumPeptides = sum(NumPeptides))

p<-ggplot(summary_df, aes(x = AS_type, y = TotalNumPeptides, fill = AS_type)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +  # 控制条柱的宽度
  geom_text(aes(label = TotalNumPeptides), vjust = -0.3, position = position_dodge(width = 0.7), 
            color = "black", size = 4,face = "bold") +  # 添加数值标签
  labs(
    x = "AS Type",
    y = "Total Number of Peptides") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")  +
  theme_minimal() +
  theme(legend.title = element_blank())  +
  theme_minimal() +  # 使用更简洁的主题
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )+ # 设置图例标题
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # 移除灰色背景
    axis.text = element_text(face = "bold"), # 加粗坐标轴数字
    axis.title = element_text(face = "bold"),
    panel.border = element_blank(), # 移除黑色边框线
    axis.line = element_line(size = 1, colour = 'black'),
    axis.ticks = element_line(size = 1, colour = 'black'),
    legend.position = "right" # 根据需要调整图例位置
  ) +
  theme(
    axis.text = element_text(size = 18, colour = 'black'),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  )

p
ggsave(plot=p,"./plot/CM_pep.pdf",width = 10,height = 6)






##no cell member的-----------------------------------------------------

net2<-read.table("./netMHC_out.txt",header = F,sep = "\t")
b<-c()
d<-c()
list_net1<-list()
for (i in grep("HLA-",net2$V4)) {
  num<-which(net2$V9[(i+2):(i+68661)]<=0.5)
  net2$V9[i]<-length(num)
  b<-c(b,i,i+1,num+1+i)
  data<-net2[c(i+1,num+1+i),]
  colnames(data)<-data[1,]
  list_net1[net2$V4[i]]<-list(data<-data[-1,])
}
net3<-net2[b,]
save(list_net1,file = "./noCM_list.RData")
write.csv(net3,"./MHC/noCM_pep.csv")
unique_pep <- list_net1 %>% 
  map(~ .x$Peptide) %>%     
  unlist() %>%              
  unique()                
print(unique_pep)
length(unique_pep)#13600


###将肽段mapping回去
noym1$seq2_list <- strsplit(as.character(noym1$seq2), split = ",", fixed = TRUE)


# 定义查找函数
find_peptide_positions <- function(peptide, dataframe) {
  # 使用sapply检查每行的seq2_list是否包含特定的肽段
  positions <- sapply(dataframe$seq2_list, function(peptide_list) {
    peptide %in% peptide_list
  })
  return(which(positions))
}

empty_lists <- replicate(length(unique_pep), list(), simplify = FALSE)

# 创建数据框
results <- data.frame(Peptide = unique_pep, Positions = I(empty_lists))

# 填充results数据框中的Positions列
for (i in seq_along(unique_pep)) {
  results$Positions[[i]] <- find_peptide_positions(unique_pep[i], noym1)
}
write.csv(results,"./noCM_pep_pos.csv")

expanded_results <- results %>%
  unnest(c(Positions))

grouped_results <- expanded_results %>%
  group_by(Positions) %>%
  summarise(Peptides = toString(Peptide),
            NumPeptides = n()) %>%
  arrange(Positions)

noym1 <- noym1 %>%
  mutate(RowName = row_number())  

grouped_results$Positions <- as.numeric(grouped_results$Positions)
final_results <- left_join(grouped_results, noym1, by = c("Positions" = "RowName"))
final_results <- final_results[,c(1,2,3,7,9)]
write.csv(final_results,"./noCM_pep_info.csv")

##画图展示
summary_df <- final_results %>%
  group_by(AS_type) %>%
  summarise(TotalNumPeptides = sum(NumPeptides))

p<-ggplot(summary_df, aes(x = AS_type, y = TotalNumPeptides, fill = AS_type)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +  # 控制条柱的宽度
  geom_text(aes(label = TotalNumPeptides), vjust = -0.3, position = position_dodge(width = 0.7), 
            color = "black", size = 4,face = "bold") +  # 添加数值标签
  labs(
    x = "AS Type",
    y = "Total Number of Peptides") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")  +
  theme_minimal() +
  theme(legend.title = element_blank())  +
  theme_minimal() +  # 使用更简洁的主题
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )+ # 设置图例标题
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # 移除灰色背景
    axis.text = element_text(face = "bold"), # 加粗坐标轴数字
    axis.title = element_text(face = "bold"),
    panel.border = element_blank(), # 移除黑色边框线
    axis.line = element_line(size = 1, colour = 'black'),
    axis.ticks = element_line(size = 1, colour = 'black'),
    legend.position = "right" # 根据需要调整图例位置
  ) +
  theme(
    axis.text = element_text(size = 18, colour = 'black'),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  )

p
ggsave(plot=p,"./noCM_pep.pdf",width = 10,height = 6)

###Organize the diagrams of HLA-----------------------
###cellmember-----------------------
rm(list=ls())
load(file="./CM_list.RData")
##HLA和肽段数目的关系
HLA<-read.table("./MHC/ss3.txt",sep = "\t")

hla_counts <- HLA %>%
  group_by(V1) %>%
  summarise(Count = n())
total_count <- sum(hla_counts$Count)
hla_counts <- hla_counts %>%
  mutate(Percentage = Count / total_count * 100)

pep_num <- sapply(list_net, nrow)
pep_num_df <- data.frame(V1 = names(pep_num), pep_num = pep_num)
hla_counts <- left_join(hla_counts, pep_num_df, by = "V1")
hla_counts <- hla_counts %>%
  mutate(HLA_type = str_extract(V1, "^HLA-[A-Z]+"))
write.csv(hla_counts,"./CM_pep_HLA.csv")
##画图
p <- ggplot(data = hla_counts, aes(x = Percentage, y = pep_num, color = HLA_type)) +
  geom_point(size = 4, alpha = 0.9) +  # 点的大小和透明度
  
  # 只为满足条件的数据添加标签
  geom_text(data = subset(hla_counts, Percentage > 2 | pep_num > 38), 
            aes(label = V1), vjust = 1.6, color = "black", size = 5,face = "bold") +
  
  # 自定义颜色
  scale_color_manual(values = c("HLA-A" = "#e86d67", "HLA-B" = "#3ab05b", "HLA-C" = "#6489c1")) +
  
  # 设置图形的标题和轴标签
  labs(
    x = "HLA Frequency (%)",
    y = "Peptide Binding Number",
    color = "HLA Type") +
  
  # 使用最小化主题，并调整图例位置
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.background = element_blank(),  # 移除灰色背景
    panel.border = element_rect(colour = "black", fill = NA, size = 2),  # 加粗四周边框
    legend.background = element_rect(fill = "white", colour = "black"),  # 图例背景调整
    axis.text = element_text(color = "black", size = 12),  # 调整坐标轴文本
    axis.title = element_text(size = 14)  # 调整坐标轴标题
  )+
  theme(
    axis.text = element_text(size = 18, colour = 'black'),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  )+
  scale_y_continuous(limits = c(0, 55), breaks = seq(0, 55 ,by = 5))
p
ggsave(plot=p,"./CM_pep_HLA.pdf",width = 15,height = 12)

##top10出现次数最多的肽段与HLA类型的关系
peptide_counts <- lapply(names(list_net), function(hla_name) {
  df <- list_net[[hla_name]]
  hla_type <- ifelse(grepl("HLA-A", hla_name), "A",
                     ifelse(grepl("HLA-B", hla_name), "B", "C"))
  data.frame(Peptide = df$Peptide, HLA_type = hla_type)
}) %>%
  bind_rows() %>%
  group_by(Peptide, HLA_type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  ungroup()

# 汇总总数，并选出Top 20肽段
top_peptides <- peptide_counts %>%
  group_by(Peptide) %>%
  summarise(TotalCount = sum(Count)) %>%
  arrange(desc(TotalCount)) %>%
  slice_head(n = 20)

# 合并Top 20肽段详细数据
top_peptides_detailed <- top_peptides %>%
  left_join(peptide_counts, by = "Peptide") %>%
  pivot_wider(names_from = HLA_type, values_from = Count, values_fill = list(Count = 0))


long_data <- top_peptides_detailed %>%
  pivot_longer(cols = c("A", "B", "C"), names_to = "HLA_type", values_to = "Num")

long_data <- long_data %>%
  arrange(desc(TotalCount)) %>%
  mutate(Peptide = factor(Peptide, levels = unique(Peptide)))

long_data <- long_data %>%
  mutate(HLA_type = paste("HLA-", HLA_type, sep = ""))  

write.csv(long_data,"./CM_pep_hla_count.csv")

p<-ggplot(long_data, aes(x = Peptide, y = Num, fill = HLA_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("HLA-A" = "#e86d67", "HLA-B" = "#3ab05b", "HLA-C" = "#6489c1")) +
  labs( x = "Peptide", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    legend.position = "right",
    panel.background = element_blank(),  # 移除灰色背景
    panel.border = element_rect(colour = "black", fill = NA, size = 2),  # 加粗四周边框
    legend.background = element_rect(fill = "white", colour = "black"),  # 图例背景调整
    axis.text = element_text(color = "black", size = 12),  # 调整坐标轴文本
    axis.title = element_text(size = 14)  # 调整坐标轴标题
  )+
  theme(
    axis.text = element_text(size = 18, colour = 'black'),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  ) +
  scale_y_continuous(limits = c(0, 55), breaks = seq(0, 55 ,by = 5))
p
ggsave(plot=p,"./v2_CM_pep_HLA.pdf",width = 15,height = 12)


###cellmember not screened------------------------------
rm(list=ls())
load(file="./noCM_list.RData")
##HLA和肽段数目的关系
HLA<-read.table("./MHC/ss3.txt",sep = "\t")

hla_counts <- HLA %>%
  group_by(V1) %>%
  summarise(Count = n())
total_count <- sum(hla_counts$Count)
hla_counts <- hla_counts %>%
  mutate(Percentage = Count / total_count * 100)

pep_num <- sapply(list_net1, nrow)
pep_num_df <- data.frame(V1 = names(pep_num), pep_num = pep_num)
hla_counts <- left_join(hla_counts, pep_num_df, by = "V1")
hla_counts <- hla_counts %>%
  mutate(HLA_type = str_extract(V1, "^HLA-[A-Z]+"))
write.csv(hla_counts,"./noCM_pep_HLA.csv")
##画图
p <- ggplot(data = hla_counts, aes(x = Percentage, y = pep_num, color = HLA_type)) +
  geom_point(size = 4, alpha = 0.9) +  # 点的大小和透明度
  
  # 只为满足条件的数据添加标签
  geom_text(data = subset(hla_counts, Percentage > 3| pep_num > 1000), 
            aes(label = V1), vjust = 1.6, color = "black", size = 5,face = "bold") +
  
  # 自定义颜色
  scale_color_manual(values = c("HLA-A" = "#e86d67", "HLA-B" = "#3ab05b", "HLA-C" = "#6489c1")) +
  
  # 设置图形的标题和轴标签
  labs(
    x = "HLA Frequency (%)",
    y = "Peptide Binding Number",
    color = "HLA Type") +
  
  # 使用最小化主题，并调整图例位置
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.background = element_blank(),  # 移除灰色背景
    panel.border = element_rect(colour = "black", fill = NA, size = 2),  # 加粗四周边框
    legend.background = element_rect(fill = "white", colour = "black"),  # 图例背景调整
    axis.text = element_text(color = "black", size = 12),  # 调整坐标轴文本
    axis.title = element_text(size = 14)  # 调整坐标轴标题
  )+
  theme(
    axis.text = element_text(size = 18, colour = 'black'),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  )+
  scale_y_continuous(limits = c(0, 1300), breaks = seq(0, 1300 ,by = 260))
p
ggsave(plot=p,"./noCM_pep_HLA.pdf",width = 16,height = 15)

##The relationship between the top10 most frequently occurring peptide segments and HLA types
peptide_counts <- lapply(names(list_net1), function(hla_name) {
  df <- list_net1[[hla_name]]
  hla_type <- ifelse(grepl("HLA-A", hla_name), "A",
                     ifelse(grepl("HLA-B", hla_name), "B", "C"))
  data.frame(Peptide = df$Peptide, HLA_type = hla_type)
}) %>%
  bind_rows() %>%
  group_by(Peptide, HLA_type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  ungroup()

# 汇总总数，并选出Top 20肽段
top_peptides <- peptide_counts %>%
  group_by(Peptide) %>%
  summarise(TotalCount = sum(Count)) %>%
  arrange(desc(TotalCount)) #%>%
#slice_head(n = 20)

# 合并Top 20肽段详细数据
top_peptides_detailed <- top_peptides %>%
  left_join(peptide_counts, by = "Peptide") %>%
  pivot_wider(names_from = HLA_type, values_from = Count, values_fill = list(Count = 0))


long_data <- top_peptides_detailed %>%
  pivot_longer(cols = c("A", "B", "C"), names_to = "HLA_type", values_to = "Num")

long_data <- long_data %>%
  arrange(desc(TotalCount)) %>%
  mutate(Peptide = factor(Peptide, levels = unique(Peptide)))

long_data <- long_data %>%
  mutate(HLA_type = paste("HLA-", HLA_type, sep = ""))  

write.csv(long_data,"./no_CM_pep_hla_count.csv")

p<-ggplot(long_data, aes(x = Peptide, y = Num, fill = HLA_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("HLA-A" = "#e86d67", "HLA-B" = "#3ab05b", "HLA-C" = "#6489c1")) +
  labs( x = "Peptide", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    legend.position = "right",
    panel.background = element_blank(),  # 移除灰色背景
    panel.border = element_rect(colour = "black", fill = NA, size = 2),  # 加粗四周边框
    legend.background = element_rect(fill = "white", colour = "black"),  # 图例背景调整
    axis.text = element_text(color = "black", size = 12),  # 调整坐标轴文本
    axis.title = element_text(size = 14)  # 调整坐标轴标题
  )+
  theme(
    axis.text = element_text(size = 18, colour = 'black'),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  ) 
#scale_y_continuous(limits = c(0, 55), breaks = seq(0, 55 ,by = 5))
p
ggsave(plot=p,"./plot/v2_no_CM_pep_HLA.pdf",width = 15,height = 12)

###See if the final gene is related to aging---------------------
rm(list=ls())
cell<-read.csv("./CM_pep_info.csv")
nocell<-read.csv("./noCM_pep_info.csv")
aging_gene<-read.csv("./genage_human.csv")
result_df <- inner_join(nocell, aging_gene, by = c("geneSymbol" = "symbol"))
result_df <- inner_join(cell, aging_gene, by = c("geneSymbol" = "symbol"))
##Count the first 28 peptides presented by all 80 elderly people-------------------------------
selected_peptides <- peptide_freq$Peptides[1:28]

# 创建一个空列表来存储交集行
intersected_rows <- list()

# 遍历图2的数据框
for (i in 1:nrow(df)) {
  # 按逗号分隔图2的Peptides列
  peptides_in_row <- unlist(strsplit(df$Peptides[i], ", "))
  
  # 找到与图1前28个肽段的交集
  intersection <- intersect(peptides_in_row, selected_peptides)
  
  # 如果存在交集，将该行添加到结果中
  if (length(intersection) > 0) {
    intersected_rows[[length(intersected_rows) + 1]] <- df[i, ]
  }
}

# 将列表转换为数据框
result_df <- do.call(rbind, intersected_rows)
write.csv(result_df,"./protein/28_info.csv")

##看这些基因与aging、sasp交集
SASP<- read.csv("./SASP/SASP_SenMayo.csv",sep = "\t")
SASP<-subset(SASP,SASP$Annotation=="SASP")
common_genes <- intersect(result_df$geneSymbol, SASP$Symbol)#无

aging_gene<- read.csv("./file/genage_human.csv")
common_genes <- intersect(result_df$geneSymbol, aging_gene$symbol)


###Statistical chart of all splicing events and gene numbers-----------------------
rm(list=ls())
RI<-read.csv("./file/RI_vol.csv")
SE<-read.csv("./file/SE_vol.csv")
MXE<-read.csv("./file/MXE_vol.csv")
A3SS<-read.csv("./file/A3SS_vol.csv")
A5SS<-read.csv("./file/A5SS_vol.csv")
# 为每个数据框添加 'AS_type' 列
SE$AS_type <- "SE"
RI$AS_type <- "RI"
MXE$AS_type <- "MXE"
A3SS$AS_type <- "A3SS"
A5SS$AS_type <- "A5SS"

# 合并这五个数据框
df_combined <- bind_rows(SE, RI, MXE, A3SS, A5SS)

# 添加 'differ' 列
df_combined$differ <- ifelse(df_combined$status == "old", 1, 0)

df_combined<-subset(df_combined,df_combined$status %in% c("old","young"))
df<-melt(table(df_combined$status,df_combined$AS_type))
# df1<-df[order(df$Var1,df$value),]
# 
# df1$group<-paste0(df1$Var1,1:5)
# df1$type_order=factor(rev(as.integer(1:dim(df1)[1])),labels =rev(df1$group))

#基因
data<-melt(table(df_combined$status,df_combined$AS_type))
data$ga<-paste0(data$Var1,data$Var2)
df_combined$ga<-paste0(df_combined$status,df_combined$AS_type,df_combined$geneSymbol)
df_combined1<-df_combined[!duplicated(df_combined$ga),]#这里就保证了是要符合status AS gene三种情况相同，就可以做出没有重复得
data<-melt(table(df_combined1$status,df_combined1$AS_type))

##并起来
df<-df[,-1]
data<-data[,-1]

df_summarized <- df %>%
  group_by(Var2) %>%
  summarize(sum_value = sum(value))

data_summarized <- data %>%
  group_by(Var2) %>%
  summarize(sum_value = sum(value))



df_summarized$type <- "Events"
data_summarized$type <- "Genes"
df_combined <- bind_rows(df_summarized, data_summarized)

# 计算每个类型的总数并添加到数据框中
df_combined <- df_combined %>%
  group_by(type) %>%
  mutate(total_value = sum(sum_value))

# 绘制堆栈图
p <- ggplot(df_combined, aes(x = type, y = sum_value, fill = Var2)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.8, color = 'black') +
  geom_text(aes(y = total_value, label = total_value), size = 6, vjust = -0.3,colour = 'black') +
  facet_wrap(~ type, ncol = 2, scales = "free_x") +
  labs(x = NULL, y = "Number", fill = "Type", face = "bold") +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(colour = 'black', size = 16, face = "bold"),
    axis.text.y = element_text(colour = 'black', size = 16, face = "bold"),
    strip.text = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 18),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  ) +
  scale_fill_manual(values = c('#D55E00', '#0072B2', '#F0E442', '#009E73', '#9270ca', '#CC79A7')) +
  labs(title = "Total AS Events and Genes", y = "Number")

print(p)

ggsave(filename = './plot/total_as_gene_.pdf',
       p,
       width = 10,height = 8,dpi=600)

####Proteomic data were re-collected and processed---------------------
##PXD050061----------
rm(list=ls())
df<-read.table("./protein/PXD050061.xls",sep = "\t")
colnames(df)<-df[1,]
df<-df[-1,]
# 筛选出以“O”开头的样本的量化值列
columns <- colnames(df)
o_sample_columns <- grep("^\\[\\d+\\] O.*\\.PEP\\.Quantity$", columns, value = TRUE)

##提old列
required_columns <- c("PG.Genes", "PEP.StrippedSequence",
                      "[1] O4.raw.PEP.Quantity", "[2] O5.raw.PEP.Quantity", "[3] O6.raw.PEP.Quantity",
                      "[7] O1.raw.PEP.Quantity", "[8] O2.raw.PEP.Quantity", "[9] O3.raw.PEP.Quantity")


existing_columns <- required_columns[required_columns %in% colnames(df)]
selected_data <- df[, existing_columns]

##转数据类型
quantity_columns <- grep("Quantity$", names(selected_data), value = TRUE)

# 将这些数量列从字符型转换为数值型
selected_data[quantity_columns] <- lapply(selected_data[quantity_columns], function(x) as.numeric(as.character(x)))

# 再次查看结构确认转换
str(selected_data[quantity_columns])

##删掉全为NaN的行
selected_data[3:8] <- lapply(selected_data[3:8], function(x) {
  x <- as.character(x)   # 确保数据是字符类型
  x[x == "NaN"] <- NA    # 替换 "NaN" 为 NA
  as.numeric(x)          # 将字符转换为数值
})
selected_data <- selected_data[rowSums(is.na(selected_data[3:8])) != 6, ]
write.csv(selected_data,"./protein/PXD050061_pep.csv")

##PXD034030----------------
#年纪60-90
rm(list=ls())
data <- read_excel("./protein/PXD034030.xlsx")
colnames(data)<-data[2,]
data<-data[-c(1:2),]
data$gene <- gsub(".*GN=([^ ]+).*", "\\1", data$Description)
data<-data[,c(69,10,15)]
result <- data %>%
  group_by(gene) %>%
  distinct(Sequence, .keep_all = TRUE)
write.csv(result,"./protein/PXD034030_pep.csv")

##File processing for re-building the database------------
df_pep<-read.csv("./protein/df_pep.txt")
a<-read.csv("./protein/PXD050061_pep.csv")

generate_fasta <- function(gene, sequence) {
  fasta_header <- paste(">", gene, sep="")
  fasta_sequence <- sequence
  return(c(fasta_header, fasta_sequence))
}

# 应用这个函数到数据框的每一行
fasta_entries <- apply(a, 1, function(row) generate_fasta(row['PG.Genes'], row['PEP.StrippedSequence']))

# 将生成的列表转换为单个字符向量
fasta_text <- unlist(fasta_entries)
writeLines(fasta_text, "./protein/v2/PXD050061.txt")

b<-read.csv("./protein/PXD034030_pep.csv")
fasta_entries <- apply(b, 1, function(row) generate_fasta(row['gene'], row['Sequence']))

# 将生成的列表转换为单个字符向量
fasta_text <- unlist(fasta_entries)
writeLines(fasta_text, "./protein/v2/PXD034030.txt")

###统计一番
blast_data <- read.table("./v2_old_you/protein/v2/name.blast", header = F, sep = "\t")
sorted_df <- blast_data[order(blast_data$V3, decreasing = TRUE), ]
filtered_df <- sorted_df[sorted_df$V3 == 100, ]
filtered_df$pep_gene <- sapply(strsplit(as.character(filtered_df$V1), "_"), function(x) tail(x, 1))
# ATP5B    C9orf3     CAPN2     IGSF3     KDM1B     KDM6B      LCP1   PIKFYVE      PLEK 
# 2         1         3         1         1         1         9         8        27 
# SELL   TMEM116 TMPRSS11B    UBE2D2 
# 1         4         1         1 
#UBE2D2有研究和衰老相关吧
write.csv(filtered_df,"./v2_old_you/protein/v2/ture.csv")

##统计被80个老人都呈递的那20个肽段由什么HLA分子呈递最多---------------------
top_peptides <- peptide_freq$Peptides[1:20]

# 在df2中筛选这些peptides
filtered_df <- long_data %>% 
  dplyr::filter(Peptide %in% top_peptides)

p<-ggplot(filtered_df, aes(x = Peptide, y = Num, fill = HLA_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("HLA-A" = "#e86d67", "HLA-B" = "#3ab05b", "HLA-C" = "#6489c1")) +
  labs( x = "Peptide", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    legend.position = "right",
    panel.background = element_blank(),  # 移除灰色背景
    panel.border = element_rect(colour = "black", fill = NA, size = 2),  # 加粗四周边框
    legend.background = element_rect(fill = "white", colour = "black"),  # 图例背景调整
    axis.text = element_text(color = "black", size = 12),  # 调整坐标轴文本
    axis.title = element_text(size = 14)  # 调整坐标轴标题
  )+
  theme(
    axis.text = element_text(size = 18, colour = 'black'),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20, colour = 'black', face = "bold"), # 改变图例文字的字体大小
    legend.title = element_text(size = 20, colour = 'black', face = "bold")
  ) 
#scale_y_continuous(limits = c(0, 55), breaks = seq(0, 55 ,by = 5))
p
ggsave(plot=p,"./plot/80_no_CM_pep_HLA.pdf",width = 15,height = 12)


##The table of peptide segments matched by the proteome was statistically analyzed-------------------------
genes_to_extract <- c("ATP5B", "C9orf3", "CAPN2", "IGSF3", "KDM1B", "KDM6B", 
                      "LCP1", "PIKFYVE", "PLEK", "SELL", "TMEM116", "TMPRSS11B", "UBE2D2")


filtered_df <- final_results %>%
  dplyr::filter(geneSymbol %in% genes_to_extract)


filtered_df <- filtered_df %>%
  mutate(NumPeptides = str_count(Peptides, ",") + 1)


write.csv(filtered_df,"/data3/lishuhan/Article/附表/protein.csv")

##Redraw Figure 6 hours----------------------------------------------
df_long <- HLA_old1 %>%
  gather(key = "HLA_type", value = "Count", HLA_A_Count:HLA_C_Count)
df_long$name<-factor(df_long$sample,levels = paste0("old",c(1:80)))

colors <- c("HLA_A_Count" = "#E41A1C", "HLA_B_Count" = "#377EB8", "HLA_C_Count" = "#4DAF4A")

p<-ggplot(df_long, aes(x = name, y = Count, fill = HLA_type)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Sample", y = "Count", fill = "HLA_type")+
  scale_fill_manual(values = colors)+
  # theme_minimal() +
  theme(
    # panel.border = element_rect(color = "black", size = 1.5),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, face = "bold",color = "black"),
    axis.title = element_text(size = 16, face = "bold",color = "black"),
    legend.title = element_text(size = 14, face = "bold",color = "black"),
    legend.text = element_text(size = 12, face = "bold",color = "black")
  )+scale_y_continuous(expand = c(0,0))+
  scale_y_continuous(limits = c(0, 6000), breaks = seq(0, 6000, by = 2000))
p
ggsave("./HLA_old/plot/nocm_3.pdf", plot = p, width = 18, height = 8,dpi=600)


df_long <- df_long %>%
  group_by(name) %>%
  mutate(total_count = sum(Count)) %>%
  ungroup()

# 计算每个HLA_type占比
df_long <- df_long %>%
  mutate(percentage = Count / total_count)

# 绘制堆叠图，Y轴用占比来表示
colors <- c("HLA_A_Count" = "#E41A1C", "HLA_B_Count" = "#377EB8", "HLA_C_Count" = "#4DAF4A")

p <- ggplot(df_long, aes(x = name, y = percentage, fill = HLA_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_blank(),  # 隐藏横坐标的标签
        axis.ticks.x = element_blank()) +  # 隐藏横坐标的刻度
  labs(x = "Sample", y = "Proportion", fill = "HLA_type") +
  scale_fill_manual(values = colors) +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, face = "bold", color = "black"),
    axis.title = element_text(size = 16, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, face = "bold", color = "black")
  ) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent)

p
ggsave("./HLA_old/plot/nocm_4.pdf", plot = p, width = 15, height = 8,dpi=600)

##统计
df_C_dominant <- df_long %>%
  group_by(name) %>%
  # 先过滤出每个name中三种HLA_type的Count
  dplyr::filter(HLA_type == "HLA_C_Count" & 
                  Count > max(Count[HLA_type == "HLA_A_Count"], Count[HLA_type == "HLA_B_Count"])) %>%
  ungroup()

# 计算符合条件的 name 数量
count_C_dominant <- nrow(df_C_dominant)#17个


df_A_dominant <- df_long %>%
  group_by(name) %>%
  # 先过滤出每个name中三种HLA_type的Count
  dplyr::filter(HLA_type == "HLA_A_Count" & 
                  Count > max(Count[HLA_type == "HLA_B_Count"], Count[HLA_type == "HLA_C_Count"])) %>%
  ungroup()

# 计算符合条件的 name 数量
count_C_dominant <- nrow(df_A_dominant)#4个

##


df_B_dominant <- df_long %>%
  group_by(name) %>%
  # 先过滤出每个name中三种HLA_type的Count
  dplyr::filter(HLA_type == "HLA_B_Count" & 
                  Count > max(Count[HLA_type == "HLA_A_Count"], Count[HLA_type == "HLA_C_Count"])) %>%
  ungroup()

# 计算符合条件的 name 数量
count_B_dominant <- nrow(df_B_dominant)#59个

##Redraw Figure 6I---------------------------------
peptide_freq<-read.csv("./HLA_old/nocm_peptide_freq.csv")

freq_table <- as.data.frame(table(peptide_freq$n))

colnames(freq_table) <- c("n", "Frequency")

freq_table$n <- as.numeric(as.character(freq_table$n))

p <- ggplot(freq_table, aes(x = n, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#3e7eb6", color = "black", width = 0.8) + # 调整宽度去掉间隔
  scale_x_reverse(breaks = seq(1, 80, by = 1), expand = c(0, 0)) + # 反转x轴并移除间隔
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1200), breaks = seq(0, 1200, by = 240)) + # 移除y轴的空白间隔
  theme_minimal() +
  labs(x = "Number of Individuals Presenting Peptides", y = "Peptides Num") +
  theme(
    axis.text.x = element_blank(), # 隐藏横坐标的标签
    axis.ticks.x = element_blank(), # 隐藏横坐标的刻度线
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 15, face = "bold", color = "black"),
    axis.title.y = element_text(size = 15, face = "bold", color = "black"),
    panel.grid.major = element_blank(), # 移除主网格线
    panel.grid.minor = element_blank(), # 移除次网格线
    panel.background = element_blank(),  # 移除背景
    axis.line = element_line(size = 0.8, color = "black") # 加粗坐标轴线
  )

# 移除 geom_text
# geom_text(aes(label = Frequency), vjust = -0.5, size = 3, color = "black", face = "bold", check_overlap = TRUE) 

p
ggsave("./HLA_old/plot/v2_nocm_pep_80.pdf", plot = p, width = 10, height = 6, dpi = 600)












