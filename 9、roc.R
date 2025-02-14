
####################################################数据处理####################################################################

##数据处理1：师姐的文件命名方式没有办法直接匹配AUC，所以得系统的给文件改个文件名：固定字符+模型+次数+AUC

setwd("F:/0A其他工作/3、刘睿娜师姐热项目——多界/7、多界预测2.0/feature_selection最终")

cl_list <- c(
             "archaea_fungi_bacteria_virus","archaea_fungi_virus","archaea_virus",
             "fungi_bacteria_virus","virus")

for(cl in cl_list){

setwd(paste0("F:/0A其他工作/3、刘睿娜师姐热项目——多界/7、多界预测2.0/feature_selection最终/",cl))
fileout <- list.files(pattern = "rocdata*")

for(filename00 in fileout){
  
  
  filename <- gsub("\\.txt$", "", filename00)
  # 使用 "_" 分割字符串
  split_parts <- strsplit(filename, "_")[[1]]
  
  # 提取分割后的部分
  prefix <- split_parts[1]
  value <- split_parts[2]
  part1 <- split_parts[3]
  part2 <- split_parts[4]
  
  # 按照指定顺序拼接字符串
  new_filename <- paste(prefix, part1, part2, value, sep = "_")
  new_filename <- paste0(new_filename, ".txt")
  
  old_path <- file.path(filename00)
  new_path <- file.path(new_filename)
  
  # 使用 file.rename 函数进行文件名更改
  file.rename(old_path, new_path)
}
}

##数据处理2：因为输出了所有的模型结果，但是画图我们就要前20个最优的，所以要写个代码把最优模型数据挑出来
library(stringr)
library(tools)
library(fs)
setwd("F:/0A其他工作/14、MTBI项目0812/3、多界预测模型/mTBI_code_Results_最终")
cl_list <- c("A","B","F","V",
             "AB","AF","AV","BV","FV","BF",
             "ABV","ABF","AFV","BFV",
             "ABFV")

for(cl in cl_list){

setwd(paste0("F:/0A其他工作/14、MTBI项目0812/3、多界预测模型/mTBI_code_Results_最终/",cl))
max_folder <- paste0("max/")     # 目标文件夹路径
newdir <- (paste0("F:/0A其他工作/14、MTBI项目0812/3、多界预测模型/mTBI_code_Results_最终/",cl,"/max"))

dir.create(newdir)

setwd(paste0("F:/0A其他工作/14、MTBI项目0812/3、多界预测模型/mTBI_code_Results_最终/",cl,"/roc"))
# 获取文件夹中所有符合模式的文件
files <- list.files(pattern = "roc*")

# 初始化一个空列表，用于存储每个模型的最高AUC值文件
max_files <- list()

# 遍历每个模型
for (i in 1:20) {
  # 获取当前模型的文件
  model_files <- grep(paste0("roc_", i, "_"), files, value = TRUE)
  
  # 初始化最高AUC值和对应的文件路径
  max_auc <- -Inf
  max_file <- NULL
  
  # 遍历该模型的所有文件
  for (file in model_files) {
    # 从文件名中提取AUC值,如果有报错，考虑是不是文件后缀名出错了，是txt还是csv
    auc <- as.numeric(str_match(file, "(\\d+\\.\\d+)\\.csv$")[, 2])
    
    # 检查是否为当前模型的最高AUC值
    if (auc > max_auc) {
      max_auc <- auc
      max_file <- file
    }
  }
  
  # 将最高AUC值文件添加到列表中
  max_files[[i]] <- max_file
  # 将最高AUC值文件拷贝到目标文件夹中
  file.copy(max_files[[i]], newdir)
}

}

#########################平均AUC#############################

setwd("F:/0A其他工作/14、MTBI项目0812/3、多界预测模型/mTBI_code_Results_最终")
cl_list <- c("A","B","F","V",
             "AB","AF","AV","BV","FV","BF",
             "ABV","ABF","AFV","BFV",
             "ABFV")
# 创建一个空的数据框用于存储结果
result_df <- data.frame(File = character(), AUC_Mean = numeric(), stringsAsFactors = FALSE)
# 循环处理每个文件
for (cl in cl_list) {
  # 读取文件
  data <- read.csv(file = paste0(cl,'/mTBI_Discovery_selfcv_repeat_max_',cl,'.csv',sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
  #sig <- read.csv(file = paste0(cl,'/CRC_Discovery_model_A3_sig_',cl,'.csv',sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
  # 计算 AUC 列均值
  auc_mean <- mean(data$AUC, na.rm = TRUE)
  # 提取文件名
  file_name <-cl
  #p_value <- sig[1,3]
  # 将文件名和 AUC 均值添加到结果数据框
  result_df <- rbind(result_df, data.frame(File = file_name, AUC_Mean = auc_mean,stringsAsFactors = FALSE))#,p=p_value
}
# 将结果写入新文件
write.csv(result_df, file = "mean_AUC_最终.csv", row.names = FALSE)



###################################  以上，数据处理完毕，开始画图  我们用的第二版代码 ##################################################




####ROC(第一版：使用原始数据大小绘图)

setwd("F:/0A其他工作/3、刘睿娜师姐项目/2、ROC曲线/feature_selection")

cl_list <- c("archaea","archaea_bacteria","archaea_bacteria_fungi","archaea_bacteria_virus","archaea_fungi",
             "archaea_fungi_bacteria_virus","archaea_fungi_virus","archaea_virus","bacteria","bacteria_virus",
             "fungi_bacteria","fungi_bacteria_virus","fungi_virus","virus")

for(cl in cl_list){
 
  fileout <- list.files(paste0(cl,"/"),pattern = "rocdata*")
  file1 <- data.frame()
  #file2 <- list()
  merged_df <- data.frame(matrix(ncol = 1, nrow = 55))#这个是用来调整行数的，nrow一般设置为所有文件的最大行数

  #mean <- data.frame(nrow(37))
  max <-  read.csv(paste0(cl,"/heatstress_Discovery_selfcv_repeat_metric_",cl,".csv"))#读取预测模型中，最佳预测结果文件

  
  for (m in 1:length(fileout)) {
  
  roc <- read.table(paste0(cl,"/",fileout[m]),sep = "\t", comment.char="",check.names=F,quot="",stringsAsFactors = F)#读取相应模型的roc表
  roc$group <- paste0(cl,"_",m)#组别
  roc$AUC <- max$AUC[m]#auc值

  file1 <- rbind(file1,roc)#纵向合并每一个roc表
  #file2[[m]] <- roc
  #mean <- merge
  mean <- roc
  if(nrow(mean)<55){#统一每一个表的长度，因为后面会删除空值，所以设大一点
    empty_rows <- 55 - nrow(mean)
    empty_data <- data.frame(matrix(ncol = ncol(mean), nrow = empty_rows))
    colnames(empty_data) <- colnames(mean)
    mean <- rbind(mean, empty_data)
  }
  merged_df <- cbind(merged_df,mean)
  }

  # 使用逻辑向量选择列
  sensitivities <- merged_df[, colnames(merged_df)%in%"sensitivities"]
  sensitivities$mean <- rowMeans(sensitivities)
  
  specificities <- merged_df[, colnames(merged_df)%in%"1-specificities"]
  specificities$mean <- rowMeans(specificities)
  
  AUC <- merged_df[, colnames(merged_df)%in%"AUC"]
  AUC$mean <- rowMeans(AUC)
  
  mean_pic <- cbind(sensitivities$mean,specificities$mean)
  mean_pic <- as.data.frame(mean_pic)
  mean_pic$group <- "Mean"
  mean_pic <- cbind(mean_pic,AUC$mean)
  
  colnames(mean_pic) <- c("sensitivities","specificities","group","AUC")
  
  data <- rbind(file1,mean_pic)
  
  data <- na.omit(data)
  
  write.csv(data,file = paste0("D:/project/SCZ/个性化图/ROC曲线/ROC_plot/roc_data_",cl,".csv"),quote = F,row.names = F)
  
  
  getPalette_family = colorRampPalette(brewer.pal(20, "Set3"))
  Family_col <- getPalette_family(20)#Family_col <- getPalette_family(10)
  Family_col <- c(Family_col,"red")
  
  roc_pic1 <- ggplot(data, aes(x = 1-specificities, y=sensitivities)) +
    geom_path(aes(color = group), linewidth=1.5, alpha = 0.7) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color="darkgrey", linetype=4)+
    scale_color_manual(values = Family_col)+
    #scale_color_npg(alpha=0.7)+
    annotate("text", x=0.3, y=0.3,label=paste0("Mean AUC = ",round(mean_pic$AUC[1],3)),size=5)+
    theme_bw() + 
    theme(element_text(color = "black",size=1),
          plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          legend.title=element_blank(), 
          legend.background = element_rect(fill=NULL, linewidth=0.5, 
                                           linetype="solid",colour ="black")
    )
  
  roc_pic1
  
  ggsave(path = paste0("D:/project/SCZ/个性化图/ROC曲线/ROC_plot"), filename =paste0("roc_",cl,".pdf"), roc_pic1, width = 120, height = 120, units="mm")
  ggsave(path = paste0("D:/project/SCZ/个性化图/ROC曲线/ROC_plot"), filename =paste0("roc_",cl,".png"), roc_pic1, width = 120, height = 120, units="mm")
  
}
# 


####ROC(第二版：取最小行数做标准，修剪每个文件大小使行数与最小行数相同，绘制曲线)


setwd("F:/0A其他工作/14、MTBI项目0812/3、多界预测模型/mTBI_code_Results_最终")
cl_list <- c("A","B","F","V",
             "AB","AF","AV","BV","FV","BF",
             "ABV","ABF","AFV","BFV",
             "ABFV")

for(i in 1:length(cl_list)){
  
  fileout <- list.files(paste0(cl_list[i],"/max/"),pattern = "roc*")#各个预测模型中的所有roc文件名
  file1 <- data.frame()
  
  min_lines <- Inf#设置最小行数
  
  max <-  read.csv(paste0(cl_list[i],"/mTBI_Discovery_selfcv_repeat_max_",cl_list[i],".csv"))#读取max文件，后面需要用到AUC值
  for (file_name in fileout) {#第一个内循环，先提取最小行数
    # 读取文件
    file_path <- file.path(paste0(cl_list[i],"/max/"), file_name)
    file_lines <- length(readLines(file_path))
    
    # 更新最小行数
    if (file_lines < min_lines) {
      min_lines <- file_lines
    }
  }
  min_lines=min_lines-1#最小行数减1是为了去掉列名行
  
  merged_df <- data.frame(matrix(ncol = 1, nrow = min_lines))#建立空的数据框
  
  for (m in 1:length(fileout)) {#第二个内循环，也是主体，用来生成绘制roc曲线的数据框
    
    roc <- read.csv(paste0(cl_list[i],"/max/",fileout[m]),comment.char="",check.names=F,quot="",stringsAsFactors = F)#读取相应模型的roc表
    roc$group <- paste0(cl_list[i],"_",m)#组别
    roc$AUC <- max$AUC[m]#auc值
    if(nrow(roc)>min_lines){#删除除首尾5行的任意几行，是行数满足最小行数
      rows_to_delete <- nrow(roc) - min_lines
      indices_to_delete <- sample(5:(nrow(roc) - 5),rows_to_delete)#生成随机行数
      roc <- roc[-indices_to_delete, ]
    }
    
    file1 <- rbind(file1,roc)#纵向合并每一个roc表
    
    merged_df <- cbind(merged_df,roc)#横向合并各个roc表，方面计算均值
  }
  #merged_df <- merged_df[,-1]
  #计算平均值
  # 使用逻辑向量选择列
  sensitivities <- merged_df[, colnames(merged_df)%in%"sensitivities"]
  sensitivities$mean <- rowMeans(sensitivities)
  
  specificities <- merged_df[, colnames(merged_df)%in%"specificities"]
  specificities$mean <- rowMeans(specificities)
  
  AUC <- merged_df[, colnames(merged_df)%in%"AUC"]
  AUC$mean <- rowMeans(AUC)
  
  mean_pic <- cbind(specificities$mean,sensitivities$mean)
  mean_pic <- as.data.frame(mean_pic)
  mean_pic$group <- "Mean"
  mean_pic <- cbind(mean_pic,AUC$mean)
  
  colnames(mean_pic) <- c("specificities","sensitivities","group","AUC")
  
  data <- rbind(file1,mean_pic)
  
  levels <- paste0(cl_list[i], "_", 1:20)#设置因子
  levels <- c(levels, "Mean")
  
  data$group <- factor(data$group, levels = levels, ordered = TRUE)
  
  
  write.csv(data,file = paste0("1、ROC_plot/roc_data_",cl_list[i],".csv"),quote = F,row.names = F)
  
  # getPalette_family = colorRampPalette(pal_nejm()(20))
  # # getPalette_family = colorRampPalette(brewer.pal(20, "Set3"))
  # Family_col <- getPalette_family(20)#Family_col <- getPalette_family(10)
  Family_col <- c(  "#6A3D9A", "#7570B3", "#BEBADA", "#BC80BD", "#DECBE4", "#D69A87", "#ED7953", "#FB9F3A",
                    "#F7D038", "#F0F921", "#D9EE42", "#ABED6E", "#6CE448", "#33A02C", "#A6CEE3", "#6DFDF2", 
                    "#8DD3C7", "#C2F7FF", "#FFFFCC", "#E5D8BD","red")
  
  data$legend <- "Mean"  #这一步的作用就是添加图例的时候可以将Mean作为单独的部分进行添加
  
  roc_pic1 <- ggplot(data, aes(x = 1-specificities, y=sensitivities)) +
    geom_path(aes(color = group,linetype=legend), linewidth=1.2, alpha = 0.7) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color="darkgrey", linetype=4)+
     scale_color_manual(values = Family_col)+
    # scale_color_futurama(a=0.7)+
    #scale_color_npg(alpha=0.7)+
    annotate("text", x=0.4, y=0.7,label=paste0("Mean AUC = ",round(mean_pic$AUC[1],3)),size=4)+
    guides(color=FALSE)+  #这一步是为了不展示每条roc曲线的图例
    guides(linetype = guide_legend(override.aes = list(linetype = 1, color = "red"))) + #这一步是为了设置Mean图例的颜色
    ggtitle(paste0("ROC of ",cl_list[i]))+
    theme_bw() + 
    theme(axis.text = element_text(color = "black"),
          legend.key.size = unit(0.4, "cm"),
          plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), 
          legend.position=c(.95, .25),
          legend.title=element_blank(), 
          legend.background = element_rect(fill=NULL, linewidth=0.2, 
                                           linetype="solid",colour ="black")
    )
  
  roc_pic1
  
  ggsave(path = paste0("1、ROC_plot"), filename =paste0("roc_",cl_list[i],".pdf"), roc_pic1, width = 110, height = 110, units="mm")
  ggsave(path = paste0("1、ROC_plot"), filename =paste0("roc_",cl_list[i],".png"), roc_pic1, width = 110, height = 110, units="mm")
  
}



####ROC(第三版：读取第二版的输出文件，将单一界和四界的ROC画在一张图里)
#1、首先读取单一界和四界的文件，只保留mean，然后合并

setwd("F:/0A其他工作/3、刘睿娜师姐热项目——多界/7、多界预测2.0/feature_selection最终/1、ROC_plot")

cl_list <- c("archaea","archaea_fungi_bacteria_virus","bacteria","virus")

pl_list <- c("A","ABFV","B","V")


result_df <- data.frame('1-specificities' = numeric(), sensitivities = numeric(), group=character(),AUC=numeric(),stringsAsFactors = FALSE)

for(i in 1:length(cl_list)){
  roc <- read.csv(paste0("roc_data_",cl_list[i],".csv"),sep = ",", comment.char="",check.names=F,quot="",stringsAsFactors = F)#读取相应模型的roc表
  roc <- subset(roc,group=="Mean")
  roc$group <- pl_list[i]
  result_df <- rbind(result_df,roc)
  
}

write.csv(result_df,file = "F:/0A其他工作/3、刘睿娜师姐热项目——多界/7、多界预测2.0/feature_selection最终/2、ROC_14/roc_1界_4界.csv",row.names = F)
#2、画图

result_df$group <- factor(result_df$group ,
                       levels = c("A","B","V","ABFV"),ordered = TRUE)


first_row_df <- result_df %>% 
  group_by(group) %>% 
  slice(1) 

legend_labels <- paste0(first_row_df$group, ": ", round(first_row_df$AUC,3))

getPalette_family = colorRampPalette(brewer.pal(5, "Set1"))
Family_col <- getPalette_family(5)#Family_col <- getPalette_family(10)

roc_pic2 <- ggplot(result_df, aes(x = `1-specificities`, y = sensitivities)) +
  geom_path(aes(color = group), linewidth = 1.2, alpha = 0.7) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "darkgrey", linetype = 4) +
  scale_color_manual(values = Family_col, labels = legend_labels) +
  #ggtitle(paste0("ROC of single and four kindoms")) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5),
    legend.justification = c(1, 0),
    legend.position = c(0.95, 0.25),
    legend.title = element_blank(),
    legend.background = element_rect(fill = NULL, linewidth = 0.2, linetype = "solid", color = "black")
  )
roc_pic2
ggsave(path = paste0("2、ROC_14"), filename =paste0("roc_1界_4界_notitle.pdf"), roc_pic2, width = 110, height = 110, units="mm")
ggsave(path = paste0("2、ROC_14"), filename =paste0("roc_1界_4界_notitle.png"), roc_pic2, width = 110, height = 110, units="mm")

