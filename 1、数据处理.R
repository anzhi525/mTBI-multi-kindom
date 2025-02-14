##多界分析前数据处理

#####第一步，合并菌表------

setwd("F:/0A其他工作/14、MTBI项目0812/3、多界预测模型")


cl_list <- c("archaea","bacteria","fungi","virus")

# 创建一个空的数据框 hebing 用于存储合并后的数据
hebing <- NULL

# 循环遍历 cl_list 中的文件
for (file_name in cl_list) {

  data <- read.delim(paste0("out_kraken2_",file_name,"_species.txt"),  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
  data$kindom <- file_name
  
  # 如果 hebing 数据框为空，直接将当前文件数据添加到 hebing 中
  if (is.null(hebing)) {
    hebing <- data
  } else {
    # 检查当前文件数据的列名与 hebing 中的列名是否相同
    if (all(colnames(data) %in% colnames(hebing))) {
      # 如果列名相同，则按行合并数据
      hebing <- rbind(hebing, data)
    } else {
      
    }

    
    write.csv(hebing,file = "out_kraken2_four_kindom.csv",quote = F,row.names = F)
  }
}

#####第二步，提取差异菌的菌表---------------

diff_zong <- read_excel("差异菌表-zong.xlsx")
diff_zong <- diff_zong[,1:2]
colnames(diff_zong)[2] <- "clade_name"


out_species <- read.delim("F:/0A其他工作/14、MTBI项目0812/3、多界预测模型/out_species_x100_0.001.txt")


otutab <- out_species[out_species$clade_name%in%diff_zong$clade_name,]

otutab00 <- inner_join(diff_zong,otutab,by="clade_name")

####第三步，输出各个预测模型对应菌表---------------

  pl_list <- c("archaea","bacteria","fungi","virus",
               "archaea  bacteria","archaea  fungi","archaea  virus","bacteria  virus","fungi  virus","bacteria  fungi",
               "archaea  bacteria  virus","archaea  bacteria  fungi","archaea  fungi  virus","bacteria  fungi  virus",
               "archaea  bacteria  fungi  virus")  
  
  pl_name <- c("A","B","F","V",
               "AB","AF","AV","BV","FV","BF",
               "ABV","ABF","AFV","BFV",
               "ABFV")  
               
  for(i in 1:length(pl_list)){
    pl2 <- as.character(unlist(strsplit(pl_list[i], split = "  ")))
    otutab01 <- subset(otutab00, kindom %in% pl2)
    otutab01 <- otutab01[, !grepl("kindom", names(otutab01))]
    #otutab01$name <- make.names(otutab01$name)
    write.table(otutab01,file =paste0("./otutab/otu_",pl_name[i],".txt"),quote = F,row.names = F,col.names = T,sep = "\t")
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
  write.csv(result_df, file = "mean_AUC.csv", row.names = FALSE)
  
  
  
#############################将合适模型的结果结果提取到CRC_code_Results数据集中----------------------
#1、首先拷贝所有的feature_com.csv,都在CRC_code_Results2021文件夹下

  setwd("F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results")
  cl_list <- c("A","B","F","V",
               "AB","AF","AV","BV","FV","BF",
               "ABV","ABF","AFV","BFV",
               "ABFV")
  for(cl in cl_list){
    dir.create(cl,recursive = TRUE)
    file.copy(paste0( "F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results2021/",cl,"/feature_com_",cl,".csv"), cl,overwrite=T)
    file.copy(paste0( "F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results2021/",cl,"/CRC_Discovery_feature_selection_importance_com_",cl,".pdf"), cl,overwrite=T)
    
  }

  
 # 2、拷贝预测模型结果
  
  setwd("F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results")
  cl_list2021 <- c("ABFV","V","ABV","AFV","ABF")
  for(cl in cl_list2021){
    dir.create(cl,recursive = TRUE)
    file.copy(paste0( "F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results2021/",cl,"/CRC_Discovery_selfcv_repeat_max_",cl,".csv"), cl,overwrite=T)
    file.copy(paste0( "F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results2021/",cl,"/CRC_Discovery_selfcv_repeat_max_para_",cl,".csv"), cl,overwrite=T)
    file.copy(paste0( "F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results2021/",cl,"/CRC_Discovery_selfcv_repeat_metric_",cl,".csv"), cl,overwrite=T)
    source_folder <- paste0("F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results2021/",cl,"/roc")
    destination_folder <- paste0("F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results/",cl)
    file.copy(from = source_folder, to = destination_folder, recursive = TRUE)
    
    
  }
  
  

  setwd("F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results")
  cl_list666 <- c("A","B","F",
                  "AB","AF","AV","BV","FV","BF",
                  "BFV")
  for(cl in cl_list666){
    dir.create(cl,recursive = TRUE)
    file.copy(paste0( "F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results666/",cl,"/CRC_Discovery_selfcv_repeat_max_",cl,".csv"), cl,overwrite=T)
    file.copy(paste0( "F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results666/",cl,"/CRC_Discovery_selfcv_repeat_max_para_",cl,".csv"), cl,overwrite=T)
    file.copy(paste0( "F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results666/",cl,"/CRC_Discovery_selfcv_repeat_metric_",cl,".csv"), cl,overwrite=T)
    source_folder <- paste0("F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results666/",cl,"/roc")
    destination_folder <- paste0("F:/0A其他工作/11、MTBI项目/4、差异菌的四界分析/CRC_code_Results/",cl)
    file.copy(from = source_folder, to = destination_folder, recursive = TRUE)
    
    
  }
  
  
  
  
  #####################配色方案####################
  2
  "#440154", "#482878", "#3E4989", "#31688E", "#26828E", "#1F9E89", "#35B779", "#6DCD59",
  "#B4DD2B", "#FDE725", "#F9F871", "#F3E881", "#E6C88A", "#D69A87", "#CB6B60", "#D53E4F",
  "#EA0F6B", "#F5317F", "#FF6E54", "#FFC594"
  
  3
  "#6A3D9A", "#7570B3", "#BEBADA", "#BC80BD", "#DECBE4", "#D69A87", "#ED7953", "#FB9F3A",
  "#F7D038", "#F0F921", "#D9EE42", "#ABED6E", "#6CE448", "#33A02C", "#A6CEE3", "#6DFDF2", 
  "#8DD3C7", "#C2F7FF", "#FFFFCC", "#E5D8BD"
  
  1
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", 
  "#6A3D9A", "#B15928", "#A6CEE3", "#B2DF8A"
  
  
  4
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
  "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"
  
  5
  "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
  "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666"