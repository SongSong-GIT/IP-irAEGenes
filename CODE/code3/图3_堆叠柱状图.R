setwd("E:\\iraes-pathway工作总结\\irAEs6\\result2\\result_finallll2_rdata")
###############################
  # 初始化变量
  cancer <- c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "LIHC",
              "LUAD", "LUSC", "PAAD", "PRAD", "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC")
  
  name <- unique(statistic10$pathway)
  all_da <- data.frame(pathway = character(), num = numeric(), cancer = character())
  # 设置包含 RData 文件的目录
  dir <- "E:\\iraes-pathway工作总结\\irAEs6\\result2\\result_finallll2_rdata/"
  file_list <- list.files(path = dir, pattern = "*.RData", recursive = TRUE, full.names = TRUE)
  # 处理每个文件
  for (f in 1:length(file_list)) {
    load(file_list[f])
    # 为当前癌症类型初始化数据框
    data <- data.frame(pathway = name, num = 0, cancer = cancer[f])
    # 计算每个通路的出现次数
    for (i in name) {
      data$num[data$pathway == i] <- sum(statistic10$pathway == i)
    }
    data <- data[order(data$num, decreasing = TRUE), ]
    data <- head(data, 3)
    all_da <- rbind(all_da, data)
  }
  write.csv(all_da,"图3_堆叠柱状图.csv",row.names = F)
####################  画图    ######################################
  
  library(ggplot2)  # 导入ggplot2包，用于绘图
  library(dplyr)    # 导入dplyr包，用于数据处理
  library(readxl)   # 导入readxl包，用于读取Excel文件
  library(RColorBrewer) # 导入RColorBrewer包，用于设置颜色板
  
  # 设置主题为黑白主题
  theme_set(theme_bw())
  
  # 设置工作目录
  setwd("E:\\iraes-pathway工作总结\\irAEs\\图3")  
  
  # 读入Excel文件数据，并进行数据清洗和准备
  data <- all_da
  
  data %>%
    # 名称列转换成因子变量并按照指定的水平顺序排列，便于绘图
    mutate(Sample = factor(cancer, levels=c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "LIHC",
                                            "LUAD", "LUSC", "PAAD", "PRAD", "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC")),
           # 等外显子数列转换成因子变量并按照指定的水平顺序排列，便于绘图
          ) -> data01
  
  # 开始绘图，使用ggplot()函数设置映射
  ggplot(data01, aes(x = Sample, y = num, fill = pathway)) +
    
    # 绘制柱形图，使用geom_bar()函数，统计标识符为"identity"
    # 参数position="stack"表示将不同Exon_num对应的数量堆叠在一起，参数width=0.5表示柱子的宽度
    geom_bar(stat = "identity", position = "fill", width = 0.5) +
    
    # 设置颜色板，使用scale_fill_brewer()函数
    scale_fill_manual(values = c("Antigen_Processing_and_Presentation"="#5B9BD5","Antimicrobials"="#E64B35",
                                 "BCRSignalingPathway"="#9A60B4","Chemokine_Receptors"="#91CC75",
                                 "Chemokines"="#73C0DE","Cytokine_Receptors"="#FC8452",
                                 "NaturalKiller_Cell_Cytotoxicity"="#7ECBD1","TGFb_Family_Member_Receptor"="#D05E15",
                                 "Interleukins_Receptor"="#D88778","Cytokines"="#EA7CCC","TCRsignalingPathway"="#0E6F6F",
                                 "Interferons"="#FAC858","TNF_Family_Members"="#239873","TNF_Family_Members_Receptors"="#4E87B5",
                                 "TGFb_Family_Member"="#7C2245","Interleukins"="#3BA272"))+
    
    # 将x轴和y轴翻转，便于展示,将它注释掉，将垂直展示
    #coord_flip() +
    
    # 添加数值标签，使用geom_text()函数和position_stack()函数设置标签的位置和在图形中的表现形
    
    # 设置坐标轴标签
    labs(x = "cancer", y = "pathway", title = "") + 
    
    # 设置图表主题，包括去除网格线和边框，设置坐标轴线的颜色
    theme(panel.grid.major = element_blank(),  # 去除主网格线
          panel.grid.minor = element_blank(),  # 去除次网格线
          panel.border = element_blank(),      # 将图表的外边框（也称为panel）设置为空白
          axis.line = element_line(colour = "black", linewidth = 0.5)) # 坐标轴线条的颜色为黑色，粗细为0.5  
  

















