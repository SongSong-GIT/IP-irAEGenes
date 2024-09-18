# 读取第一组CSV数据
data1 <- read.csv("E:/iraes-pathway工作总结/test/luad.csv")

# 读取第二组CSV数据
data2 <- read.csv("E:/iraes-pathway工作总结/test/lusc.csv")

# 提取两组数据中的基因集合
luad <- as.character(data1$gene)
lusc <- as.character(data2$gene)

# 计算交集数量
intersection <- length(intersect(luad, lusc))

# 计算并集数量
union <- length(union(luad, lusc))

# 计算Jaccard相似度
jaccard_similarity <- intersection / union

# 打印结果
print(jaccard_similarity)



###################################################################################################


# 创建一个列表来存储20组子数据
#data_list <- list(esetData1[[1]], esetData1[[2]], ..., esetData1[[20]])

# 提取每组子数据的行名
row_names_list <- lapply(esetData1, function(data) row.names(data))

# 初始化一个空的相似性矩阵
similarity_matrix <- matrix(0, nrow = 20, ncol = 20)


# 计算Jaccard相似性并填充矩阵
for (i in 1:20) {
  for (j in 1:20) {
    if (i != j) {
      intersect_count <- length(intersect(row_names_list[[i]], row_names_list[[j]]))
      union_count <- length(union(row_names_list[[i]], row_names_list[[j]]))
      jaccard_similarity <- intersect_count / union_count
      similarity_matrix[i, j] <- jaccard_similarity
    }
  }
}

# 将结果保存为数据框
similarity_df <- as.data.frame(similarity_matrix)

# 设置行和列名
colnames(cross_table)
rownames(similarity_df) <-colnames(cross_table)
colnames(similarity_df) <- colnames(cross_table)
##绘图##
p4 <- pheatmap(similarity_df,  
               color = colorRampPalette(c('#E3E3E3','#98F5FF' ,'#FF6347'))(400), 
               border_color = "NA",  
               scale = "none", 
               cluster_rows = FALSE, 
               cluster_cols = FALSE, 
               legend = TRUE, 
               #legend_breaks = c(-1, 0, 1), 
               #legend_labels = c("low","","heigh"), 
               show_rownames = TRUE, 
               show_colnames = TRUE, 
               fontsize = 8,
               display_numbers = TRUE,  #是否显示每个色块对应的数值(经归一化后的数值)
               number_format = "%.2f",#数值格式，%.2f表示保留小数点后两位，
               #%.1e表示使用科学计数法并保留小数点后一位
               number_color = "black",  #设置数值颜色
               fontsize_number = 10    #设置数值的字体大小
)
pdf("热图4.pdf",width = 8,height = 8)
p4
dev.off()



















































