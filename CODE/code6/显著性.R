top10gene=import("E:\\iraes-pathway工作总结\\图\\生存\\图\\top10gene.csv")

dir = "E:\\iraes-pathway工作总结\\图\\生存\\生存结果"
file_list = list.files(path = dir, pattern = "*",recursive = TRUE,full.names = TRUE)  
result_all <- data.frame()
for (f in 1:20) {
  exp_data <- read.csv(file_list[f], row.names = 1)
  
  if (f == 1) {
    file_path <- file_list[f]
    file_name <- gsub(".*/([^/]+)_.*", "\\1", file_path)
    file_name <- gsub("_.*", "\\", file_name)
    
    k <- data.frame(gene=rownames(exp_data),
                    p=exp_data$PVal,
                    type=file_name) 

    result_all <- rbind(result_all, k)
    
  } else {
    file_path <- file_list[f]
    file_name <- gsub(".*/([^/]+)_.*", "\\1", file_path)
    file_name <- gsub("_.*", "\\", file_name)
    k <- data.frame(gene=rownames(exp_data),
                    p=exp_data$PVal,
                    type=file_name) 
    
    result_all <- rbind(result_all, k)
  }
}
length(unique(result_all$gene))

data=merge(result_all,top10gene,by="gene")
type_counts <- table(data$type)
# 找出出现次数大于4的值
selected_values <- names(type_counts[type_counts > 4])
# 选出在 type 列中值为 selected_values 的行
selected_rows <- data[data$type %in% selected_values, ]


library(ggplot2)
library(reshape2)

selected_rows$value <- ifelse(selected_rows$p < 0.01, "***", ifelse(selected_rows$p < 0.05, "*", ""))

ggplot(data = selected_rows, aes(type,gene ,fill = p))+
  # geom_rect(color = "white")+
  geom_point(aes(size = count), shape = 21)+
  scale_fill_gradientn(colors = c("white","#d94c4c"),limits=c(0,0.05)) +  
  scale_size(range = c(10,12),
             breaks = c(5,6),
             labels = c(5,6),
             limits = c(5,6))+
  theme_minimal()+ 
  #geom_text(aes(label=value))+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()

write.csv(selected_rows,"显著性图.csv",row.names = F,quote = F)