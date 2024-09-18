
setwd("E:\\iraes-pathway工作总结\\irAEs\\图1")
dir = "E:/iraes-pathway工作总结/irAEs6/result2/result_finallll_rdata/" 
file_list = list.files(path =dir, pattern = "*.RData",recursive = TRUE,full.names = TRUE) 
cancer=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KIRC","LIHC","LUAD","LUSC","PAAD","PRAD","READ",
         "SARC","SKCM","STAD","THCA","UCEC")
statistic <- data.frame(cancer = character(), PathwayAppear = character(), GeneNum = numeric())
for (i in 1:length(file_list)) {
  #i=2
load(file_list[i])
tablePathway <- table(TCGA_symbolCodCor_filter[,4])
countPathway <- data.frame(Symbol = names(tablePathway),
                           PathwayNum = as.vector(tablePathway))   
countPathway <- countPathway[order(countPathway$PathwayNum,decreasing = T),]
statistic1 <- data.frame(cancer = rep(cancer[i],length(countPathway$Symbol)),
                         PathwayAppear = countPathway$Symbol,
                         GeneNum = as.vector(countPathway$PathwayNum))
statistic <- rbind(statistic, statistic1)
}
write.csv(statistic, file = "E:\\iraes-pathway工作总结\\irAEs\\图1/statisticALL.csv", quote = F)

dot_data <- read.csv("./statisticALL.csv")
summary(dot_data$GeneNum)
(p1<-ggplot(dot_data, aes(x=cancer, y=PathwayAppear)) +
    geom_point(aes(size = GeneNum, color = GeneNum), shape = 16) +
    scale_color_gradient(low = '#7FFFD4',high = '#00BFFF', name = 'Number of irAEs-related Gene')+
    scale_y_discrete(limits = sortss[order(as.numeric(sortss[,2]), decreasing = T),1],
                     breaks = sortss[order(as.numeric(sortss[,2]), decreasing = T),1],
                     labels = sortss[order(as.numeric(sortss[,2]), decreasing = T),1])+
    scale_size(range = c(1,8),  breaks = c(10,20,40,60,80), name = NULL)+
    theme_test() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = 'bold', color = 'black'),
          axis.text = element_text(face = 'bold', color = 'black'),
          axis.title.y = element_text(face = 'bold', color = 'black', size = 10),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(face = 'bold', color = 'black',size = 10),
          legend.title = element_text(face = 'bold', color = 'black',size = 10, hjust = 0.5),
          legend.position = 'top',
          panel.border = element_blank()))

result <- statistic %>%
  group_by(PathwayAppear) %>%
  summarize(
    cancer = paste(unique(cancer), collapse = ", "),  # 合并 cancer 列的值
    GeneNum = sum(GeneNum)  # 计算 GeneNum 的总和
  )
result <- result[order(result$GeneNum),]

p4 <- ggplot(result) +
  geom_col(aes(x = seq(0.1, 2, length.out = nrow(result)), y = GeneNum, fill = GeneNum),
           color = 'black', 
           width = 0.095) + 
  scale_fill_gradientn(colors = c("#7FFFD4", "#00BFFF")) +  # 颜色渐变从蓝色到红色
  labs(x = NULL, y = 'Gene count') +
  scale_y_continuous(breaks = seq(0, max(result$GeneNum), length.out = 6)) +
  scale_x_continuous(limits = c(0.05, 2.05), breaks = seq(0.1, 2, length.out = nrow(result))) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(face = 'bold', color = 'black', angle = 45, size = 8, vjust = 0.5),
        axis.title.x = element_text(face = 'bold', color = 'black', size = 10),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, 'cm')) +
  coord_flip()

# 显示图形
print(p4)


























