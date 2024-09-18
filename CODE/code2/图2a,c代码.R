types = read.table("E:/iraes-pathway工作总结/test/type.txt")
colnames(types) <- types[1, ]
types = colnames(types)
count1 <- 0
count2 <- 0
count3 <- 0
count <-data.frame(type=character(),target=numeric(),value=numeric(),gene=character())
load(file = "E:/iraes-pathway工作总结/irAEs4//UCEC_esetData1.RData")
for(i in 1:length(esetData1)){
  #i=1
  # name1 = paste0("E:/iraes-pathway工作总结/test/group/", types[i], "_group")
  violinData_g123 = esetData1[[i]]
  violinData_g123$gene = rownames(violinData_g123)
  violinData_g123 = data.frame(type = rep(types[i], dim(violinData_g123)[1]),
                               violinData_g123)  

count = rbind(count,violinData_g123)
   
  
}
count1 <- count[count$target==1,]
count2 <- count[count$target==2,]
count3 <- count[count$target==3,]

l1=length(count1$gene)
l2=length( count2$gene)
l3=length(count3$gene)

countall=l1+l2+l3

percent1 <- (l1 / countall) * 100
percent2 <- (l2 / countall) * 100
percent3 <- (l3 / countall) * 100

# 创建一个数据框
data <- data.frame(Category = c("group1", "group2", "group3"),
                   Percentage = c(percent1, percent2, percent3))
library(ggplot2)

# 创建饼图
pdf(file = "E:/iraes-pathway工作总结/test/group123g.pdf", width = 6, height = 4)
ggplot(data, aes(x = "", y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values=c("#91CBCF","#EBD7A5","#F8CCCD"))+
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(Percentage), "%")), position = position_stack(vjust = 0.5)) +
  labs(fill = "Immune regulons") +
  theme_void()
dev.off()
setwd("E:\\iraes-pathway工作总结\\irAEs\\图2")
write.table(data,"图2c.txt",row.names = F,quote = F,sep = "\t")

########################### 1    #######################################
library(ggplot2)
library(cowplot)

barData_G1=data.frame(Types=names(table(count1$type)),Count=as.numeric(table(count1$type)))
barData_G2=data.frame(Types=names(table(count2$type)),Count=as.numeric(table(count2$type)))
barData_G3=data.frame(Types=names(table(count3$type)),Count=as.numeric(table(count3$type)))

#color=c("#4BB4BB", "#DDAC2A", "#FF6E6E")
color=c("#91CBCF","#EBD7A5","#F8CCCD")
pdf(file = "E:\\iraes-pathway工作总结\\irAEs\\图2/braplot_G1.pdf", width = 6, height = 4)
ggplot(barData_G1, aes(Types, weight = Count, fill = color[1])) +
  geom_bar(fill = color[1], color = "black", width = .7, position = 'stack', linewidth = .1) +
  labs( y = 'Number') +
  scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black", vjust = 0.5),
        axis.title.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black", linewidth = .1),
        axis.line = element_line(color = "black", linewidth = .1),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        plot.margin = unit(c(1,1,1,1),'cm'),
        legend.position = 'none')
dev.off()
pdf(file = "E:\\iraes-pathway工作总结\\irAEs\\图2/braplot_G2.pdf", width = 6, height = 4)
ggplot(barData_G2, aes(Types, weight = Count, fill = color[2])) +
  geom_bar(fill = color[2], color = "black", width = .7, position = 'stack', linewidth = .1) +
  labs( y = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = .1),
        axis.line = element_line(color = "black", linewidth = .1),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        plot.margin = unit(c(1,1,1,1),'cm'),
        legend.position = 'none')
dev.off()
pdf(file = "E:\\iraes-pathway工作总结\\irAEs\\图2/braplot_G3.pdf", width = 6, height = 4)
ggplot(barData_G3, aes(Types, weight = Count, fill = color[3])) +
  geom_bar(fill = color[3], color = "black", width = .7, position = 'stack', linewidth = .1) +
  labs( y = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = .1),
        axis.line = element_line(color = "black", linewidth = .1),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        plot.margin = unit(c(1,1,1,1),'cm'),
        legend.position = 'none')
dev.off()

write.table(barData_G1,"图2a.txt",row.names = F,quote = F,sep = "\t")
write.table(barData_G2,"图2a1.txt",row.names = F,quote = F,sep = "\t")
write.table(barData_G3,"图2a2.txt",row.names = F,quote = F,sep = "\t")

write.table(count1,"Mild.txt",quote = F,row.names = F,sep = "\t")
write.table(count2,"Active.txt",quote = F,row.names = F,sep = "\t")
write.table(count3,"Silent.txt",quote = F,row.names = F,sep = "\t")

# 加载ggplot2包

