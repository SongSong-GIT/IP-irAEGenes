
hh=import("E:\\iraes-pathway工作总结\\irAEs\\图3\\kegg.result.csv")

h=hh
setwd("E:\\iraes-pathway工作总结\\irAEs\\图3")
h$order=factor(rev(as.integer(rownames(h))),labels = rev(h$Description))
ggplot(h,aes(y=order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "red",high ="blue" )+#颜色自己可以换
  labs(title = "KEGG UP Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
ggsave("./KEGG.bar.pdf")

df=hh[,c("ID","Description","p.adjust","Count","GeneRatio")]
df=df[order(df$Count,decreasing = T),]
df$Description=str_to_sentence(df$Description)
df$Description=factor(df$Description,df$Description)

## 气泡图
ggplot(df,aes(x=Count, y=Description)) + 
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+  #点的大小根据Count数变化，颜色根据P值变化
  scale_colour_gradient(low="blue", high="red")+  #设置图例
  labs(
    color=expression(-log[10](p.adjust)),
    size="Num of Genes",
    x="Num of Genes")+
  theme_bw()+  #设置背景
  theme(
    axis.text.y = element_text(size = rel(1)),
    axis.title.x = element_text(size=rel(1.0)),
    axis.title.y = element_blank())
ggsave("./KEGG.point.pdf")
