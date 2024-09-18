install.packages("fmsb")

max_min <- data.frame(
  BLCA= c(1, 0), BRCA = c(1, 0), CESC = c(1, 0), CHOL = c(1, 0), COAD = c(1, 0), ESCA = c(1, 0),GBM =  c(1, 0), HNSC =  c(1, 0), KIRC = c(1, 0),LIHC = c(1, 0),
  LUAD= c(1, 0) ,LUSC= c(1, 0), PAAD= c(1, 0),PRAD= c(1, 0), READ= c(1, 0) ,SARC= c(1, 0) ,SKCM= c(1, 0) ,STAD= c(1, 0) ,THCA= c(1, 0), UCEC= c(1, 0)
  )
rownames(max_min) <- c("Max", "Min")

# 合并数据
df <- rbind(max_min, data)
df
library(fmsb)
#student1_data <- df[c("Max", "Min", "Student.1"), ]
radarchart(
  df, axistype = 1,
  # Customize the polygon
  pcol = "#4985B3", pfcol = scales::alpha("#4985B3", 0.5), plwd = 2, plty = 1,
  # Customize the grid
  cglcol = "grey", cglty = 1, cglwd = 1,
  # Customize the axis
  axislabcol = "grey", 
  # Variable labels
  vlcex = 0.7, vlabels = colnames(df),
  caxislabels = c(0, 0.01, 0.05, 0.07, 1))


# 创建一个包含20组数据名称的向量
group_names <- c('BLCA', 'BRCA', 'CESC' ,'CHOL' ,'COAD'  ,'ESCA'  ,'GBM' ,'HNSC' , 'KIRC' ,'LIHC' , 'LUAD' ,'LUSC', 'PAAD','PRAD', 'READ' ,'SARC' ,'SKCM' ,'STAD' ,'THCA', 'UCEC')




