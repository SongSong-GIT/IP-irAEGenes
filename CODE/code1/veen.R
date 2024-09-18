
dir1 = "E:\\iraes-pathway工作总结\\test" 
file_list1 = list.files(path = dir1, pattern = "*..symbols.gmt",recursive = TRUE,full.names = TRUE) 
p <- NULL
venn_data <- list()
for(f in 1:length(file_list1)){
  knownGene <- clusterProfiler::read.gmt(file_list1[f])
  knownGene <- unique(knownGene[,2])
  venn_data[[f]] <- knownGene
  interestGene <- rownames(F_Rank)
  q <- length(intersect(interestGene, knownGene))-1
  m <- length(knownGene)
  n <- 25000-length(knownGene)
  k <- length(interestGene)
  p <- c(p, phyper(q, m, n, k, lower.tail = F))
}
res_phyper <- data.frame(set = c('c5.hpo', 'c6.all', 'c7.all'), PValue = p)
venn_data[[4]] <- rownames(F_Rank)
names(venn_data) <- c("Phenotype","Oncogenic","Immunologic","IMM-CAN-irAEs")


#venn
library(grid)
library(VennDiagram)
pdf(file = "E:\\iraes-pathway工作总结\\test/VennPlot.pdf")
venn.plot <- venn.diagram(venn_data,
                          filename = "E:\\iraes-pathway工作总结\\test/VennPlot.pdf",  
                          category =  c("Phenotype","Oncogenic","Immunologic","IMM-CAN-irAEs"),
                          fill = c('#F7F18D', '#C1DFF0', '#7BB6E8','#F5CD30'),
                          lty = "blank",
                          cex = 2,
                          cat.cex = 1.2,
                          cat.col = c('#F7F18D', '#C1DFF0', '#7BB6E8','#F5CD30'))
dev.off()
grid.draw(venn.plot)
