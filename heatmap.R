library(ComplexHeatmap)
library(circlize)
stringAsFactors=FALSE
data<- read.csv("data/LX_dp4_dp6_ht.csv",header = TRUE,stringsAsFactors = F)
mtx = data[,-c(1,2,3)]
rownames(mtx) = data$peak
col_dp = c(rep('4',18),rep('6',36))
col_split_c = col_dp
col_gap_c = c(3)
row_split_c = data$LC
row_gap_c = c(3)

cola = columnAnnotation(
  dp = col_dp,
  col = list(
    DP = c('4'='#ff4e00','6'='#01cd74')
  ),
  show_annotation_name = T,
  na_col='white',
  simple_anno_size=unit(0.4, "cm")
  )
rowa = rowAnnotation(
  LC = data$LC,
  col = list(
    LC = c('L'='#00c4cc','X'='#6a3be4')
  ),
  show_annotation_name = T
)
col_fun=colorRamp2(c(1, 10),c("#fff200", "#003b64")) 

ht = Heatmap(mtx,
             col=col_fun,
             name = "-log10 p-value",
             na_col = '#eeeeee',
             show_row_dend = F,
             show_column_dend = F,
             cluster_columns = F,
             cluster_rows = F,
             show_row_names = T,
             show_column_names = T,
             column_split = col_split_c,
             column_gap=unit(col_gap_c, "mm"),
             row_split = row_split_c,
             row_gap=unit(row_gap_c, "mm"),
             row_names_gp = gpar(fontsize = 10),
             bottom_annotation=cola,
             left_annotation = rowa,
             row_names_side = 'left', 
             column_names_gp = gpar(fontsize = 10, rotation=90),
             rect_gp = gpar(col = "white", lwd = 2)
)
ht


data = read.csv("simulated_data/positive/noise_1_10.csv",header = F,stringsAsFactors = F)

data = data[,-1]

# col_fun=colorRamp2(c(0, 1),c("#fff200", "#003b64")) 
sum(data>0.1)
mat = as.matrix(data[1:100,])
mat = as.matrix(data)
ht = Heatmap(mat,
             # col=col_fun,
             name = "test",
             na_col = '#eeeeee',
             show_row_dend = F,
             show_column_dend = F,
             cluster_columns = F,
             cluster_rows = F,
             show_row_names = F,
             show_column_names = F
)
ht

