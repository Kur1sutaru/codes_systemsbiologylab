setwd("~/DRUG REPOSITIONING MMRD")
# import package
library("ggplot2")



# plot: dot plot
ggplot(data = mismatch.repair.inducer.drugs, aes(x = Cell.Line, y = Perturbagen, 
                        color = q.value..BH., size = Dose)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8)) +
  ggtitle("Mismatch repair inducer drug enrichment analysis")