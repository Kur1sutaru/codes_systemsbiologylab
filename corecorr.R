setwd("~/Daniel Rna seq data/06-07 correlation plots")

## To read the data
corrall <- read.csv("~/Daniel Rna seq data/06-07 correlation plots/corrall.csv")
## To visualize data
library("ggpubr")
ggscatter(corrall, x = "genesigall", y = "datarsrd", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Genes signature", ylab = "RSRD")

## 
## Shapiro-Wilk normality test
shapiro.test(corrall$datarsrd)
shapiro.test(corrall$genesigall)
shapiro.test(corrall$genesighigh)
shapiro.test(corrall$genesiglow)
shapiro.test(corrall$hrd)
shapiro.test(corrall$trex1)
shapiro.test(corrall$rnaseh2a)


## Graphical normality test
ggqqplot(corrall$rnaseh2a, ylab = "RNASEH2A score")

## CORRELATION MATRIX - ALL DATA
## Pearson correlation
cor(corrall$genesighigh,corrall$trex1, method = "pearson")
cor.test(corrall$genesighigh,corrall$trex1, method = "pearson")

## Spearman correlation
cor(dados$PROD,dados$Trat, method = "spearman")
cor.test(dados$PROD,dados$Trat, method = "spearman")

cor(dados,dados, method = "spearman")
cor.test(dados,dados, method = "spearman")

# The spearman correlation makes rank of the lowset to biggest values
# Only make sense to use for rank values

# To estimate the  #p-value
#  cor.mtest <- function(mat, ...) {
#    mat <- as.matrix(mat)
#    n <- ncol(mat)
#    p.mat<- matrix(NA, n, n)
#    diag(p.mat) <- 0
#    for (i in 1:(n - 1)) {
#      for (j in (i + 1):n) {
#        tmp <- cor.test(mat[, i], mat[, j], ...)
#        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
#      }
#    }
#    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
#    p.mat
#  }
# p.mat
#p.mat <- cor.mtest(test)
# dim(p.mat)


# Pearson - applied to parametric data, and Spearman to non parametric data.



## To visualize the data
ggscatter(dados1, x = "Trat", y = "PROD", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tratamento", ylab = "PROD")


###############################################################################
# Other correlation analysis                                                  #
###############################################################################
#Correlation plots - tutorial retrieved from 
#https://r-coder.com/correlation-plot-r/#:~:text=Correlation%20plots%2C%20also%20known%20as,the%20correlation%20between%20continuous%20variables.
setwd("~/Daniel Rna seq data/06-07 correlation plots")

library(psych)

# Sample data
# Numerical variables
data <- corrall
# Factor variable (groups)
groups <- data[colnames(data)] 


# Plot correlation matrix
pairs(data)


# Equivalent but using the plot function
plot(data)


# Function to add histograms
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  his <- hist(x, plot = FALSE)
  breaks <- his$breaks
  nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = rgb(0, 1, 1, alpha = 0.5), ...)
  # lines(density(x), col = 2, lwd = 2) # Uncomment to add density lines
}

# Creating the scatter plot matrix
pairs(data,
      upper.panel = NULL,         # Disabling the upper panel
      diag.panel = panel.hist)    # Adding the histograms


# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

# Plotting the correlation matrix
pairs(data,
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth) # Smoothed regression lines


# The corPlot function creates a graph of a correlation matrix, 
# coloring the regions by the level of correlation.

corPlot(data, cex = 1.2)

# Recall to type ?corPlot for additional arguments and details.


# Correlogram with corrgram and corrplot packages

####  corrgram function ######
### On the one hand, the corrgram package calculates the correlation of the data
# and draws correlograms. The function of the same name allows customization via
# panel functions. As an example, you can create a correlogram in R where the 
# upper panel shows pie charts and the lower panel shows shaded boxes with the 
# following code:

install.packages("corrgram")
library(corrgram)

corrgram(data,
         order = TRUE,              # If TRUE, PCA-based re-ordering
         upper.panel = panel.pie,   # Panel function above diagonal
         lower.panel = panel.shade, # Panel function below diagonal
         text.panel = panel.txt,    # Panel function of the diagonal
         main = "Correlogram")      # Main title

# There are several panel functions that you can use.
# Using the apropos function you can list all of them:

apropos("panel.")


install.packages("corrplot")
library(corrplot)

corrplot(cor(data),        # Correlation matrix
         method = "shade", # Correlation plot method
         type = "full",    # Correlation plot style (also "upper" and "lower")
         diag = TRUE,      # If TRUE (default), adds the diagonal
         tl.col = "black", # Labels color
         bg = "white",     # Background color
         title = "",       # Main title
         col = NULL)       # Color palette

### To trade the color palette use the colorRampPalette function

# This function also allows clustering the data. The clustering methods 
# according to the documentation are: "original" (default order), "AOE" 
# (angular order of eigenvectors), "FPC" (first principal component order), 
# "hclust" (hierarchical clustering order) and "alphabet" (alphabetical order).
# If you chose hierarchical clustering you can select between the following
# methods: "ward", "ward.D", "ward.D2", "single", "complete", "average", 
# "mcquitty", "median" and "centroid". In this case, you can also create 
# clustering with rectangles. An example is shown in the following block of code:

corrplot(cor(data),
         method = "circle",       
         order = "hclust",         # Ordering method of the matrix
         hclust.method = "ward.D", # If order = "hclust", is the cluster method to be used
         addrect = 2,              # If order = "hclust", number of cluster rectangles
         rect.col = 3,             # Color of the rectangles
         rect.lwd = 3)             # Line width of the rectangles





