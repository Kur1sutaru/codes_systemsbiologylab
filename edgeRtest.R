library("edgeR")


setwd("D:/aline/GSE4290")


files <- c("epilepGSM97800.CEL",
           "epilepGSM97803.CEL",
           "epilepGSM97804.CEL",
           "epilepGSM97805.CEL",
           "epilepGSM97807.CEL",
           "epilepGSM97809.CEL",
           "epilepGSM97811.CEL",
           "epilepGSM97812.CEL",
           "epilepGSM97816.CEL",
           "epilepGSM97817.CEL",
           "epilepGSM97820.CEL",
           "epilepGSM97825.CEL",
           "epilepGSM97827.CEL",
           "epilepGSM97828.CEL",
           "epilepGSM97833.CEL",
           "epilepGSM97834.CEL",
           "epilepGSM97840.CEL",
           "epilepGSM97846.CEL",
           "epilepGSM97848.CEL",
           "epilepGSM97849.CEL",
           "epilepGSM97850.CEL",
           "epilepGSM97853.CEL",
           "epilepGSM97855.CEL",
           "GSM300168.CEL",
           "GSM300169.CEL",
           "GSM300170.CEL",
           "GSM300171.CEL",
           "GSM300172.CEL",
           "GSM300185.CEL",
           "GSM300187.CEL",
           "GSM300190.CEL",
           "GSM300193.CEL",
           "GSM300197.CEL",
           "GSM300205.CEL",
           "GSM300209.CEL",
           "GSM300215.CEL",
           "GSM300219.CEL",
           "GSM300223.CEL",
           "GSM300227.CEL",
           "GSM300231.CEL",
           "GSM300235.CEL",
           "GSM300239.CEL",
           "GSM300243.CEL",
           "GSM300262.CEL",
           "GSM300268.CEL",
           "GSM300272.CEL",
           "GSM300280.CEL")

DG <- readDGE(files,header=TRUE)

#To create a data frame now

#To load each sample individually

epilep97800 = read.table('epilepGSM97800.CEL', sep = '\t',
                        header = T)
epilepGSM97803 = read.table('epilepGSM97803.CEL', sep = '\t',
                           header = F)
epilepGSM97804 = read.table('epilepGSM97804.CEL', sep = '\t',
                            header = F)
epilepGSM97805 = read.table('epilepGSM97805.CEL', sep = '\t',
                            header = F)
epilepGSM97807 = read.table('epilepGSM97807.CEL', sep = '\t',
                            header = F)
epilepGSM97809 = read.table('epilepGSM97809.CEL', sep = '\t',
                            header = F)
epilepGSM97811 = read.table('epilepGSM97811.CEL', sep = '\t',
                         header = F)
epilepGSM97812 = read.table('epilepGSM97812.CEL', sep = '\t',
                           header = F)
epilepGSM97816 = read.table('epilepGSM97816.CEL', sep = '\t',
                            header = F)
epilepGSM97817 = read.table('epilepGSM97817.CEL', sep = '\t',
                            header = F)
epilepGSM97820 = read.table('epilepGSM97820.CEL', sep = '\t',
                            header = F)
epilepGSM97825 = read.table('epilepGSM97825.CEL', sep = '\t',
                           header = F)
epilepGSM97827 = read.table('epilepGSM97827.CEL', sep = '\t',
                            header = F)
epilepGSM97828 = read.table('epilepGSM97828.CEL', sep = '\t',
                            header = F)
epilepGSM97833 = read.table('epilepGSM97833.CEL', sep = '\t',
                            header = F)
epilepGSM97834 = read.table('epilepGSM97834.CEL', sep = '\t',
                            header = F)
epilepGSM97840 = read.table('epilepGSM97840.CEL', sep = '\t',
                            header = F)
epilepGSM97846 = read.table('epilepGSM97846.CEL', sep = '\t',
                           header = F)
epilepGSM97848 = read.table('epilepGSM97848.CEL', sep = '\t',
                            header = F)
epilepGSM97849 = read.table('epilepGSM97849.CEL', sep = '\t',
                            header = F)
epilepGSM97850 = read.table('epilepGSM97850.CEL', sep = '\t',
                           header = F)
epilepGSM97853 = read.table('epilepGSM97853.CEL', sep = '\t',
                           header = F)
epilepGSM97855 = read.table('epilepGSM97855.CEL', sep = '\t',
                            header = F)
GSM300168 = read.table('GSM300168.CEL', sep = '\t',
                      header = F)
GSM300169 = read.table('GSM300169.CEL', sep = '\t',
                      header = F)
GSM300170 = read.table('GSM300170.CEL', sep = '\t',
                       header = F)
GSM300171 = read.table('GSM300171.CEL', sep = '\t',
                       header = F)
GSM300172 = read.table('GSM300172.CEL', sep = '\t',
                       header = F)
GSM300185 = read.table('GSM300185.CEL', sep = '\t',
                      header = F)
GSM300187 = read.table('GSM300187.CEL', sep = '\t',
                       header = F)
GSM300190 = read.table('GSM300190.CEL', sep = '\t',
                       header = F)
GSM300193 = read.table('GSM30093.CEL', sep = '\t',
                       header = F)
GSM300197 = read.table('GSM300197.CEL', sep = '\t',
                       header = F)
GSM300205 = read.table('GSM300205.CEL', sep = '\t',
                       header = F)
GSM300209 = read.table('GSM300209.CEL', sep = '\t',
                       header = F)
GSM300215 = read.table('GSM300215.CEL', sep = '\t',
                       header = F)
GSM300219 = read.table('GSM300219.CEL', sep = '\t',
                       header = F)
GSM300223 = read.table('GSM300223.CEL', sep = '\t',
                       header = F)
GSM300227 = read.table('GSM300227.CEL', sep = '\t',
                       header = F)
GSM300231 = read.table('GSM300231.CEL', sep = '\t',
                       header = F)
GSM300235 = read.table('GSM300235.CEL', sep = '\t',
                       header = F)
GSM300239 = read.table('GSM300239.CEL', sep = '\t',
                       header = F)
GSM300243 = read.table('GSM300243.CEL', sep = '\t',
                       header = F)
GSM300262 = read.table('GSM300262.CEL', sep = '\t',
                       header = F)
GSM300268 = read.table('GSM300268.CEL', sep = '\t',
                       header = F)
GSM300272 = read.table('GSM300272.CEL', sep = '\t',
                       header = F)
GSM300280 = read.table('GSM300280.CEL', sep = '\t',
                       header = F)

# To combine data sets into a matrix
geneCounts = data.frame(epilep97800 [,1], epilepGSM97803 [,2], epilepGSM97804 [,2], epilepGSM97805 [,2], epilepGSM97807 [,2], epilepGSM97809 [,2], epilepGSM97811 [,2], epilepGSM97812 [,2], epilepGSM97816 [,2], epilepGSM97817 [,2], epilepGSM97820 [,2], epilepGSM97825 [,2], epilepGSM97827 [,2], epilepGSM97828 [,2], epilepGSM97833 [,2], epilepGSM97834 [,2], epilepGSM97840 [,2], epilepGSM97846 [,2], epilepGSM97848 [,2], epilepGSM97849 [,2], epilepGSM97850 [,2], epilepGSM97853 [,2], epilepGSM97855 [,2], GSM300168 [,2], GSM300169 [,2], GSM300170 [,2], GSM300171 [,2], GSM300172 [,2], GSM300185 [,2], GSM300187 [,2], GSM300190 [,2], GSM300193 [,2], GSM300197 [,2],  GSM300205 [,2], GSM300209 [,2], GSM300215 [,2], GSM300219 [,2], GSM300223 [,2], GSM300227 [,2], GSM300231 [,2], GSM300235 [,2], GSM300239 [,2], GSM300243 [,2], GSM300262 [,2], GSM300268 [,2], GSM300272 [,2], GSM300280 [,2])
row.names(geneCounts) = MOCK1[,1]
sizeGeneCounts = dim(geneCounts)
geneCounts = geneCounts[1:(sizeGeneCounts[1]-5),]
condition = c(rep('epilep97800',2), rep('GSM300272', 2))
sampleNames = c('epilep97800', 'epilepGSM97803', 'epilepGSM97804', 'epilepGSM97805', 'epilepGSM97807', 'epilepGSM97809', 'epilepGSM97811', 'epilepGSM97812', 'epilepGSM97816', 'epilepGSM97817', 'epilepGSM97820', 'epilepGSM97825', 'epilepGSM97827', 'epilepGSM97828', 'epilepGSM97833', 'epilepGSM97834', 'epilepGSM97840', 'epilepGSM97846', 'epilepGSM97848', 'epilepGSM97849', 'epilepGSM97850', 'epilepGSM97853', 'epilepGSM97855', 'GSM300168', 'GSM300169', 'GSM300170', 'GSM300171', 'GSM300172', 'GSM300185', 'GSM300187', 'GSM300190', 'GSM300193', 'GSM300197',  'GSM300205', 'GSM300209', 'GSM300215', 'GSM300219', 'GSM300223', 'GSM300227', 'GSM300231', 'GSM300235', 'GSM300239', 'GSM300243', 'GSM300262', 'GSM300268', 'GSM300272', 'GSM300280')
colnames(geneCounts) = sampleNames
View(geneCounts)

#To build the generalized linear model that will be used for differential expression 
dge <- DGEList(counts=geneCounts, group=condition)
design <- model.matrix(~condition+0, data=dge$samples)
colnames(design) = gsub("condition","",colnames(design))

# To perform the normalization 
dge <- calcNormFactors(dge)
dge$samples
plotMDS(dge)

#To write a plot (PCA)
jpeg(file="epilepsy_vs_normal_hippocampus.jpeg", width=5000, height=5000, units="px", res=300)
plotMDS(dge)
dev.off()

#To estimate the dispersion
disp <- estimateGLMCommonDisp(dge, design)
disp <- estimateGLMTrendedDisp(disp, design)
disp <- estimateGLMTagwiseDisp(disp, design)
plotBCV(disp)

#To analyze the differential expression and likelihood ratio
fit <- glmQLFit(disp,design)
qlf <- glmQLFTest(fit,coef=2)

fit <- glmFit(disp, design)
lrt <- glmLRT(fit,coef=2)
top<-topTags(lrt)
is.data.frame(top)
#To confirm if conditions are correct
colnames(design)

#To perform the exact test
dgeTest <- exactTest(disp)
dgeTest

#To see the histogram of P-Values
hist(dgeTest$table[,"PValue"], breaks=50)


#To tell edgeR what comparison you want to perform
epixnorm = makeContrasts(EPI-NORM, levels=design)


##### EPILEPSY x NORMAL ####


#To perform the differential expression testing for that comparison
lrt.epixnorm = glmLRT(fit, contrast=epixnorm)
res_epixnorm<-topTags(lrt.epixnorm, n=500000, sort.by = "p.value")
write.csv(res_epixnorm, file="epixnorm.csv")

#To write a table with the significative p-values
sig_epixnorm<-topTags(lrt.epixnorm, n=500000, p.value=.05, sort.by = "p.value")
write.csv(sig_epixnorm, file="epixnorm_sig.csv")
