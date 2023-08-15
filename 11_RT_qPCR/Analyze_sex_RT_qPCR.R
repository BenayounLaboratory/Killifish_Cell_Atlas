setwd('/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/CODE/11_RT_qPCR')
options(stringsAsFactors = F)

library(readxl)
library(beeswarm)
library(outliers)

# 2023-04-25
# Analyze RT-qPCR

######################################################################
# read excel sheet
my.data         <- read_xlsx('~/Downloads/2023-06_22_Master_qPCR_Dynamic_Threshold.xlsx', sheet = 1)

# subset for stats
my.data.liver   <- my.data[my.data$Tissue == "Liver",]
my.data.muscle  <- my.data[my.data$Tissue == "Muscle",]
my.data.spleen  <- my.data[my.data$Tissue == "Spleen",]


######################################################################
# ncFem1

# one male ncFem spleen might be outlier (> 4 times the mean malee value), test with Grubbs test
grubbs.test(na.omit(my.data.spleen$ncRNA_exp[my.data.spleen$Sex == "Male"]), type = 10)
# Grubbs test for one outlier
# data:  na.omit(my.data.spleen$ncRNA_exp[my.data.spleen$Sex == "Male"])
# G = 2.431954, U = 0.034384, p-value = 5.148e-05
# alternative hypothesis: highest value 4.64322268890957 is an outlier
##### ==> observed high male data point in spleen is an outlier

# Do stats
nc.liver.pval  <- wilcox.test(na.omit(my.data.liver$ncRNA_exp[my.data.liver$Sex == "Female"]),
                             na.omit(my.data.liver$ncRNA_exp[my.data.liver$Sex == "Male"]))
nc.muscle.pval <- wilcox.test(na.omit(my.data.muscle$ncRNA_exp[my.data.muscle$Sex == "Female"]),
                             na.omit(my.data.muscle$ncRNA_exp[my.data.muscle$Sex == "Male"]))
nc.spleen.pval <- wilcox.test(na.omit(my.data.spleen$ncRNA_exp[my.data.spleen$Sex == "Female"]),
                             na.omit(my.data.spleen$ncRNA_exp[my.data.spleen$Sex == "Male"])[-1]) # removes the 4.64 outlier

# make data.frame for plotting
my.data.ncFem        <- my.data[,c("Tissue","Sex","ncRNA_exp")]
my.data.ncFem        <- data.frame(my.data.ncFem[!is.na(my.data.ncFem$ncRNA_exp),]) # remove NAs
my.nc.outlier        <- which(my.data.ncFem$ncRNA_exp == max(na.omit(my.data.spleen$ncRNA_exp[my.data.spleen$Sex == "Male"])))
my.data.ncFem.no_out <- my.data.ncFem[-my.nc.outlier,]
  
pdf(paste0(Sys.Date(),"_ncFem1_tissue_RT_qPCR_outlier_removed.#pdf"), height = 4, width = 5)
boxplot( ncRNA_exp ~ Sex + Tissue   , data = my.data.ncFem.no_out, outline = F, las = 2,
         col = c("deeppink", "deepskyblue"), xlab = '', main = "ncFem1", ylim = c(0,3),
         ylab = "Relative ncFem1 expression (A.U.)")
beeswarm(ncRNA_exp ~ Sex + Tissue   , data = my.data.ncFem.no_out, add = T, pch = 16)
text(1.5,2.8,signif(nc.liver.pval$p.value ,2), cex = 0.8)
text(3.5,2.8,signif(nc.muscle.pval$p.value,2), cex = 0.8)
text(5.5,2.8,signif(nc.spleen.pval$p.value,2), cex = 0.8)
dev.off()

######################################################################
# LINE-R2

# Do stats
line.liver.pval  <- wilcox.test(na.omit(my.data.liver$LINE_exp[my.data.liver$Sex == "Female"]),
                              na.omit(my.data.liver$LINE_exp[my.data.liver$Sex == "Male"]))
line.muscle.pval <- wilcox.test(na.omit(my.data.muscle$LINE_exp[my.data.muscle$Sex == "Female"]),
                              na.omit(my.data.muscle$LINE_exp[my.data.muscle$Sex == "Male"]))
line.spleen.pval <- wilcox.test(na.omit(my.data.spleen$LINE_exp[my.data.spleen$Sex == "Female"]),
                              na.omit(my.data.spleen$LINE_exp[my.data.spleen$Sex == "Male"])) 

# plotting
pdf(paste0(Sys.Date(),"_LINER2_tissue_RT_qPCR.#pdf"), height = 4, width = 5)
boxplot( LINE_exp ~ Sex + Tissue   , data = my.data, outline = F, las = 2,
         col = c("deeppink", "deepskyblue"), xlab = '', main = "LINE-R2", ylim = c(0,3),
         ylab = "Relative LINE-R2 expression (A.U.)")
beeswarm(LINE_exp ~ Sex + Tissue   , data = my.data, add = T, pch = 16)
text(1.5,2.8,signif(line.liver.pval$p.value ,2), cex = 0.8)
text(3.5,2.8,signif(line.muscle.pval$p.value,2), cex = 0.8)
text(5.5,2.8,signif(line.spleen.pval$p.value,2), cex = 0.8)
dev.off()


######################################################################

######################################################################
# Hpx

# Do stats
hpx.liver.pval  <- wilcox.test(na.omit(my.data.liver$Hpx_exp[my.data.liver$Sex == "Female"]),
                                na.omit(my.data.liver$Hpx_exp[my.data.liver$Sex == "Male"]))
hpx.muscle.pval <- wilcox.test(na.omit(my.data.muscle$Hpx_exp[my.data.muscle$Sex == "Female"]), 
                                na.omit(my.data.muscle$Hpx_exp[my.data.muscle$Sex == "Male"]))  
hpx.spleen.pval <- wilcox.test(na.omit(my.data.spleen$Hpx_exp[my.data.spleen$Sex == "Female"]),
                                na.omit(my.data.spleen$Hpx_exp[my.data.spleen$Sex == "Male"])) 

# Plotting
pdf(paste0(Sys.Date(),"_Hpx_tissue_RT_qPCR.#pdf"), height = 4, width = 5)
boxplot( Hpx_exp ~ Sex + Tissue   , data = my.data, outline = F, las = 2, log = "y",
         col = c("deeppink", "deepskyblue"), xlab = '', main = "Hpx", ylim = c(0.1,5000),
         ylab = "Relative Hpx expression (A.U.)")
beeswarm(Hpx_exp ~ Sex + Tissue   , data = my.data, add = T, pch = 16)
text(1.5,4000,signif(hpx.liver.pval$p.value ,2), cex = 0.8)
text(3.5,4000,signif(hpx.muscle.pval$p.value,2), cex = 0.8)
text(5.5,4000,signif(hpx.spleen.pval$p.value,2), cex = 0.8)
dev.off()
######################################################################

######################################################################
#zp3
# Do stats
zp3.liver.pval  <- wilcox.test(na.omit(my.data.liver$zp3_exp[my.data.liver$Sex == "Female"]),
                               na.omit(my.data.liver$zp3_exp[my.data.liver$Sex == "Male"]))
zp3.muscle.pval <- wilcox.test(na.omit(my.data.muscle$zp3_exp[my.data.muscle$Sex == "Female"]), 
                               na.omit(my.data.muscle$zp3_exp[my.data.muscle$Sex == "Male"]))  
zp3.spleen.pval <- wilcox.test(na.omit(my.data.spleen$zp3_exp[my.data.spleen$Sex == "Female"]),
                               na.omit(my.data.spleen$zp3_exp[my.data.spleen$Sex == "Male"])) 

# Plotting
pdf(paste0(Sys.Date(),"_zp3_tissue_RT_qPCR.#pdf"), height = 4, width = 5)
boxplot( zp3_exp ~ Sex + Tissue   , data = my.data, outline = F, las = 2, log = "y",
         col = c("deeppink", "deepskyblue"), xlab = '', main = "zp3", ylim = c(0.001,100),
         ylab = "Relative zp3 expression (A.U.)")
beeswarm(zp3_exp ~ Sex + Tissue   , data = my.data, add = T, pch = 16)
text(1.5,100,signif(zp3.liver.pval$p.value ,2), cex = 0.8)
text(3.5,100,signif(zp3.muscle.pval$p.value,2), cex = 0.8)
text(5.5,100,signif(zp3.spleen.pval$p.value,2), cex = 0.8)
dev.off()

##################################
#apoa1a
# Do stats
apoa1a.liver.pval  <- wilcox.test(na.omit(my.data.liver$apoa1a_exp[my.data.liver$Sex == "Female"]),
                               na.omit(my.data.liver$apoa1a_exp[my.data.liver$Sex == "Male"]))
apoa1a.muscle.pval <- wilcox.test(na.omit(my.data.muscle$apoa1a_exp[my.data.muscle$Sex == "Female"]), 
                               na.omit(my.data.muscle$apoa1a_exp[my.data.muscle$Sex == "Male"]))  
apoa1a.spleen.pval <- wilcox.test(na.omit(my.data.spleen$apoa1a_exp[my.data.spleen$Sex == "Female"]),
                               na.omit(my.data.spleen$apoa1a_exp[my.data.spleen$Sex == "Male"])) 

# Plotting
pdf(paste0(Sys.Date(),"_apoa1a_tissue_RT_qPCR.#pdf"), height = 4, width = 5)
boxplot( apoa1a_exp ~ Sex + Tissue   , data = my.data, outline = F, las = 2, log = "y",
         col = c("deeppink", "deepskyblue"), xlab = '', main = "apoa1a", ylim = c(0.001,4000),
         ylab = "Relative apoa1a expression (A.U.)")
beeswarm(apoa1a_exp ~ Sex + Tissue   , data = my.data, add = T, pch = 16)
text(1.5,4000,signif(apoa1a.liver.pval$p.value ,2), cex = 0.8)
text(3.5,4000,signif(apoa1a.muscle.pval$p.value,2), cex = 0.8)
text(5.5,4000,signif(apoa1a.spleen.pval$p.value,2), cex = 0.8)
dev.off()





