setwd('/Users/berenice/Dropbox/Manuscripts_and_Publications/2023/2023_Killifish_single_cell_Atlas_Bryan/CODE/9_GRZ_Microscopy/')
options(stringsAsFactors = F)

library(readxl)
library(beeswarm)
library(outliers)


######################################################################
# read excel sheet
my.data         <- read_xlsx('2023-04-27_Killi_oil_red_area_2cohorts.xlsx', sheet = 1)

# subset for stats
my.data.GRZ   <- my.data[my.data$Strain == "GRZ",]
my.data.ZMZ   <- my.data[my.data$Strain == "ZMZ",]


######################################################################
# Do stats
oro.GRZ.pval  <- wilcox.test(my.data.GRZ$ORO_norm_to_F[my.data.GRZ$Sex   == "F"],
                              my.data.GRZ$ORO_norm_to_F[my.data.GRZ$Sex == "M"])
oro.ZMZ.pval <- wilcox.test(my.data.ZMZ$ORO_norm_to_F[my.data.ZMZ$Sex    == "F"],
                              my.data.ZMZ$ORO_norm_to_F[my.data.ZMZ$Sex == "M"])

my.data.plotting <- list("GRZ_F" = my.data.GRZ$ORO_norm_to_F[my.data.GRZ$Sex   == "F"],
                         "GRZ_M" = my.data.GRZ$ORO_norm_to_F[my.data.GRZ$Sex   == "M"],
                         "ZMZ_F" = my.data.ZMZ$ORO_norm_to_F[my.data.ZMZ$Sex   == "F"],
                         "ZMZ_M" = my.data.ZMZ$ORO_norm_to_F[my.data.ZMZ$Sex   == "M"])

# make data.frame for plotting
pdf(paste0(Sys.Date(),"_Liver_ORO_normaalized_to_FCohortMedian.pdf"), height = 5, width = 6)
boxplot( rev(my.data.plotting), outline = F, las = 1,
         col = c("deepskyblue","deeppink"), 
         ylab = '', main = "ORO liver signal", ylim = c(0,10),
         xlab = "Relative ORO Signal (A.U.)", horizontal=TRUE)
beeswarm(rev(my.data.plotting), add = T, pch = 16, cex = 1.2, horizontal=TRUE)
text(9, 3.5,signif(oro.GRZ.pval$p.value,2))
text(9, 1.5,signif(oro.ZMZ.pval$p.value,2))
abline(v = 1, col = "red", lty = "dashed")
dev.off()

######################################################################
#######################
sink(file = paste(Sys.Date(),"_RT_qPCR_session_Info.txt", sep =""))
sessionInfo()
sink()
