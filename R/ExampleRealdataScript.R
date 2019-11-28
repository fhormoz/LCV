# LCV example script written by Katie Siewert
if(!require(optparse)){
    install.packages("optparse")
    library(optparse)
}
option_list = list(
  make_option("--sum_stat_name1", type="character", default="s1.txt",
              help="First GWAS summary statistics file name", metavar="character"),
  make_option("--sum_stat_name2", type="character", default="s2.txt",
              help="Second GWAS summary statistics file name", metavar="character"),
  make_option("--ld_score_file_name", type="character", default="ldscorefile",
              help="the file of ldscore", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Start with data munged using the ldsc package
trait1File = opt$sum_stat_name1;
trait2File = opt$sum_stat_name2;
ldscoresFile = opt$ld_score_file_name

print(trait1File)
print(trait2File)

#Load trait 1 data and calculate Zs
d1 = na.omit(read.table(gzfile(trait1File),header=TRUE,sep="\t",stringsAsFactors = FALSE))

#Load trait 2 data and calculate Zs
d2 = na.omit(read.table(gzfile(trait2File),header=TRUE,sep="\t",stringsAsFactors = FALSE))

#Load LD scores
d3=read.table(gzfile(ldscoresFile), header=TRUE, sep='\t', stringsAsFactors=FALSE)

#Merge
m = merge(d3,d1,by="SNP")
data = merge(m,d2,by="SNP")

#Sort by position
data = data[order(data[,"CHR"],data[,"BP"]),]

#Flip sign of one z-score if opposite alleles-shouldn't occur with UKB data
#If not using munged data, will have to check that alleles match-not just whether they're opposite A1/A2
mismatch = which(data$A1.x!=data$A1.y,arr.ind=TRUE)
data[mismatch,]$Z.y = data[mismatch,]$Z.y*-1
data[mismatch,]$A1.y = data[mismatch,]$A1.x
data[mismatch,]$A2.y = data[mismatch,]$A2.x


#Run LCV-need to setwd to directory containing LCV package
source("RunLCV.R")
print(sum(data$Z.x))
print(sum(data$Z.y))
LCV = RunLCV(data$L2,data$Z.x,data$Z.y)
sprintf("Estimated posterior gcp=%.2f(%.2f), log10(p)=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, log(LCV$pval.gcpzero.2tailed)/log(10), LCV$rho.est, LCV$rho.err)

