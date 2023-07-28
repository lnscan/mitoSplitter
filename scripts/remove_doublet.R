## Distinguishing between singlets and doubles according to Favg Gaussian fitted distribution

library(data.table)
library(mixtools) ## R >= 4.0.0
library(stringr)

args <- commandArgs(T)
merfile <- args[1]
outdir <- args[2]
removebcf <- args[3]

merdf0 <- fread(merfile, header = T, sep = "\t")
merdf <- as.data.frame(merdf0)
merdf <- merdf[,-1]
merdf[is.na(merdf)] <- 0

if (file.size(removebcf) == 0){ 
    removebcs <- c() 
}else{
    removebcs00 <- read.table(removebcf, header = F, stringsAsFactors = F)
    removebcs0 <- removebcs00$V1
    removebcs <- str_replace(removebcs0, "-", ".")
}

## cal Favg 
subdf <- merdf[, !(colnames(merdf) %in% removebcs)]    
sums <- colSums(subdf)
counts <- apply(subdf != 0, 2, sum)
favg <- sums / counts
        
## fit Gaussian distribution
out<-normalmixEM(favg, k=2, fast = TRUE)
pdf(paste0(outdir, "/Favg_Gaussian_fitted_distribution.pdf"), width = 5, height = 5)
print(plot(out, 2)) 
dev.off()
        
## get uniroot
d1 <- function(x) dnorm(x, out$mu[1], out$sigma[1]) * out$lambda[1]
d2 <- function(x) dnorm(x, out$mu[2], out$sigma[2]) * out$lambda[2]

intersection <- uniroot(function(x) dnorm(x, out$mu[1], out$sigma[1]) * out$lambda[1] - dnorm(x, out$mu[2], out$sigma[2]) * out$lambda[2], interval = c(0, 1)) 
thre <- intersection$root

singletlist0 <- names(favg[favg>thre])
singletlist <- str_replace(singletlist0, "\\.", "-")
write.table(singletlist, paste0(outdir, "/Favg_Gaussian_fitted_singlet.list"), quote = F, row.names = F, col.names = F)

doubletlist0 <- names(favg[favg<=thre])
doubletlist <- str_replace(doubletlist0, "\\.", "-")
write.table(doubletlist, paste0(outdir, "/Favg_Gaussian_fitted_doublet.list"), quote = F, row.names = F, col.names = F)
