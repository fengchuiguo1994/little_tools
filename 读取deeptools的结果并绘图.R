gf = gzfile("patski.regions.gz", 'rt')
data = read.table(gf, header = F, sep = '\t',skip=1)
mydatflt = data[,c(7:ncol(data))]
quantile(unlist(mydatflt), probs = c(0.1,0.9,0.95,0.98,0.999), na.rm = FALSE)