## Takes the 012 matracies with SNP info (first five columns) and retains only SNPs with 10 or more samples in at least two genotypes. 
library(data.table)
f.path = "/STORAGE/George/GTEx/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/split_by_chr/Matrix012/appended"
out.path = "/WORKSPACE/George/GTEx_genotype/Filtered"
f = list.files(path=f.path ,pattern=".txt$", full.names =T)

min_samples = 10

lapply(f, function(x) {df <- fread(x)
gr <- apply(df[,-c(1:5)],1, function(y) sum(table(y[y != -1]) > min_samples) >=2)
df <- df[gr,]
outfile = gsub(".+phg", "phg", x)
outfile = gsub("\\.txt", "_filter\\.txt", outfile)
outfile = paste(out.path, outfile, sep="/")
fwrite(df, outfile, row.names =F, sep ="\t")	
})
