setwd("/home/thwang12/ME-440/taeyoung") # modify to your project directory ("Module_6")

# Read GTF files
gtf <- read.table("./Genome/gencode.vM30.annotation.gtf", header=F, sep="\t", stringsAsFactors = F, colClasses = c("character", "NULL", "character", "numeric", "numeric", "NULL", "character", "NULL", "character"))
colnames(gtf) <- c("chrom", "type", "start", "end", "strand", "annot")

# Hold only gene features
table(gtf$type)
gtf <- subset(gtf, type=="gene")
gtf$type <- NULL

# Parse "annot" column
gtf.annot <- strsplit(gtf$annot, split=";")
select.annot <- function(x) {
  res <- sapply(gtf.annot, function(t) {y<-grep(x, t, value=T); ifelse(is.null(y), NA, y)})
  res <- gsub(x, "", res)
  res <- gsub(" ", "", res)
  return(res)
}
gtf$gene_id <- select.annot("gene_id")
gtf$gene_type <- select.annot("gene_type")
gtf$gene_name <- select.annot("gene_name")
gtf$mgi_id <- select.annot("mgi_id")
gtf$annot <- NULL

# Read sampleList (a list of sample names)
# First, copy the "/data/lgoff2/ME-440/taeyoung/sampleList" text file to your project directory.
sampleList <- scan("./sampleList", what="character", sep="\n")

# Read featureCounts output
fCounts <- data.frame()
for (s in sampleList) {
  s.file <- paste0("./3.Count/", s, ".fCounts")
  temp <- read.table(s.file, header=T, sep="\t", stringsAsFactors = F)
  if (length(fCounts)==0) {
    fCounts <- temp[, -ncol(temp)]
  }
  fCounts[,s] <- temp[, ncol(temp)]
}
colnames(fCounts)[1:6] <- c("gene_id", "chrome", "start", "end", "strand", "length")
all(fCounts$gene_id == gtf$gene_id)

# Save
save(gtf, fCounts, file="./Results/results_count.rda")
