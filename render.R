files <- list.files(".", recursive = TRUE)
files <- files[grepl("^(modules|prereqs).*\\.Rmd$", files, ignore.case = TRUE)]
print(files)
lapply(files, rmarkdown::render)
