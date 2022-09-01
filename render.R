files <- list.files(".", recursive = TRUE)
files <- files[grepl("^(modules|prereqs).*\\.Rmd$", files, ignore.case = TRUE)]
lapply(files, rmarkdown::render)
