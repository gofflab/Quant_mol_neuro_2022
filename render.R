files <- list.files(".", recursive = TRUE)
files <- files[grepl("^(modules|prereqs).*\\.Rmd$", files)]
lapply(files, rmarkdown::render)
