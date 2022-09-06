files <- list.files(".", recursive = TRUE)
files <- files[grepl("^(modules|prereqs).*\\.Rmd$", files, ignore.case = TRUE)]
for (file in files) {
    rmarkdown::render(file)
}
