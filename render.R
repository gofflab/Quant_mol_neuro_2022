files <- list.files(".", pattern = "*.Rmd$", recursive = TRUE)
for (file in files) {
    rmarkdown::render(file)
}
