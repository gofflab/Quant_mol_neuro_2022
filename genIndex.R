files <- list.files(".", recursive = TRUE)
files <- files[grepl("^(Module|prereqs).*\\.html", files)]

links <- paste(
    lapply(files, function(x) paste("- [", x, "](/Quant_mol_neuro_2022/", x, ")", sep = "")),
    collapse = "\n"
)

write(paste(
    readChar("README.md", file.info("README.md")$size),
    "\n### Course Materials",
    links,
    sep = "\n"
), file = "README.Rmd")

rmarkdown::render("README.Rmd", output_file = "index.html")
unlink("README.Rmd")
