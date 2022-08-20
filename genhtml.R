files <- list.files(".", pattern = "*.html$", recursive = TRUE)
links <- paste(
    lapply(files, function(x) paste('<a href="/Quant_mol_neuro_2022/', x, '">', x, "</a>", sep = "")),
    collapse = "\n<br>\n"
)
html <- paste(
    "<html><head><title>Quant_mol_neuro_2022</title></head><body><h1>Quantitative Molecular Neuroscience 2022</h1><h2>Course Materials></h2>",
    links,
    "</body></html>",
    sep = ""
)
writeLines(html, "index.html")
