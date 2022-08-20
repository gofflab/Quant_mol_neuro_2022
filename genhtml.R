files <- list.files(".", pattern = "*.html$", recursive = TRUE)
links <- paste(
    lapply(files, function(x) paste('<a href="/Quant_mol_neuro_2022/', x, '">', x, "</a>", sep = "")),
    collapse = "\n<br>\n"
)
html <- paste(
    "<html><head><title>Quant</title></head><body>",
    links,
    "</body></html>",
    sep = ""
)
writeLines(html, "index.html")
