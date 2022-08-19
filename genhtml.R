files <- Sys.glob("**/*.html")
links <- paste(
   lapply(files, function (x) paste('<a href="/', x, '">', x, '</a>', sep="")),
   collapse="\n<br>\n"
)
html <- paste(
   "<html><head><title>Quant</title></head><body>",
   links,
   "</body></html>",
   sep=""
)
writeLines(html, "index.html")