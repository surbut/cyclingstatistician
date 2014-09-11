# Set working directory
setwd("/Users/sarahurbut/Dropbox/Thesis!")

# Load packages
require(knitr)
require(markdown)

# Create slides
knit("My_Analysis.Rmd")
system("pandoc -s -t slidy My_Analysis.md -o My_Analysis.html")
"/Users/sarahurbut/Dropbox/Thesis!"
setwd("C:/Documents and Settings/name")

# Load packages
require(knitr)
require(markdown)

# Create .md, .html, and .pdf files
knit("My_Analysis.Rmd")
markdownToHTML('My_Analysis.md', 'My_Analysis.html', options=c("use_xhml"))
system("pandoc -s My_Analysis.html -o My_Analysis.pdf")

###From terminal

pandoc -s --mathjax -o output.html input.Rmd
