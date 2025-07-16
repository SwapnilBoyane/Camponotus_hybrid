install.packages("pdftools")
install.packages("magick")
library(pdftools)
library(magick)
library(grid)
install.packages("qpdf")
library(qpdf)

# Define the directory where your PDF files are located
pdf_directory <- "/Volumes/camp_hybrid/pdf"
setwd("/Volumes/camp_hybrid/pdf")

# List all PDF files in the directory
pdf_files <- list.files(pdf_directory, pattern = "\\.pdf$", full.names = TRUE)


# Specify the output file name for the combined PDF
output_pdf <- "combined_output.pdf"

# Combine the PDFs
pdf_combine(input = pdf_files, output = output_pdf)

# Message to indicate completion
message("PDF files combined successfully into ", output_pdf)
