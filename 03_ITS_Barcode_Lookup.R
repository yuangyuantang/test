# Load required libraries
library(ShortRead)
library(dada2)
library(stringr)
library(logr)

# Set working directory
path <- getwd()

# Detect a single __1.fastq file (you can also set this manually)
r1_file_path <- list.files(path, pattern = "*__1.fastq$", full.names = TRUE)[1]

# Stop if no file is found
if (is.na(r1_file_path)) {
  stop("No __1.fastq file found in the directory.")
}

# Extract base name (without directory and suffix)
file_parts <- unlist(str_split(r1_file_path, "/"))
base_name <- unlist(str_split(file_parts[length(file_parts)], "__1.fastq"))[1]

# Define file paths
barcode_file <- file(description = paste0(base_name, "_barcodes.fastq"), open = "r", blocking = TRUE)
r1_file      <- file(description = paste0(base_name, "__1.fastq"), open = "r", blocking = TRUE)
r2_file      <- file(description = paste0(base_name, "__2.fastq"), open = "r", blocking = TRUE)

# Read contents of files
barcode_lines <- readLines(barcode_file)
r1_lines      <- readLines(r1_file)
r2_lines      <- readLines(r2_file)

# Load and clean barcode reference
barcode_df <- read.csv("~/ITS_Barcode_Mix.csv", header = FALSE)
barcode_matrix <- gsub(" ", "", as.matrix(barcode_df))

# Process each barcode
for (l in 1:nrow(barcode_matrix)) {
  
  barcode_seq <- barcode_matrix[l, 2]
  match_indices <- grep(barcode_seq, barcode_lines)

  for (n in match_indices) {
    if (n > 1 && n + 2 <= length(r1_lines) && n + 2 <= length(r2_lines)) {
      
      r1_block <- r1_lines[(n - 1):(n + 2)]
      r2_block <- r2_lines[(n - 1):(n + 2)]
      
      sample_name <- unlist(str_split(base_name, "_"))[1]
      out_r1 <- paste0(sample_name, "-", barcode_matrix[l, 1], "_R1.fastq")
      out_r2 <- paste0(sample_name, "-", barcode_matrix[l, 1], "_R2.fastq")
      
      write(r1_block, out_r1, append = TRUE, sep = "\n")
      write(r2_block, out_r2, append = TRUE, sep = "\n")
    }
  }
}

# Close files
closeAllConnections()
