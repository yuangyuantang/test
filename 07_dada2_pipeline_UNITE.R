# ================================
#        Load Required Libraries
# ================================
# dada2: Denoising and ASV inference
# phyloseq: For creating phylogenetic data structures
# Biostrings: For handling DNA sequences
# dplyr, stringr, tidyr: Tidyverse utilities for string manipulation and data wrangling

library(dada2)
library(phyloseq)
library(Biostrings)
library(dplyr)
library(stringr)
library(tidyr)


# ===================================
#        Initialize Logging Setup
# ===================================
# Set the working directory path
path <- getwd()

# Split path into components to extract the final folder name
path_parts <- str_split(path, "/")[[1]]

# Define log file using the last folder name + suffix
log_file <- paste0(tail(path_parts, 1), "_DADA_pt7.txt")

# Redirect console output to the log file (for debugging and reproducibility)
sink(log_file, append = TRUE)
sink(log_file, append = TRUE, type = "output")


# ========================================
#       Identify seqtab.rds Input Files
# ========================================
# Look for all files ending with _seqtab.rds in current working directory
seqtab_files <- sort(list.files(path, pattern = "_seqtab\\.rds$", full.names = TRUE))

# Initialize a list to keep track of successfully processed files
valid_indices <- c()


# ====================================
#     Create Output Directory ("line")
# ====================================
# This directory will store the .fasta files generated for each sample
line_dir <- file.path(path, "line")

# Create the directory if it doesn't exist
if (!dir.exists(line_dir)) dir.create(line_dir)


# ====================================
#     Process Each seqtab.rds File
# ====================================
# Loop over each RDS file to extract and export DNA sequences

for (i in seq_along(seqtab_files)) {
  seqtab_path <- seqtab_files[i]
  filename <- basename(seqtab_path)  # Extract file name
  sample_name <- str_remove(filename, "_seqtab\\.rds$")  # Extract sample ID
  output_fasta <- file.path(line_dir, paste0(sample_name, ".fasta"))  # Set output fasta path

  # Error handling in case the file somehow doesn't exist
  if (!file.exists(seqtab_path)) {
    stop("YTError: Missing file: ", seqtab_path)
  }

  # Load the sequence table from RDS
  seqtab <- readRDS(seqtab_path)

  # Skip the file if it contains no ASVs (zero columns)
  if (ncol(seqtab) == 0) {
    message("WARNING: File ", filename, " has zero columns. Skipping.")
    next
  }

  # Construct OTU table and DNAStringSet
  otu <- otu_table(t(seqtab), taxa_are_rows = TRUE)  # Transpose to have taxa as rows
  dna <- DNAStringSet(getSequences(seqtab))          # Extract DNA sequences
  names(dna) <- dna                                   # Use sequences as their own names

  # Build phyloseq object
  ps <- phyloseq(otu, dna)

  # Assign informative names to taxa: ASV number, total reads (seqs), and number of samples detected in
  taxa_names(ps) <- paste0("ASV_", seq_along(taxa_names(ps)), 
                           ";seqs=", taxa_sums(ps), 
                           ";samples=", apply(otu, 2, function(x) sum(x > 0)))

  # Write sequences to fasta file for BLAST analysis
  writeXStringSet(refseq(ps), output_fasta)

  # Keep track of successfully processed samples
  valid_indices <- c(valid_indices, i)

  message("Processed: ", filename)
}


# ============================
#   Run BLAST on FASTA Files
# ============================
# BLAST each sample's fasta file against a custom reference database

blast_db <- "/media/zhailab/Mix/DADA2/Database/1280089"  # Path to local BLAST database

for (i in valid_indices) {
  seqtab_path <- seqtab_files[i]
  filename <- basename(seqtab_path)
  sample_name <- str_remove(filename, "_seqtab\\.rds$")  # Extract sample ID
  fasta_path <- file.path(line_dir, paste0(sample_name, ".fasta"))  # Path to input fasta
  blast_out <- paste0(sample_name, ".1280089.txt")  # Define output file name for BLAST results

  # Construct the BLASTN command
  blast_cmd <- paste0(
    "blastn -query ", fasta_path,
    " -db ", blast_db,
    " -evalue 10 -max_target_seqs 50",
    " -out ", blast_out,
    ' -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send nident mismatch gapopen gaps ppos frames qframe sframe qcovs qcovhsp evalue bitscore score length pident"'
  )

  # Execute BLAST command in system shell
  system(blast_cmd)

  message("BLAST completed for: ", sample_name)
}
