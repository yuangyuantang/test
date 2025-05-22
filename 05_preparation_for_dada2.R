# Load necessary libraries
library(ShortRead)
library(dada2)
library(stringr)
library(logr)

rm(list = ls())  # Clear environment

# Set working directory
path <- getwd()
path_parts <- unlist(str_split(path, "/"))
log_file <- file(paste0(path_parts[length(path_parts)], "_DADA_pt5.txt"))

# Log console output
sink(log_file, append = TRUE)
sink(log_file, append = TRUE, type = "message")

# Define primer sequences (change if needed)
FWD <- "GCATCGATGAAGAACGCAGC"
REV <- "TCCTCCGCTTATTGATATGC"

# List paired gzipped FASTQ files
R1_files <- list.files(path, pattern = "_R1\\.fastq\\.gz$", full.names = TRUE)
R2_files <- list.files(path, pattern = "_R2\\.fastq\\.gz$", full.names = TRUE)

# Function to generate all orientations of a primer
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna,
               Complement = complement(dna),
               Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  sapply(orients, toString)
}

# Function to count primer hits
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  sum(nhits > 0)
}

# Loop over each file pair
for (i in seq_along(R1_files)) {
  fnF <- R1_files[i]
  fnR <- R2_files[i]
  sample_id <- str_split(basename(fnF), "_")[[1]][1]

  # Create directory for filtered output
  filt_dir <- file.path(path, "filtN")
  dir.create(filt_dir, showWarnings = FALSE)
  fnF_filt <- file.path(filt_dir, basename(fnF))
  fnR_filt <- file.path(filt_dir, basename(fnR))

  # Removes reads with ambiguous bases (maxN = 0) or very low quality (output will also be .gz)
  filterAndTrim(fnF, fnF_filt, fnR, fnR_filt,
                maxN = 0, truncQ = 0, multithread = TRUE, compress = TRUE)

  # Count primer hits BEFORE cutadapt
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)

  cat("Primer hits BEFORE cutadapt for", sample_id, "\n")
  print(rbind(
    FWD.Fwd = sapply(FWD.orients, primerHits, fn = fnF_filt),
    FWD.Rev = sapply(FWD.orients, primerHits, fn = fnR_filt),
    REV.Fwd = sapply(REV.orients, primerHits, fn = fnF_filt),
    REV.Rev = sapply(REV.orients, primerHits, fn = fnR_filt)
  ))

  # Run cutadapt  externally from R to remove primers from the filtered reads
  # Finds and trims primers in both read
  # Supports paired-end mode
  # Saves output in a new cutadapt/ folder
  cutadapt <- "/usr/bin/cutadapt"  # <-- Update path if needed
  cut_dir <- file.path(path, "cutadapt")
  dir.create(cut_dir, showWarnings = FALSE)
  fnF_cut <- file.path(cut_dir, basename(fnF))
  fnR_cut <- file.path(cut_dir, basename(fnR))

  R1_flags <- paste("-g", FWD, "-a", dada2:::rc(REV))
  R2_flags <- paste("-G", REV, "-A", dada2:::rc(FWD))
  args <- c(R1_flags, R2_flags, "-n", 2, "-o", fnF_cut, "-p", fnR_cut,
            fnF_filt, fnR_filt, "-j", "0")

  out <- system2(cutadapt, args = args, stdout = TRUE)
  cat(out, file = log_file, sep = "\n")

  # Count primer hits AFTER cutadapt
  cat("Primer hits AFTER cutadapt for", sample_id, "\n")
  print(rbind(
    FWD.Fwd = sapply(FWD.orients, primerHits, fn = fnF_cut),
    FWD.Rev = sapply(FWD.orients, primerHits, fn = fnR_cut),
    REV.Fwd = sapply(REV.orients, primerHits, fn = fnF_cut),
    REV.Rev = sapply(REV.orients, primerHits, fn = fnR_cut)
  ))
}

# Close logging
sink()
sink(type = "message")
closeAllConnections()
