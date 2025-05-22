# ===============================
# 1. Load Required Libraries
# ===============================
library(dada2)        # For ASV handling
library(phyloseq)     # For phylogenetic data integration
library(Biostrings)   # For DNA sequence I/O
library(dplyr)        # Data wrangling
library(stringr)      # String operations
library(tidyr)        # Data reshaping
library(gdata)        # General data tools

# ===============================
# 2. Define Working Directory
# ===============================
path <- getwd()  # Get current directory

# ===============================
# 3. Define BLAST Parsing Function
# ===============================
qblast <- function(tax.file = "asv_seqs.fasta.unite.txt", tax_table = TRUE) {
  # Load namespaces if not already
  requireNamespace("data.table", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("stringr", quietly = TRUE)

  # Load BLAST result as data frame
  t <- data.table::fread(tax.file, quote = "") %>% tibble::as_tibble()

  # Define BLAST column names (outfmt 6 with custom format)
  colnames(t) <- c("qseqid", "qlen", "sseqid", "slen", "qstart", "qend",
                   "sstart", "send", "nident", "mismatch", "gapopen", "gaps",
                   "ppos", "frames", "qframe", "sframe", "qcovs", "qcovhsp",
                   "evalue", "bitscore", "score", "length", "pident")

  # Separate taxonomy fields
  t <- t %>%
    separate(sseqid, into = c("accession_unite", "taxonomy", "species_hypothesis"), sep = "[|]") %>%
    separate(taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "[;]", remove = FALSE)

  # Clean taxonomy labels
  t <- t %>%
    mutate(
      Kingdom = str_remove(Kingdom, "k__"),
      Phylum = str_remove(Phylum, "p__"),
      Class = str_remove(Class, "c__"),
      Order = str_remove(Order, "o__"),
      Family = str_remove(Family, "f__"),
      Genus = str_remove(Genus, "g__"),
      Species = str_remove(Species, "s__"),
      Species = str_replace(Species, "_", " ")
    )

  # Extract ASV (OTU) number
  t <- t %>%
    mutate(
      otu = qseqid,
      otu.number = as.numeric(str_extract(otu, "(?<=ASV_)[0-9]+"))
    )

  # Count how often each species appears per OTU
  t <- t %>%
    group_by(otu, Species) %>%
    mutate(occurence.sp = n()) %>%
    ungroup()

  # Filter out ambiguous genus-level identifications
  t <- t %>%
    mutate(genuslevel = grepl("_sp|unidentified|uncultured ", Species)) %>%
    filter(genuslevel == FALSE)

  # Choose best match per OTU
  t <- t %>%
    group_by(otu) %>%
    arrange(desc(pident), desc(occurence.sp), evalue) %>%
    filter(!duplicated(taxonomy)) %>%
    select(otu, Phylum, Family, Species, evalue, occurence.sp, pident, length, everything())

  # Return full table or summarized tax_table
  if (!tax_table) {
    return(t %>% ungroup() %>% arrange(otu.number))
  } else {
    return(t %>%
             filter(row_number() == 1) %>%
             ungroup() %>%
             arrange(otu.number) %>%
             select(otu, evalue, pident, Kingdom, Phylum, Class, Order, Family, Genus, Species))
  }
}

# ===============================
# 4. Locate and Process All BLAST Files
# ===============================
# Get all relevant BLAST result files
T1 <- sort(list.files(path, pattern = "*1280089.txt$", full.names = TRUE))  # Adjust pattern if using different database

# Loop through each BLAST result file
for (t in seq_along(T1)) {
  if (file.info(T1[t])$size != 0) {  # Skip empty files

    # Generate full annotation table (all matches)
    tblast.full <- qblast(T1[t], tax_table = FALSE)

    # Generate top match only (per ASV)
    tblast.top <- qblast(T1[t], tax_table = TRUE)

    # Extract filename for saving output
    samplename <- unlist(str_split(T1[t], "/"))
    name <- unlist(str_split(samplename[length(samplename)], ".txt"))

    # Write both result types to CSV
    write.csv(tblast.top,  paste0(name[1], ".top.csv"),  row.names = FALSE)
    write.csv(tblast.full, paste0(name[1], ".full.csv"), row.names = FALSE)
  }
}
