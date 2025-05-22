library(stringr)

# 1. Set working directory and path
# setwd("/media/zw/DADA2/miseq0427Mock/V45")
# setwd("/media/zw/DADA2/miseq0427/ComparingITS")
# setwd("/media/zw/DADA2/NH/S16/V45/EE7Q7")
path <- getwd()
path1 <- unlist(str_split(path, "/"))

#python3 RemoveLine.py /media/zw/DADA2/K4NKV/B_EEQ/test
name1 <- paste0(path, "/line") # Path to a folder ONLY contains the .fasta files to remove lines
out <- system(paste0("python /media/zhailab/Mix/DADA2/Script/RemoveLine.py ", name1), intern = TRUE) # Remember to write the correct path to RemoveLine.py

# 2. Save all the type A and type B file names as R1 and R2
R1 <- list.files(pattern = "*.reformat.top.csv")
if (identical(R1, character()))
{
  R1 <- list.files(pattern = "*.silva138_1.top.csv")
}

if (identical(R1, character()))
{
  R1 <- list.files(pattern = "*.nt.top.csv")
}

R2 <- list.files(pattern = "line.*.fasta")


#unequal length
db <-  unlist(strsplit(R1[1], split = "\\."))
dbtail <- paste0(".", db[2], ".top.csv")
T1 <- gsub(dbtail,"",R1)

T2 <- gsub("\\line.","",R2)
T2 <- gsub("\\.fasta","",T2)

p <- which(!(T2 %in% T1))
sdiff <- setdiff(T2, T1)
print (sdiff)

# 3. Create sw
sw <- c("asv", "species")
sw <- c(sw, T1)
# for (z in 1:length(R2))
#   {
#   name <- unlist(strsplit(R2[z],split = "line."))
#   newname <- unlist(strsplit(name[2],split = ".fasta"))
#   sw <- c(sw, newname)
# }

# Check if R1 and R2 have the same length
if (length(R1) != length(R2)) {
  warning("R1 and R2 do not have the same length.")
  # break
  for (l in 1:(length(p)))
  {
    sw <- c(sw, sdiff[l])
  }
}
dim(sw) <- c(1,length(R2) + 2)


# 4. Loop through R1 and R2

for (x in 1:length(R1))
  {
  if (!(T1[x] %in% sdiff))
  {
    # print(T1[x])
  # 4a. For each line in R2, check for ">"
  fasta_lines <- readLines(paste0(path, "/line.", T1[x], ".fasta"))
  placeholder <- ""

  for (i in 1:length(fasta_lines))
    {
    if (startsWith(fasta_lines[i], ">"))
      {
      placeholder <- unlist(strsplit(fasta_lines[i],split = ">"))
      }
    else
      {
        if (fasta_lines[i] %in% sw[, 1])
        {
          pos <- which(sw==fasta_lines[i], arr.ind = TRUE)
          sw[pos[1], x+2] <- placeholder[2]

        }
        else
        {
          sw <- rbind(sw, rep(NA, ncol(sw)))
          sw[nrow(sw), 1] <- fasta_lines[i]
          sw[nrow(sw), x+2] <- placeholder[2]
        }
      }
  }

  # 4b. For each row in R1, check for matching "otu"
  csv_data <- read.csv(paste0(path, "/", R1[x]), header = TRUE)
  
  for (i in 1:nrow(csv_data))
  {
    if (toString(csv_data$otu[i]) %in% sw[,x+2])
      {
        pos <- which(sw==toString(csv_data$otu[i]), arr.ind = TRUE)
        sw[pos[1], 2] <- toString(csv_data$Species[i])
      }
  }
  # 4c. Remove extra info
  sw[, x+2] <- gsub(".*;seqs=","",sw[, x+2])
  sw[, x+2] <- gsub(";samples=.*","",sw[, x+2])
  }
}

num = 0
for (y in p)
{
  num =+ 1
  # 4a. For each line in R2, check for ">"
  fasta_lines <- readLines(paste0(path, "/", R2[y]))
  placeholder <- ""
  
  for (i in 1:length(fasta_lines))
  {
    if (startsWith(fasta_lines[i], ">"))
    {
      placeholder <- unlist(strsplit(fasta_lines[i],split = ">"))
    }
    else
    {
      if (fasta_lines[i] %in% sw[, 1])
      {
        pos <- which(sw==fasta_lines[i], arr.ind = TRUE)
        sw[pos[1], length(T1)+2+num] <- placeholder[2]
        
      }
      else
      {
        sw <- rbind(sw, rep(NA, ncol(sw)))
        sw[nrow(sw), 1] <- fasta_lines[i]
        sw[nrow(sw), length(T1)+2+num] <- placeholder[2]
      }
    }
  }
  # 4c. Remove extra info
  sw[, length(T1)+2+num] <- gsub(".*;seqs=","",sw[, length(T1)+2+num])
  sw[, length(T1)+2+num] <- gsub(";samples=.*","",sw[, length(T1)+2+num])
}

# 6. Save "sw" as "DNF_Qiime2_style.csv"
write.csv(sw, file = paste0(path1[length(path1)], "_Qiime2_style.csv"), row.names = FALSE)



###################################################################################################################################################


# # 1. Set working directory and path
# setwd("/media/zw/DADA2/miseq0427")
# path <- getwd()
# path1 <- unlist(str_split(path, "/"))
# 
# 
# # 3. Create output_df
# H1 <- read.csv("H1.csv")
# H2 <- read.csv("H2.csv")
# H3 <- read.csv("H3.csv")
# 
# # Create an empty output data frame with the same column names as the input data frame
# # output_df <- data.frame(matrix(ncol = ncol(input_df), nrow = 0))
# output_df <- data.frame(matrix(ncol = 65, nrow = 0))
# hn <- c(colnames(H1), colnames(H2), colnames(H3))
# hn <- unique(hn)
# colnames(output_df) <- hn
# 
# input_df <- H3
# 
# # Loop through each row of the input dataframe
# for (i in 1:(nrow(input_df))) {
#   
#   # Check if ID is unique in output dataframe
#   if (input_df[i,"asv"] %in% output_df$asv) {
#     
#     # If ID is not unique, find the row in the output dataframe with the same ID
#     row_index <- which(output_df$asv == as.vector(input_df[i,"asv"]))
#     
#     # Record the later column values to the columns in the output dataframe that has the same column name
#     for (k in 1: length(colnames(input_df)))
#     {
#       # If ID is unique, add a new row to the output dataframe
#       output_df[row_index, colnames(input_df)[k]] <- as.vector(input_df[i, colnames(input_df)[k]])
#     }
#     
#   } else {
#     
#     empty <- data.frame(matrix(NA, ncol = ncol(output_df)))
#     colnames(empty) <- colnames(output_df)
#     output_df <- rbind(output_df, empty)
#     # temp <- c(as.vector(input_df[i, 1]), as.vector(input_df[i, 2]))
#     
#     for (k in 1: length(colnames(input_df)))
#     {
#       # If ID is unique, add a new row to the output dataframe
#       output_df[nrow(output_df), colnames(input_df)[k]] <- as.vector(input_df[nrow(output_df), colnames(input_df)[k]])
#     }
#   }
#   
# }
# 
# write.csv(output_df, "Healthy_samples.csv")


###################################################################################################################################################
