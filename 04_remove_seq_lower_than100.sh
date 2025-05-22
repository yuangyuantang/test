#!/bin/bash

# Define the threshold for the number of sequences
threshold=100

# Create a text file to store the number of sequences
output_file="sequence_count.txt"

# Iterate over all the fastq files in the current directory
for file in *.fastq; do
    # Count the number of sequences in the fastq file
    num_sequences=$(($(wc -l < "$file") / 4))
    
    # Check if the number of sequences is lower than the threshold
    if [ $num_sequences -lt $threshold ]; then
        # Create the directory if it doesn't exist
        mkdir -p lower_100
        # Move the fastq file to the lower_100 directory
        mv "$file" lower_100/
        echo "Moved $file to lower_100 directory"
    fi
    
    # Output the number of sequences to the text file
    echo "$file: $num_sequences" >> "$output_file"
done

