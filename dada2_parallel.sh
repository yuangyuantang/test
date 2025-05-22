#!/bin/bash

# ================================
# Step 1: Organize files by sample ID
for r1_file in *_R1.fastq.gz *_R1.fq.gz; do
    [[ -e "$r1_file" ]] || continue
    sampleid=$(echo "$r1_file" | sed -E 's/_R1\.(fastq|fq)\.gz$//')
    mkdir -p "$sampleid"
    for ext in fastq fq; do
        [[ -e "${sampleid}_R1.${ext}.gz" ]] && mv "${sampleid}_R1.${ext}.gz" "$sampleid/"
        [[ -e "${sampleid}_R2.${ext}.gz" ]] && mv "${sampleid}_R2.${ext}.gz" "$sampleid/"
    done
done

# ================================
# Step 2: Define the per-sample pipeline
process_sample() {
    dir="$1"
    cd "$dir" || exit 1
    sampleid="${dir%/}"

    failed=false

    # 1: Strip barcodes
    # Generate files *__1.fastq, *__2.fastq, *_barcodes.fastq, *_scrap.fastq
    echo "[$sampleid] Running strip_addons.py"
    python /media/zhailab/Pipeline/yy_code/dada2_yy/strip_addons.py *_R1*.gz *_R2*.gz \
        -remove_bar_primer \
        -fw_primer GCATCGATGAAGAACGCAGC \
        -rev_primer TCCTCCGCTTATTGATATGC || failed=true
    mkdir -p mixed_package && mv *_R1*.gz *_R2*.gz mixed_package 2>/dev/null

    # 2: Check intermediates
    [[ -f *_barcodes.* && -f *_scrap.* && -f *__1.* && -f *__2.* ]] || failed=true

    # 3: Split samples by barcode
    # Generate individual samples with names "*-1_R1.fastq", "*-1_R2.fastq", "*-2_R1.fastq", "*-2_R2.fastq"
    [[ "$failed" == false ]] && Rscript /media/zhailab/Pipeline/yy_code/dada2_yy/03_ITS_Barcode_Lookup.R || failed=true

    # 4: Filter low-seq samples (<100 sequences)
    # DADA2 do not process samples with limited sequences
    [[ "$failed" == false ]] && bash /media/zhailab/Pipeline/yy_code/dada2_yy/04_remove_seq_lower_than100.sh || failed=true
    find . -maxdepth 1 -type f -name "*R*.fastq" -exec bgzip {} \;

    # 5: Prepare for DADA2
    # Run filterAndTrim to removes reads with ambiguous bases (maxN = 0) or very low quality 
    # Run cutadapt to remove primers
    # Generate prefiltered sequences stored in filtN and cutadapt
    [[ "$failed" == false ]] && Rscript /media/zhailab/Pipeline/yy_code/dada2_yy/05_preparation_for_dada2.R || failed=true

    # 6: DADA2 filtering
    # Generate filtered sequences with names "*_filter_1.fastq.gz", "*_filter_2.fastq.gz", "*_seqtab.rds", "*_trace_list.rds"
    [[ "$failed" == false ]] && Rscript /media/zhailab/Pipeline/yy_code/dada2_yy/process_dada2_EE8Q8.R || failed=true

    # 7: Infer ASVs
    # Generate line/ with fasta sequences of each sample
    # Generate files: "*.1280089.txt" files 
    [[ "$failed" == false ]] && Rscript /media/zhailab/Pipeline/yy_code/dada2_yy/dada2_pipeline_UNITE.R || failed=true

    # 8: Assign taxonomy
    # Exporting full and top hits of taxonomic results to CSV files
    # Generate files: "*.top.csv", "*.full.csv"
    [[ "$failed" == false ]] && Rscript /media/zhailab/Pipeline/yy_code/dada2_yy/assign_taxonomy_UNITE.hightaxa.R || failed=true

    # 9: Reformat top hits to QIIME2 style
    [[ "$failed" == false ]] && bash /media/zhailab/Pipeline/yy_code/dada2_yy/reformat_top.sh || failed=true

    cd - > /dev/null

    if $failed; then
        echo "[FAIL] $sampleid encountered an error in processing." >> ../step_failures.log
        return 1
    else
        return 0
    fi
}

export -f process_sample

# ================================
# Step 3: Run samples in parallel and log outputs per step
find . -maxdepth 1 -type d ! -name '.' | parallel -j 4 --joblog step_execution.log 'process_sample {} >> step_combined_output.log 2>&1'

# ================================
# Step 4: Collect all *.reformat.top.csv to a summary directory
mkdir -p summary
find . -type f -name "*.reformat.top.csv" -exec mv {} summary/ \;

# ================================
# Step 5: Run final Qiime2-style formatting in summary folder
cd summary
Rscript /media/zhailab/Pipeline/yy_code/dada2_yy/qiime2styleoutput.R > ../step10_qiime_summary.log 2>&1 || echo "[WARNING] Step 10 failed. Check step10_qiime_summary.log" >> ../step_failures.log
