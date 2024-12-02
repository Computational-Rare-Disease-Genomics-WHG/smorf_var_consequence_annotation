#!/bin/bash
# vep_results_filt 0.0.1

main() {
    dx-download-all-inputs

    #transcript_CDS_108.csv
    dx download project-GQ4zfYjJ55K2F7xkYjkXBX5F:file-Gb184F0J55K55VgJqxgv4Xy3
    #final_highconf.bed
    dx download project-GQ4zfYjJ55K2F7xkYjkXBX5F:file-Gb1845QJ55K55VgJqxgv4Xxv
    # R script
    dx download project-GQ4zfYjJ55K2F7xkYjkXBX5F:file-Gb2PzF0J55KBG9XZ9PgBzz05
    #CDS_exons.csv
    dx download project-GQ4zfYjJ55K2F7xkYjkXBX5F:file-Gb1845QJ55K82vxXbv66bK14
    #smorf_exon_bounds.csv
    dx download project-GQ4zfYjJ55K2F7xkYjkXBX5F:file-Gb1845QJ55KJykY0F4G4P503
    
    # requirements- VEP results, spliceAI results, info on transcript CDS location (transcript_CDS_108.csv, CDS_exons.csv), smORFs (smorfs_public_highconf.bed, smorf_exon_bounds.csv)
    #combined_data_norm_split_filt.txt  - spliceAI results for a subset of variants that aren't covered by precomputed scores

    /usr/bin/Rscript "vep_results_filt.R" \
    "/home/dnanexus/in/vep_results/${vep_results_name}" \
    "combined_data_norm_split_filt.txt" \ 
    "transcript_CDS_108.csv" \
    "smorfs_public_highconf.bed" \
    "CDS_exons.csv" \
    "smorf_exon_bounds.csv" \
    $chunk_index


    mv vep_results_filt.tsv.gz "vep_results_filt_chunk_${chunk_index}.tsv.gz"
    mkdir -p out/filtered_variants/
    mv "vep_results_filt_chunk_${chunk_index}.tsv.gz" out/filtered_variants/

    mv filtering_summary.tsv "filtering_summary_chunk_${chunk_index}.tsv"
    mkdir -p out/filtering_summary/
    mv "filtering_summary_chunk_${chunk_index}.tsv" out/filtering_summary/

    dx-upload-all-outputs
}