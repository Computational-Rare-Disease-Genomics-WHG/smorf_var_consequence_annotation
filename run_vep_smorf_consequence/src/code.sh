#!/bin/bash
# vep_500k_wgs_111 0.0.1

main() {
    dx-download-all-inputs

    # download docker image
    docker load -i "$docker_image_path"
    chmod -R a+rwx $HOME

    # mount project folder - for Cache
    FUSE_MOUNT=$HOME/projects
    DX_PROJECT_CONTEXT_ID='project-XXXXXXXXXXXXXXXXXXXXXXXX'
    mkdir -p $FUSE_MOUNT
    sudo -E dxfuse -uid $(id -u) -gid $(id -g) -verbose 2 $FUSE_MOUNT $DX_PROJECT_CONTEXT_ID

    # prepare file and run docker
    output_name="${input_vcf_name}_vep.txt"
    gunzip -dc $input_vcf_path | \
    docker run -i --rm \
        -v $FUSE_MOUNT/500k_WGS/VEP:/opt/vep/.vep\
        -v $HOME:/vep_out:Z\
        ensemblorg/ensembl-vep:release_108.2 \
        ./vep \
        --assembly GRCh38\
        --force_overwrite\
        --format vcf\
        --species homo_sapiens\
        --cache --fork 4 --offline\
        --mane \
        --mane_select\
        --minimal\
        --canonical\
        --tab \
        --fasta /vep_out/in/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa\
        --dir_cache /opt/vep/.vep/Cache\
        --output_file /vep_out/$output_name\
        --numbers \
        --af_gnomadg \
        --plugin UTRAnnotator \
        --plugin SpliceAI,snv=/opt/vep/.vep/SpliceAI/spliceai_scores.masked.snv.hg38.vcf.gz,indel=/opt/vep/.vep/SpliceAI/spliceai_scores.masked.indel.hg38.vcf.gz \
        --custom /opt/vep/.vep/PhyloP/hg38.phyloP100way.bw,PhyloP,bigwig \
        --gtf /opt/vep/.vep/$smorf_gtf_name
    
    mkdir -p out/filtered_vep
    mv $output_name out/filtered_vep/
    dx-upload-all-outputs
}