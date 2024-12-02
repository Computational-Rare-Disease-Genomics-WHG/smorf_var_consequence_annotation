# smorf_var_consequence_annotation
 A set of scripts to allow annotation of variant consequence on a custom set of smORFs

Right now this is set up with a series of DNANexus app as it was developed to find variant consequences on UK Biobank variants, however it's as applicable to any variant set and is basically just a series of R and bash scripts. 

Required files
- smORF set in bed format. A smORF set generated from public datasets is made available in data/smorfs_public_highconf.bed
- [Homo_sapiens.GRCh38.108.chr.gtf.gz](https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz) (or whatever annotation GTF was used to make the smORF set)

Steps in pipeline
1. Run Rscript convert_smorfs_gtf_format.R which converts smORF bed format to GTF format (or alternatively start with a smORF set represented in GTF format)

2. Run DNAnexus app run_vep_smorf_consequence
- run on VCF and requires smORF GTF file generated in step 1
- also requires docker image of VEP on RAP, as well as VEP Cache and fasta file
- If you have your own alternate VEP setup, and have generated a smORF GTF file, the only vital line you need to add to your VEP command is:

```
--gtf /opt/vep/.vep/$smorf_gtf_name
```

This will ensure VEP annotates the variant consequence on the smORF CDS, rather than the normal transcript CDS. 

3. Run Rscript vep_filter_helper_files.R
- This generates the files required for the next DNAnexus app. Requires reference GTF and smORF bed file as above.

4. Run DNAnexus app vep_results_filt
- This is basically just a wrapper around another R script (resources/home/dnanexus/vep_results_filt.R), which takes VEP output from step 2, and (somewhat optionally) spliceAI output, as well as the helper files from step 3. It then does a bunch of processing on the VEP output to generate a final set of annotated consequences on smORF variants. It also does some monitoring of how many variants were annotated in different categories.  