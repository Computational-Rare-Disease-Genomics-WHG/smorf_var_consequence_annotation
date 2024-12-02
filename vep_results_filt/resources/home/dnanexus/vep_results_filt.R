install.packages(c('data.table', 'R.utils', 'janitor', 'dplyr', 'XML', 'BiocManager'), repos = "http://cran.us.r-project.org")
library(data.table)
library(janitor)
library(dplyr)
BiocManager::install("rtracklayer")

args = commandArgs(trailingOnly=TRUE)
vep = fread(args[1])
spliceai = fread(args[2])
tx_cds = fread(args[3])
smorfs = fread(args[4])
cds_exons = fread(args[5])
smorf_exon_bounds = fread(args[6])
chunk_index = as.numeric(args[7])






# making ref and alt allele columns 
setnames(vep, '#Uploaded_variation', 'variant_id')
vep[, ref := sapply(strsplit(variant_id, '-'), '[[', 3)]
vep[, alt := sapply(strsplit(variant_id, '-'), '[[', 4)]



# we're not interested in these consequences / mapping to nearby genes
vep = vep[!Consequence %in% c('upstream_gene_variant', 'downstream_gene_variant')]
smorf_vars = vep[SOURCE == 'final_highconf.gtf.gz']
vep = vep[SOURCE != 'final_highconf.gtf.gz']






# processing manual spliceai scores
spliceai[, SpliceAI_pred_manual := gsub('^[A-Z]*\\|', '', SpliceAI_pred_manual)]
spliceai[, Feature := sapply(strsplit(SpliceAI_pred_manual, '|', fixed = TRUE), '[[', 1)]
spliceai[, Location := paste0(CHROM, ':', POS)]

vep = spliceai[, .(Location, Feature, Allele = ALT, SpliceAI_pred_manual)][vep, on = .(Location, Feature, Allele)]

vep[!is.na(SpliceAI_pred_manual), `:=` (SpliceAI_pred_type = 'manual', SpliceAI_pred = SpliceAI_pred_manual)]
vep[is.na(SpliceAI_pred_manual), `:=` (SpliceAI_pred_type = 'precomputed')]







# processing spliceAI output
SpliceAI_pred_list = strsplit(vep$SpliceAI_pred, split = '|', fixed = TRUE)
SpliceAI_pred_dt = as.data.table(do.call(rbind, SpliceAI_pred_list))
SpliceAI_pred_dt = suppressWarnings(cbind(SpliceAI_pred_dt[, c(1)], SpliceAI_pred_dt[, c(2:9)][, lapply(.SD, as.numeric)]))

names(SpliceAI_pred_dt) = c('Gene', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL')
SpliceAI_pred_dt[!is.na(DS_AG), spliceai_max_delta := apply(SpliceAI_pred_dt[!is.na(DS_AG), 
                                                                             .(DS_AG, DS_AL, DS_DG, DS_DL)], 1, max, na.rm=TRUE)]

vep = cbind(vep, SpliceAI_pred_dt)
rm(SpliceAI_pred_list)

vep[is.na(spliceai_max_delta), spliceai_max_delta := 0]







# add spliceAI scores to smORF VEP results
spliceai_scores = vep[SpliceAI_pred != '-', .(variant_id, transcript_id = Feature, SpliceAI_pred_type, SpliceAI_pred, DS_AG, DS_AL, DS_DG, DS_DL, 
                                              DP_AG, DP_AL, DP_DG, DP_DL, spliceai_max_delta)]

smorf_vars[, SpliceAI_pred := NULL]
smorf_vars[, `:=` (smorf_id = Feature, transcript_id = sapply(strsplit(Feature, '_'), '[[', 1))]

smorf_vars = spliceai_scores[smorf_vars, on = .(variant_id, transcript_id)]

#add on smORF bounds
smorfs = smorfs[name == cluster]
smorf_cds = smorfs[, .(Feature = newname, smorf_start = chromStart, smorf_end = chromEnd)]

smorf_vars[, start := as.numeric(sapply(strsplit(Location, ':|-'), '[[', 2))]
smorf_vars[grepl('-', Location), end :=  as.numeric(sapply(strsplit(Location, ':|-'), '[[', 3))]
smorf_vars[!grepl('-', Location), end :=  start]

smorf_vars = smorf_cds[smorf_vars, on  = .(Feature)]







# get associated exon number for intronic splicing variants
smorf_vars[INTRON != '-', intron_no := as.numeric(sapply(strsplit(INTRON, '/'), '[[', 1))]
smorf_vars[grepl('acceptor', Consequence), exon_number := intron_no + 1]
smorf_vars[grepl('donor', Consequence), exon_number := intron_no ]
smorf_vars[EXON != '-', exon_number := as.numeric(sapply(strsplit(EXON, '/'), '[[', 1))]

# add on exon number containing cds start and end in transcript
smorf_vars = cds_exons[smorf_vars, on = .(transcript_id)]


# add exon that smORF starts and ends in
smorf_vars = smorf_exon_bounds[smorf_vars, on = .(Feature, transcript_id)]



# removing variants outside smORF boundaries, EXCEPT for splicing variants on exons which contain the smORF and not any annotated CDS for that transcript
variants_outside_smorf_bounds = smorf_vars[(end < smorf_start | start > smorf_end)]
variants_outside_smorf_bounds_keep = variants_outside_smorf_bounds[grepl('splice', Consequence)& 
                                                                     (exon_number < cds_start | exon_number > cds_end) &
                                                                     (exon_number >= smorf_start_exon & exon_number <= smorf_end_exon)]

smorf_vars = rbind(smorf_vars[!(end < smorf_start | start > smorf_end)],
                   variants_outside_smorf_bounds_keep)

# some indels are still here because i didn't left-align before running VEP- i should do that (note to self)
smorf_vars = smorf_vars[!Consequence %in% c('3_prime_UTR_variant', '5_prime_UTR_variant')]

# remove variants from vep that i've removed from smorf consequences
vep = vep[variant_id %in% smorf_vars$variant_id]
# remove consequences on NMD transcripts
vep = vep[!grepl('NMD_transcript_variant', Consequence)]







# downgrade impact of UTR donor/acceptor variants to low on VEP transcript so they're not automatically filtered out
# will probably have to do some more careful looking at these, some will be on exon containing start/stop codon. maybe include only splice donor/acceptor variants on internal UTR exons? so if skipping happens (most likely) CDS won't be affected?


# find any moderate/high impact outside CDS bounds for VEP annotated transcripts
vep = tx_cds[, .(Feature = transcript, cds_start, cds_end)][vep, on = .(Feature)]

vep[, var_pos := sapply(strsplit(Location, ':'), '[[', 2)]
vep[grepl('-', Location), var_start := as.numeric(sapply(strsplit(var_pos, '-'), '[[', 1))]
vep[grepl('-', Location), var_end := as.numeric(sapply(strsplit(var_pos, '-'), '[[', 2))]
vep[!grepl('-', Location), `:=` (var_start = as.numeric(var_pos), var_end = as.numeric(var_pos))]

utr_splicing_vars = vep[IMPACT %in% c('HIGH', 'MODERATE')&(var_start < cds_start | var_end > cds_end)]

utr_splicing_vars[, intron_no := as.numeric(sapply(strsplit(INTRON, '/'), '[[', 1))]
utr_splicing_vars[grepl('acceptor', Consequence), exon_number := intron_no + 1]
utr_splicing_vars[grepl('donor', Consequence), exon_number := intron_no ]
utr_splicing_vars[EXON != '-', exon_number := as.numeric(sapply(strsplit(EXON, '/'), '[[', 1))]

# add on cds star tand end exons
utr_splicing_vars = cds_exons[, .(Feature = transcript_id, 
                                  cds_start_exon = cds_start, 
                                  cds_end_exon = cds_end)][utr_splicing_vars, 
                                                           on = .(Feature)]

# these variants are donor/acceptor variants on exons where CDS starts or ends- keep these as HIGH impact
utrsplice_keephigh = utr_splicing_vars[exon_number >= cds_start_exon & exon_number <= cds_end_exon]


# downgrade impact of UTR donor/acceptor variants to low on VEP transcript so they're not automatically filtered out
vep[IMPACT %in% c('HIGH', 'MODERATE')&(var_start < cds_start | var_end > cds_end) & !variant_id %in% utrsplice_keephigh, IMPACT := 'LOW']








# initialising data.table that keeps track of how many variants are filtered at different steps, # of variant categories per file etc
this_file_stats = data.table(filename = chunk_index, variant_count = length(unique(vep$variant_id)))




this_file_stats[, spliceai_manual := nrow(unique(vep[SpliceAI_pred_type == 'manual', 
                                                     .(variant_id, Allele)]))]






### step 1: removing variants with high moderate or predicted splicing impact on annotated CDS
# get variant IDs for variants we want to exclude based on vep consequences in any transcript
variants_to_exclude = c()

# remove high and moderate impact variants except in UTR or non-coding transcripts
variants_to_exclude = unique(c(variants_to_exclude, 
                               unique(vep[IMPACT %in% c('HIGH', 'MODERATE') & 
                                            ! (grepl('(UTR_variant|non_coding_transcript)', Consequence) & grepl('splice', Consequence)), variant_id])))






this_file_stats[, remove_cds_highmoderate := -length(variants_to_exclude)]

# remove low impact variants withing CDS start and end with spliceAI delta > 0.2
variants_to_exclude = unique(c(variants_to_exclude, 
                               vep[IMPACT == 'LOW' & spliceai_max_delta >= 0.2 & 
                                     var_pos >= cds_start & var_pos <= cds_end &
                                     !is.na(cds_start) & !is.na(cds_end)  & 
                                     ! (grepl('(UTR_variant|non_coding_transcript)', Consequence) & grepl('splice', Consequence)), variant_id]))

this_file_stats[, remove_cds_lowimpactspliceai := -(length(variants_to_exclude) + remove_cds_highmoderate)]

# remove intronic variants within CDS start and end with spliceAI delta > 0.2
vep = tx_cds[, .(Feature = transcript, cds_start, cds_end, tx = 1)][vep, on = .(Feature)]
vep[, var_pos := sapply(strsplit(Location, ':', fixed = TRUE), '[[', 2)]

variants_to_exclude = unique(c(variants_to_exclude, 
                               vep[Consequence == 'intron_variant' & 
                                     spliceai_max_delta >= 0.2 &  
                                     var_pos >= cds_start & var_pos <= cds_end &
                                     !is.na(cds_start) & !is.na(cds_end) &
                                     ! (grepl('(UTR_variant|non_coding_transcript)', Consequence) & grepl('splice', Consequence)), variant_id]))

this_file_stats[, remove_cds_intronicspliceai := -(length(variants_to_exclude) + remove_cds_highmoderate + remove_cds_lowimpactspliceai)]
this_file_stats[, remove_total := length(variants_to_exclude)]


# caveat, if it's in 'processed transcript' just flag instead of excluding for now
variants_to_flag = unique(smorf_vars[variant_id %in% variants_to_exclude & Feature %in% smorfs[type == 'processed_transcript']$newname, variant_id])

this_file_stats[, flag_processed_transcript := +length(variants_to_flag)]







smorf_vars[variant_id %in% variants_to_flag, processed_tx_flag := 1]
smorf_vars[is.na(processed_tx_flag), processed_tx_flag := 0]

vep[variant_id %in% variants_to_flag, processed_tx_flag := 1]
vep[is.na(processed_tx_flag), processed_tx_flag := 0]

# exclude variants from smORF-EP and VEP output
smorf_vars = smorf_vars[!variant_id %in% setdiff(variants_to_exclude, variants_to_flag)]
vep = vep[!variant_id %in% setdiff(variants_to_exclude, variants_to_flag)]

# add UTR variant, non coding transcript variant flag
vep[grepl('UTR_variant', Consequence), UTR_variant_flag := 1]
vep[is.na(UTR_variant_flag), UTR_variant_flag := 0]

vep[grepl('non_coding_transcript', Consequence), noncoding_transcript_variant_flag := 1]
vep[is.na(noncoding_transcript_variant_flag), noncoding_transcript_variant_flag := 0]

this_file_stats[, remainingvars := length(unique(smorf_vars$variant_id))]
this_file_stats[, remainingvars_cds := nrow(unique(vep[noncoding_transcript_variant_flag == 0 & UTR_variant_flag == 0, .(variant_id, Allele)]))]
this_file_stats[, remainingvars_noncoding := nrow(unique(vep[noncoding_transcript_variant_flag == 1, .(variant_id, Allele)]))]
this_file_stats[, remainingvars_utr := nrow(unique(vep[UTR_variant_flag == 1, .(variant_id, Allele)]))]




# step 3: match VEP and smORF-EP consequences by transcript and do some further filtering
# join vep and smorf consequences by transcript
smorf_vars[, `:=` (smorf_id = Feature, transcript_id = sapply(strsplit(Feature, '_'), '[[', 1))]

vep_join = vep[Feature %in% smorf_vars$transcript_id, .(variant_id, transcript_id = Feature, tx_exon = EXON, tx_intron = INTRON, vep_consequence = Consequence, vep_impact = IMPACT, 
                                                        `5UTR_annotation`, `5UTR_consequence`, Existing_InFrame_oORFs, Existing_OutOfFrame_oORFs, Existing_uORFs,
                                                        processed_tx_flag, UTR_variant_flag, noncoding_transcript_variant_flag)]

smorf_join = smorf_vars[, .(variant_id, transcript_id,strand = STRAND,  smorf_id, smorf_consequence = Consequence, smorf_impact = IMPACT, 
                            cDNA_position, CDS_position, Protein_position, Amino_acids, Codons, smorf_exon = EXON, smorf_intron = INTRON, smorf_NMD = NMD,
                            gnomADg_AF, gnomADg_AFR_AF, gnomADg_AMR_AF, gnomADg_AMI_AF, gnomADg_ASJ_AF, gnomADg_EAS_AF, gnomADg_FIN_AF, gnomADg_MID_AF, gnomADg_NFE_AF, gnomADg_OTH_AF, gnomADg_SAS_AF,
                            SpliceAI_pred,SpliceAI_pred_type, DS_AG, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL, spliceai_max_delta,
                            Existing_variation, CLIN_SIG, ClinVar,CADD_PHRED, CADD_RAW, LOEUF, NMD, PhyloP)]

smorf_vep = vep_join[smorf_join, on = .(variant_id, transcript_id)]

smorf_vep = smorf_vep[, .(variant_id, transcript_id, strand, smorf_id, smorf_consequence, smorf_impact, vep_consequence, vep_impact, 
                          cDNA_position, CDS_position, Protein_position, Amino_acids, Codons, tx_exon, tx_intron, smorf_exon, smorf_intron, 
                          gnomADg_AF, gnomADg_AFR_AF, gnomADg_AMR_AF, gnomADg_AMI_AF, gnomADg_ASJ_AF, gnomADg_EAS_AF, gnomADg_FIN_AF, gnomADg_MID_AF, gnomADg_NFE_AF, gnomADg_OTH_AF, gnomADg_SAS_AF,
                          Existing_variation, CLIN_SIG, ClinVar,CADD_PHRED, CADD_RAW, LOEUF,smorf_NMD, NMD, PhyloP,
                          `5UTR_annotation`, `5UTR_consequence`, Existing_InFrame_oORFs, Existing_OutOfFrame_oORFs, Existing_uORFs,
                          SpliceAI_pred,SpliceAI_pred_type, DS_AG, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL, spliceai_max_delta,
                          processed_tx_flag, UTR_variant_flag, noncoding_transcript_variant_flag)]








### step 4- disentangling 'start_lost' variants
# annotate variants that went from noncanonical to noncanonical start, or noncanonical to canonical start
smorf_vep[grepl('start_lost', smorf_consequence)& grepl('/', Codons), ref_start_codon := toupper(sapply(strsplit(Codons, '/'), '[[', 1))]
smorf_vep[grepl('start_lost', smorf_consequence)& grepl('/', Codons), alt_start_codon := toupper(sapply(strsplit(Codons, '/'), '[[', 2))]




# indels which affect the start codon but retain the same codon are annotated as 'start_lost,start_retained'- remove start_lost consequence from these, it's confusing!
smorf_vep[grepl('start_retained_variant', smorf_consequence), smorf_consequence := gsub('start_lost,','', smorf_consequence)]
smorf_vep[grepl('start_retained_variant', smorf_consequence) & ref_start_codon != 'ATG', smorf_consequence := gsub('start_retained_variant', 'noncanonical_start_retained_variant', smorf_consequence)]

smorf_vep[grepl('start_lost', smorf_consequence) & ref_start_codon != 'ATG' & alt_start_codon == 'ATG', smorf_consequence := gsub('start_lost', 'canonical_start_gain_variant', smorf_consequence)]
smorf_vep[grepl('start_lost', smorf_consequence) & ref_start_codon != 'ATG' & alt_start_codon != 'ATG', smorf_consequence := gsub('start_lost', 'noncanonical_start_retained_variant', smorf_consequence)]





### step 5- reporting consequences / impacts of filtered variant sets
high_impact_consequences = c('transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 
                             'frameshift_variant', 'stop_lost', 'start_lost', 
                             'transcript_amplification', 'feature_elongation', 'feature_truncation', 'SpliceAI_0.8')
moderate_impact_consequences = c('inframe_insertion', 'inframe_deletion', 'missense_variant', 
                                 'protein_altering_variant', 'canonical_start_gain_variant', 'SpliceAI_0.2')
low_impact_consequences = c('splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant',
                            'splice_polypyrimidine_tract_variant', 'incomplete_terminal_codon_variant', 
                            'start_retained_variant', 'noncanonical_start_retained_variant', 'stop_retained_variant', 'synonymous_variant')

# integrating spliceAI into smorf variant impacts 
smorf_vep[spliceai_max_delta >= 0.8 & ! grepl('splice_donor_variant|splice_acceptor_variant', smorf_consequence), smorf_consequence := paste0(smorf_consequence, ',SpliceAI_0.8')]
smorf_vep[spliceai_max_delta >= 0.2 & ! grepl('SpliceAI_0.8|splice_donor_variant|splice_acceptor_variant', smorf_consequence), smorf_consequence := paste0(smorf_consequence, ',SpliceAI_0.2')]


# variants with downgraded impacts (start lost)
# downgraded to moderate
smorf_vep[!grepl(paste0(high_impact_consequences, collapse = '|'), smorf_consequence) & smorf_impact == 'HIGH' &
            grepl(paste0(moderate_impact_consequences, collapse = '|'), smorf_consequence), smorf_impact := 'MODERATE']
# downgraded to low
smorf_vep[!grepl(paste0(high_impact_consequences, collapse = '|'), smorf_consequence) & smorf_impact == 'HIGH' &
            grepl(paste0(low_impact_consequences, collapse = '|'), smorf_consequence), smorf_impact := 'LOW']

# variants with upgraded impacts (splicing)
smorf_vep[grepl(paste0(high_impact_consequences, collapse = '|'), smorf_consequence) & smorf_impact != 'HIGH', smorf_impact := 'HIGH']
smorf_vep[grepl(paste0(moderate_impact_consequences, collapse = '|'), smorf_consequence) &! smorf_impact %in% c('MODERATE', 'HIGH'), smorf_impact := 'MODERATE']


# annotate variants as 'processed transcript high impact': VEP
smorf_vep[processed_tx_flag == 1, smorf_impact := 'processed_transcript']



# get number of variants mapped to multiple smorfs
smorfs_per_var1 = unique(smorf_vep[, .(variant_id,smorf_id)]) %>% group_by(variant_id) %>% tally() %>% ungroup() %>%
  mutate(smorfs_per_var = n) %>% mutate(smorfs_per_var = ifelse(smorfs_per_var > 2, 'over2', smorfs_per_var)) %>% 
  group_by(smorfs_per_var) %>% tally() %>% 
  mutate(smorfs_per_var = paste0('smorfs_per_var_', smorfs_per_var)) %>%
  data.table::transpose() 

smorfs_per_var = unique(smorf_vep[, .(variant_id,smorf_id)]) %>% group_by(variant_id) %>% tally() %>% ungroup() %>%
  mutate(smorfs_per_var = n) %>% mutate(smorfs_per_var = ifelse(smorfs_per_var > 2, 'over2', smorfs_per_var)) %>% 
  group_by(smorfs_per_var) %>% tally() %>% 
  mutate(smorfs_per_var = paste0('smorfs_per_var_', smorfs_per_var)) %>%
  data.table::transpose() %>% row_to_names(row_number = 1)

this_file_stats = cbind(this_file_stats, smorfs_per_var)
this_file_stats[, total_rows := nrow(smorf_vep)]

# get smorf consequence counts
smorf_consequences = data.frame(table(smorf_vep$smorf_consequence))%>%
  mutate(Var1 = paste0('smorf_conseq_', Var1))  %>%
  data.table::transpose() %>% row_to_names(row_number = 1)

this_file_stats = cbind(this_file_stats, smorf_consequences)


smorf_impact = data.frame(table(tolower(smorf_vep$smorf_impact))) %>%
  mutate(Var1 = paste0('smorf_impact_', Var1)) %>%
  data.table::transpose() %>% row_to_names(row_number = 1)

this_file_stats = cbind(this_file_stats, smorf_impact)


fwrite(smorf_vep[smorf_impact != 'MODIFIER'], 'vep_results_filt.tsv.gz', sep = '\t')

fwrite(this_file_stats, 'filtering_summary.tsv', sep = '\t')



