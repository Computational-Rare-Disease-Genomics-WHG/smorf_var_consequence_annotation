
library(data.table)
library(tidyverse)

# get exon number where CDS starts for each transcript
gff = rtracklayer::import('data/Homo_sapiens.GRCh38.108.chr.gtf.gz')
gff_dt = as.data.table(gff)


CDS_exons = gff_dt[type == 'CDS']
CDS_exons[, `:=` (CDS_start = min(start), CDS_end = max(end)), by = .(transcript_id)]
CDS_bounds = CDS_exons[CDS_start == start | CDS_end == end]
setorder(CDS_bounds, transcript_id, exon_number)

CDS_bounds = CDS_bounds[, .(transcript_id, exon_number = as.numeric(exon_number))]

CDS_bounds = CDS_bounds[, .(exon_number, cds = ifelse(exon_number == min(exon_number), 'cds_start', 'cds_end')), 
                        by = .(transcript_id)]
CDS_bounds = CDS_bounds %>% pivot_wider(id_cols = transcript_id, names_from= cds, values_from = exon_number) %>%
  mutate(cds_end = ifelse(is.na(cds_end), cds_start, cds_end) ) 

fwrite(CDS_bounds,'vep_results_filt/resources/home/dnanexus/CDS_exons.csv')



smorfs = fread('data/smorfs_public_highconf.bed')
fwrite(smorfs, 'vep_results_filt/resources/home/dnanexus/smorfs_public_highconf.bed')
smorfs = smorfs[name == cluster]

smorf_cds = smorfs[, .(Feature = smorf_id, transcript_id = tx_id, start = chromStart, end = chromEnd)]


ensembl_exons = gff_dt[type == 'exon', .(transcript_id, start, end, exon_number)]

setkey(ensembl_exons, transcript_id, start, end)

smorf_cds_start = smorf_cds[, .(Feature, transcript_id, start, end = start, cds = 'smorf_start_exon')]
smorf_cds_start = foverlaps(smorf_cds_start, ensembl_exons[transcript_id %in% smorf_cds$transcript_id, .(transcript_id, start, end, exon_number)], by.x = c('transcript_id', 'start', 'end'), type = 'within')


smorf_cds_end = smorf_cds[, .(Feature, transcript_id, start = end, end, cds = 'smorf_end_exon')]
smorf_cds_end = foverlaps(smorf_cds_end, ensembl_exons[transcript_id %in% smorf_cds$transcript_id, .(transcript_id, start, end, exon_number)], by.x = c('transcript_id', 'start', 'end'), type = 'within')

smorf_cds_bounds = rbind(smorf_cds_start, smorf_cds_end)
setorder(smorf_cds_bounds, Feature, exon_number)

smorf_cds_bounds = smorf_cds_bounds[, .(Feature, transcript_id, exon_number, cds)]

smorf_cds_bounds = smorf_cds_bounds %>% pivot_wider(id_cols = c('Feature', 'transcript_id'), names_from= cds, values_from = exon_number) 
setDT(smorf_cds_bounds)
smorf_cds_bounds[smorf_start_exon > smorf_end_exon, `:=` (smorf_start_exon = smorf_end_exon, smorf_end_exon = smorf_start_exon )]

fwrite(smorf_cds_bounds, 'vep_results_filt/resources/home/dnanexus/smorf_exon_bounds.csv')



#making transcript_CDS_108.csv file
#transcript	chr	strand	tran_start	tran_end	cds_start	cds_end

ens_tx_dt = gff_dt[type == 'transcript']
transcript_CDS_108 = ens_tx_dt[, .(transcript = transcript_id, chr = seqnames, strand, tran_start = start, tran_end = end)]

cds_coords = CDS_exons[, .(cds_start = min(start), cds_end = max(end)), by = .(chr = seqnames, strand, transcript = transcript_id)]

transcript_CDS_108 = cds_coords[transcript_CDS_108, on = .(transcript, chr, strand)]

fwrite(transcript_CDS_108[, .(transcript, chr, strand, tran_start, tran_end, cds_start, cds_end)], 'vep_results_filt/resources/home/dnanexus/transcript_CDS_108.csv')

