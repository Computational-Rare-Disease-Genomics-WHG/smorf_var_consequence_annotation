
library(data.table)
library(GenomicRanges)
library(dplyr)
ens_gtf <- rtracklayer::import('data/Homo_sapiens.GRCh38.108.chr.gtf.gz')
ens_gtf_dt = as.data.table(ens_gtf)
smorfs_bed = fread('data/smorfs_public_highconf.bed')


smorfs_gff = ens_gtf_dt[type %in% c('gene', 'transcript', 'exon') & gene_id %in% smorfs_bed$gene_id,
                        .(seqnames, source = 'smORF-riboseq', feature = type, start, end, score = '.', strand, frame = '.', gene_id, transcript_id, exon_number, exon_id)]

genes_gff = smorfs_gff[feature == 'gene']

# add smORF-IDs to gtf to replace transcript ids - will duplicate rows for transcripts with multiple smORFs mapped
smorf_tx = unique(smorfs_bed[tx_id != 'unmapped', .(gene_id, transcript_id = tx_id, smorf_id = smorf_id)])
smorfs_gff = smorfs_gff[smorf_tx, on = .(gene_id, transcript_id), allow.cartesian=TRUE]
smorfs_gff = rbind(smorfs_gff, genes_gff, fill = TRUE)

# intersect gtf exons w/ smORF start & end coordinates to get CDS rows for smORF
smorf_gr = GRanges(smorfs_bed[, .(seqnames = smorf_id, start = chromStart, end = chromEnd, strand)])
tx_exons_gr = GRanges(smorfs_gff[feature == 'exon', .(seqnames = smorf_id, start, end, strand, exon_number)])

smorf_cds_gr <- GenomicRanges::intersect(smorf_gr, tx_exons_gr)

# Add the metadata from gr2
m <- findOverlaps(smorf_gr, tx_exons_gr)
mcols(smorf_cds_gr) <- cbind.data.frame(mcols(tx_exons_gr[subjectHits(m)]))

smorf_cds = as.data.table(smorf_cds_gr)
gtf_info = unique(smorfs_gff[feature == 'transcript', .(gene_id, transcript_id, smorf_id, seqnames, source, score, frame)])
smorf_cds = smorf_cds[, .(smorf_id = seqnames, start, end, strand, feature = 'CDS', exon_number)][gtf_info, on = .(smorf_id)]

smorfs_gff = rbind(smorfs_gff, smorf_cds, fill = TRUE)




#gene
smorfs_genes = smorfs_gff[feature == 'gene', .(seqnames, source, feature, start, end, score, strand, frame, 
                                               attribute = paste0('gene_id "', gene_id, '";'))]

#transcript
smorfs_tx = smorfs_gff[feature == 'transcript', .(seqnames, source, feature, start, end, score, strand, frame, 
                                                  attribute = paste0('gene_id "', gene_id, '"; ',
                                                                     'transcript_id "', smorf_id, '"; ',
                                                                     'transcript_biotype "', "protein_coding", '";'))]
#exon
smorfs_exons = smorfs_gff[feature == 'exon', .(seqnames, source, feature, start, end, score, strand, frame, 
                                               attribute = paste0('gene_id "', gene_id, '"; ',
                                                                  'transcript_id "', smorf_id, '"; ',
                                                                  'exon_number "', exon_number, '"; ',
                                                                  'exon_id "', exon_id, '"; '))]

#cds
smorfs_cds = smorfs_gff[feature == 'CDS', .(seqnames, source, feature, start, end, score, strand, frame, 
                                            attribute = paste0('gene_id "', gene_id, '"; ',
                                                               'transcript_id "', smorf_id, '"; ',
                                                               'exon_number "', exon_number, '"; ',
                                                               'ccds_id "', smorf_id, '"; '))]

smorfs_gff = rbind(smorfs_genes, smorfs_tx, smorfs_exons, smorfs_cds)






# adding unmapped smORFs

unmapped_smorfs_gff = rbind(smorfs_bed[grepl('unmapped', smorf_id), .(seqnames = gsub('chr', '', chrom), source = 'smORF-riboseq', feature = 'gene', start = chromStart-5, end = chromEnd+5, score = '.', strand, frame = '.', 
                                                                     attribute = paste0('gene_id "', smorf_id, '_gene";'))],
                            smorfs_bed[grepl('unmapped', smorf_id), .(seqnames = gsub('chr', '', chrom), source = 'smORF-riboseq', feature = 'transcript', start = chromStart, end = chromEnd, score = '.', strand, frame = '.', 
                                                                     attribute = paste0('gene_id "', smorf_id, '_gene"; ',
                                                                                        'transcript_id "', smorf_id, '"; ',
                                                                                        'transcript_biotype "', "protein_coding", '";'))])

# making smORF exons
exons_dt = smorfs_bed[grepl('unmapped', smorf_id), .(blockSizes = as.numeric(unlist(strsplit(blockSizes, '\\,'))),
                                                    blockStarts = as.numeric(unlist(strsplit(blockStarts, '\\,')))), 
                      by = setdiff(names(smorfs_bed), c('blockSizes', 'blockStarts'))]
exons_dt[, `:=` (exon_start = chromStart + blockStarts,
                 exon_end = chromStart + blockStarts + blockSizes)]
exons_dt = rbind(arrange(exons_dt[strand == '+'], smorf_id, exon_start),
                 arrange(exons_dt[strand == '-'], smorf_id, -exon_start))
exons_dt[, exon_number := rowid(smorf_id)] 

exons = exons_dt[, .(seqnames = gsub('chr', '', chrom), source = 'smORF-riboseq', feature = 'exon', 
                     start = exon_start, end = exon_end, score = '.', strand, frame = '.',
                     attribute = paste0('gene_id "', smorf_id, '_gene"; ',
                                        'transcript_id "', smorf_id, '"; ',
                                        'exon_number "', exon_number, '"; ',
                                        'exon_id "', paste0(smorf_id, '_', exon_number), '"; '))]

cds = exons_dt[, .(seqnames = gsub('chr', '', chrom), source = 'smORF-riboseq', feature = 'CDS', 
                   start = exon_start, end = exon_end, score = '.', strand, frame = '.',
                   attribute = paste0('gene_id "', smorf_id, '_gene"; ',
                                      'transcript_id "', smorf_id, '"; ',
                                      'exon_number "', exon_number, '"; ',
                                      'ccds_id "', smorf_id, '"; '))]

unmapped_smorfs_gff = rbind(unmapped_smorfs_gff, exons, cds)
setorder(unmapped_smorfs_gff, seqnames, start)





smorfs_gff = rbind(smorfs_gff, unmapped_smorfs_gff)
setorder(smorfs_gff, seqnames, start)
fwrite(smorfs_gff, 'data/smorfs_public_highconf.gtf', sep = '\t', col.name = FALSE, quote=F)

system('grep -v "#" data/smorfs_public_highconf.gtf | sort -k1,1 -k4,4n -k5,5n -t$"\t" | bgzip -c > data/smorfs_public_highconf.gtf.gz')
system('tabix data/smorfs_public_highconf.gtf.gz')
system('rm data/smorfs_public_highconf.gtf')
















