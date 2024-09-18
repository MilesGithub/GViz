library(Gviz)
library(GenomicRanges)

all_gene_granges <- readRDS('data/all_gene_granges.rds')
all_gistic_granges <- readRDS('data/all_gistic_granges.rds')
all_CNAs <- readRDS('data/all_TCGA_CNAs.rds')
clinical <- readRDS('data/TCGA_clinical.rds')

all_CNAs <- all_CNAs[all_CNAs$width < 100000, ]
all_gistic_granges <- all_gistic_granges[all_gistic_granges$peak_width >= 1000, ]
indexes <- all_gistic_granges$index

for (i in 1:length(indexes)) {
  print(i)

  index <- indexes[i]
  gistic_grange <- all_gistic_granges[all_gistic_granges$index == index]
  start <- max(as.numeric(start(gistic_grange)) - 500000, 0)
  end <- as.numeric(end(gistic_grange)) + 500000
  chr <- as.character(seqnames(gistic_grange))
  cancer_type <- as.character(gistic_grange$cancer_type)
  CNA_type <- as.character(gistic_grange$CNA_type)
  
  grange <- GRanges(Rle(chr), ranges = IRanges(start = start, end = end))
  
  # Tracks
  gtrack <- GenomeAxisTrack()
  gistic_track <- AnnotationTrack(gistic_grange, genome = "hg19", name = "gistic", rotation.title = 0, background.title ='#808080', fill = "#ffcc00")
  biomartTrack <- BiomartGeneRegionTrack(
    rotation.title = 0, background.title ='#808080', genome = "hg19", chromosome = chr,
    start = start, end = end, stacking = "squish", collapseTranscripts = "meta",
    name = "ENSEMBL", fill = "#ebb52e", filters = list(biotype = "protein_coding")
  )
  itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

  all_CNA_grange <- GRanges(Rle(all_CNAs$Chromosome), ranges = IRanges(start = all_CNAs$Start, end = all_CNAs$End))
  all_CNA_grange$patient_id <- all_CNAs$patient_id
  all_CNA_grange$Segment_Mean <- all_CNAs$Segment_Mean
  
  clinical_sub <- clinical[clinical$acronym %in% cancer_type, ]
  all_CNA_grange_sub <- all_CNA_grange[all_CNA_grange$patient_id %in% clinical_sub$bcr_patient_barcode, ]
  
  overlaps <- as.data.frame(findOverlaps(gistic_grange, all_CNA_grange_sub))
  
  overlap_CNAs <- all_CNA_grange_sub[overlaps$subjectHits]
  overlap_CNAs_df <- data.frame(overlap_CNAs)
  overlap_CNAs_df <- overlap_CNAs_df[order(-overlap_CNAs_df$width), ]
  
  patient_ids <- unique(overlap_CNAs_df$patient_id)
  
  if (length(patient_ids) >= 10) {
    tracks <- c(itrack, gtrack, biomartTrack, gistic_track)
    
    for (patient_id in patient_ids) {
      overlap_CNAs_sub <- overlap_CNAs[overlap_CNAs$patient_id == patient_id, ]
      overlap_CNAs_sub_reduced <- reduce(overlap_CNAs_sub)
      segment_mean <- mean(overlap_CNAs$Segment_Mean)

      if(segment_mean > 0){
        fill_color <- '#aa0000'
      }else if(segment_mean < 0){
        fill_color <- '#003380'
      }
      
      track <- AnnotationTrack(
        overlap_CNAs_sub_reduced, name = patient_id, rotation.title = 0, 
        genome = "hg19", fill = fill_color, fontcolor.group = "black", background.title ='#4d4d4d', showTitle = TRUE
      )
      
      tracks <- c(tracks, track)
    }

    filename <- paste0('results/', index, '_', cancer_type, '_', CNA_type, '_', chr, '_', start, '_', end, '.pdf')
    height<-5
    if(nrow(overlap_CNAs_df)>10){
      height<-ceiling(nrow(overlap_CNAs_df)/2)
    }
    
    pdf(filename, width = 10, height = height)
    plotTracks(tracks, from = start, to = end, main=paste0(cancer_type, ' ', CNA_type, ' ', chr, ':', start,'-', end),  title.width = 4, cex.title = 0.3, transcriptAnnotation = "symbol", cex.main = 1)
    dev.off()
  }
}
