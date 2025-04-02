# --- Load Libraries ---
# List of packages to check and install
packages <- c("BiocManager", "Gviz", "GenomicRanges", "BiocParallel")

# Function to check if a package is installed
load_packages <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% rownames(installed.packages())) {
      print(paste0("Package already installed: ", pkg))
    } else if (pkg %in% BiocManager::available()) {
      BiocManager::install(pkg)
      print(paste0("Installed from Bioconductor: ", pkg))
    } else {
      install.packages(pkg)
      print(paste0("Installed from CRAN: ", pkg))
    }
  }
  library(pkg, character.only = TRUE)
  print(paste0("Loaded: ", pkg))
}

lapply(packages, load_packages)

# --- Constants ---
# Input files
all_gene_granges_file <- 'data/all_gene_granges.rds'
all_gistic_granges_file <- 'data/all_gistic_granges.rds'
all_cnas_file <- 'data/all_TCGA_CNAs.rds'
clinical_file <- 'data/TCGA_clinical.rds'

# Output
output_dir <- 'results_improved'
min_patients_for_plot <- 10
plot_width <- 10
default_plot_height <- 5
height_scaling_factor <- 0.5

# Parameters
window_padding <- 500000
min_gistic_peak_width <- 1000
max_cna_width_filter <- 100000

# Plotting Aesthetics
gistic_fill <- "#ffcc00"
ensembl_fill <- "#ebb52e"
amp_fill <- '#aa0000'
del_fill <- '#003380'
track_bg_title <- '#808080'
patient_track_bg_title <- '#4d4d4d'
genome_version <- "hg19"

# --- Load Data ---
message("Loading data...")

all_gistic_granges <- readRDS(all_gistic_granges_file)
all_CNAs <- readRDS(all_cnas_file)
clinical <- readRDS(clinical_file)
message("Data loaded.")

# --- Pre-processing ---
message("Preprocessing data...")

# Filter input CNAs

initial_cna_count <- nrow(all_CNAs)
all_CNAs <- all_CNAs[all_CNAs$width < max_cna_width_filter, ]
message("Filtered CNAs by width (< ", max_cna_width_filter, "bp): ",
        initial_cna_count, " -> ", nrow(all_CNAs))

# Filter GISTIC peaks
initial_gistic_count <- length(all_gistic_granges)
all_gistic_granges <- all_gistic_granges[all_gistic_granges$peak_width >= min_gistic_peak_width, ]
message("Filtered GISTIC peaks by width (>= ", min_gistic_peak_width, "bp): ",
        initial_gistic_count, " -> ", length(all_gistic_granges))

message("Creating GRanges for all CNAs...")
all_CNA_grange <- GRanges(
  seqnames = Rle(all_CNAs$Chromosome),
  ranges = IRanges(start = all_CNAs$Start, end = all_CNAs$End),
  patient_id = all_CNAs$patient_id,
  Segment_Mean = all_CNAs$Segment_Mean
)

# Standardize chromosome naming convention
seqlevelsStyle(all_CNA_grange) <- "UCSC"
message("Finished creating CNA GRanges.")

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  message("Creating output directory: ", output_dir)
  dir.create(output_dir, recursive = TRUE)
}

# --- Main Loop ---
message("Starting GISTIC peak processing loop...")
num_peaks <- length(all_gistic_granges)

for (i in seq_along(all_gistic_granges)) {
  
  current_gistic_grange <- all_gistic_granges[i]
  index <- current_gistic_grange$index
  
  message("Processing peak ", i, "/", num_peaks, " (Index: ", index, ")")
  
  chr <- as.character(seqnames(current_gistic_grange))

  if (!startsWith(chr, "chr")) {
    chr <- paste0("chr", chr)
  }
  peak_start <- start(current_gistic_grange)
  peak_end <- end(current_gistic_grange)
  cancer_type <- as.character(current_gistic_grange$cancer_type)
  cna_type <- as.character(current_gistic_grange$CNA_type)
  
  # Define plotting window
  plot_start <- max(peak_start - window_padding, 1)
  plot_end <- peak_end + window_padding

  # Create GRanges for the plotting window
  plotting_window_grange <- GRanges(
    seqnames = Rle(chr),
    ranges = IRanges(start = plot_start, end = plot_end)
  )
  seqlevelsStyle(plotting_window_grange) <- "UCSC"
  
  # --- Prepare Tracks ---
  
  # Static tracks for this iteration
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = genome_version, chromosome = chr)
  
  # GISTIC peak track
  seqlevelsStyle(current_gistic_grange) <- "UCSC"
  gistic_track <- AnnotationTrack(
    current_gistic_grange,
    genome = genome_version,
    chromosome = chr,
    name = "GISTIC Peak",
    rotation.title = 0,
    background.title = track_bg_title,
    fill = gistic_fill,
    feature = paste(cna_type, "Peak"),
    showFeatureId = TRUE
  )
  
  # Gene track (BioMart query)
  biomartTrack <- tryCatch({
    BiomartGeneRegionTrack(
      genome = genome_version,
      chromosome = chr,
      start = plot_start,
      end = plot_end,
      stacking = "squish",
      collapseTranscripts = "meta",
      name = "ENSEMBL Genes",
      fill = ensembl_fill,
      rotation.title = 0,
      background.title = track_bg_title,
      filters = list(biotype = "protein_coding"),
      transcriptAnnotation = "symbol",
      cex.title = 0.8,
      cex.group = 0.7
    )
  }, error = function(e) {
    message("Warning: Biomart query failed for ", chr, ":", plot_start, "-", plot_end, ". Error: ", e$message)
    AnnotationTrack(name="BioMart Error", genome=genome_version, chromosome=chr)
  })
  
  
  # --- Find Overlapping CNAs for relevant patients ---
  
  # Filter clinical data for the specific cancer type
  clinical_sub <- clinical[clinical$acronym %in% cancer_type, ]
  
  if (nrow(clinical_sub) == 0) {
    message("  No clinical data found for cancer type: ", cancer_type, ". Skipping.")
    next
  }
  
  # Filter the pre-computed CNA GRanges for patients in this cancer type
  relevant_patient_ids <- clinical_sub$bcr_patient_barcode
  all_CNA_grange_sub <- all_CNA_grange[all_CNA_grange$patient_id %in% relevant_patient_ids]
  
  if (length(all_CNA_grange_sub) == 0) {
    message("  No CNA data found for patients in cancer type: ", cancer_type, ". Skipping.")
    next
  }
  
  # Find overlaps between the GISTIC peak and the subsetted CNAs
  overlaps <- findOverlaps(current_gistic_grange, all_CNA_grange_sub)
  
  if (length(overlaps) == 0) {
    message("  No CNAs overlap the GISTIC peak region. Skipping.")
    next
  }
  
  # Get the actual CNA segments that overlap the GISTIC peak
  overlap_CNAs <- all_CNA_grange_sub[subjectHits(overlaps)]
  
  # Get unique patient IDs with overlapping CNAs
  patient_ids_with_overlaps <- unique(overlap_CNAs$patient_id)
  num_overlapping_patients <- length(patient_ids_with_overlaps)
  
  message("  Found ", num_overlapping_patients, " patients with CNAs overlapping the GISTIC peak.")
  
  # --- Generate Plot if Enough Patients Overlap ---
  if (num_overlapping_patients >= min_patients_for_plot) {
    message("  Sufficient patients (>= ", min_patients_for_plot, "). Generating plot...")
    
    # Initialize track list
    tracks_to_plot <- c(itrack, gtrack, biomartTrack, gistic_track)
    
    # Calculate the *overall* mean segment value for CNAs overlapping this GISTIC peak
    overall_segment_mean <- mean(overlap_CNAs$Segment_Mean, na.rm = TRUE)
    message("    Overall mean Segment_Mean for overlapping CNAs: ", round(overall_segment_mean, 3))
    
    # Determine fill color based on the overall mean (or based on GISTIC cna_type if preferred)
    if (!is.na(overall_segment_mean) && overall_segment_mean > 0) {
      patient_track_fill <- amp_fill
      message("    Using Amplification color (", patient_track_fill, ") based on positive overall mean.")
    } else if (!is.na(overall_segment_mean) && overall_segment_mean < 0) {
      patient_track_fill <- del_fill
      message("    Using Deletion color (", patient_track_fill, ") based on negative overall mean.")
    } else {
      patient_track_fill <- "grey" # Neutral color if mean is zero or NA
      message("    Using Neutral color (", patient_track_fill, ") based on zero/NA overall mean.")
    }
    
    # Create tracks for each patient
    patient_ids_with_overlaps <- sort(patient_ids_with_overlaps)
    plotted_patient_count <- 0
    
    for (patient_id in patient_ids_with_overlaps) {
      
      # Get CNAs for the current patient that overlap the GISTIC peak
      patient_CNAs <- overlap_CNAs[overlap_CNAs$patient_id == patient_id]

      patient_CNAs_in_window <- subsetByOverlaps(patient_CNAs, plotting_window_grange)
      
      if (length(patient_CNAs_in_window) > 0) {

        patient_CNAs_reduced <- reduce(patient_CNAs_in_window)
        
        # Create the track for this patient
        patient_track <- AnnotationTrack(
          patient_CNAs_reduced,
          name = patient_id,
          rotation.title = 0,
          genome = genome_version,
          chromosome = chr,
          fill = patient_track_fill,
          background.title = patient_track_bg_title,
          showTitle = TRUE,
          cex.title = 0.6,
          showFeatureId = FALSE
        )
        tracks_to_plot <- c(tracks_to_plot, patient_track)
        plotted_patient_count <- plotted_patient_count + 1
      } else {
        message("    Skipping patient track for ", patient_id, " - no segments within the actual plot window.")
      }
    }
    
    # --- Save Plot ---

    plot_height <- default_plot_height + max(0, plotted_patient_count - min_patients_for_plot) * height_scaling_factor
    
    original_chr <- as.character(seqnames(all_gistic_granges[i]))
    plot_filename <- sprintf(
      "%s_%s_%s_%s_%d_%d.pdf",
      index, cancer_type, cna_type, original_chr, plot_start, plot_end
    )
    plot_filepath <- file.path(output_dir, plot_filename)
    
    # Generate plot
    message("    Saving plot to: ", plot_filepath)
    pdf(plot_filepath, width = plot_width, height = plot_height)
    tryCatch({
      plotTracks(
        tracks_to_plot,
        from = plot_start,
        to = plot_end,
        chromosome = chr,
        main = sprintf("%s %s (%s:%d-%d)", cancer_type, cna_type, chr, plot_start, plot_end),
        title.width = 2.5,
        cex.main = 0.9,
        cex.title = 0.7
      )
    }, error = function(e) {
      message("ERROR: Failed to generate plot for index ", index, ". Error: ", e$message)
      plot.new()
      title(main=paste("Error plotting:", index), sub=e$message)
    }, finally = {
      if(length(dev.list()) > 0) {
        dev.off()
      }
    })

    if(exists("plot_filepath") && file.exists(plot_filepath) && length(dev.list()) > 0 && names(dev.cur()) == "pdf") {
      dev.off()
    }
    
    
  } else {
    message("  Skipping plot generation: Not enough patients with overlaps (", num_overlapping_patients, " < ", min_patients_for_plot, ").")
  }
  
} # End of main loop

message("Processing complete.")