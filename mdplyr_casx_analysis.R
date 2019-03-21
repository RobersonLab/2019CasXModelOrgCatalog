############
# packages #
############
library( tidyverse )
library( knitr )
library( here )
library( GenomicRanges )
library( GenomicFeatures )
library( multidplyr )

###################
# setup multiplyr #
###################
r_cpus <- read_csv( here::here( 'changeable_parameters', 'num_r_cpus.txt' ), col_names = 'cpu' ) %>%
	pull( 'cpu' ) %>%
	as.integer( . )

cluster <- create_cluster( r_cpus )
set_default_cluster( cluster )

#############
# functions #
#############
read_fasta_index <- function( filename )
{
	read_tsv( file = filename, col_names = c( 'name', 'size', 'offset', 'basesPerLine', 'bytesPerLine' ), col_types = c( col_character(), col_integer(), col_integer(), col_integer(), col_integer() ), progress = FALSE, skip_empty_rows = TRUE ) %>%
		return( . )
}

read_motif_csv <- function( filename ) {
  read_csv( file = filename, col_names = TRUE, col_types = 'ciiccc', progress = FALSE, skip_empty_rows = TRUE ) %>%
	dplyr::select( .data = ., -Motif ) %>%
	mutate( .data = ., PAM = str_extract( string = Sequence, pattern = "^[ACGT]{4}" ) ) %>%
	mutate( .data = ., Target = str_extract( string = Sequence, pattern = "[ACGT]{20}$" ) ) %>%
	return( . )
}

motif_to_GRanges <- function( motif_obj, fai_obj, genome_name )
{
	# deal with mito being circular
	fai_obj <- fai_obj %>%
	mutate( .data = ., circ = case_when(
		name %in% c( "MT", "Mt", "M", 'MtDNA', 'Mito', 'chrMT', "chrMt", "chrM", "chrMtDNA", "chrMito" ) ~ TRUE,
		TRUE ~ FALSE
	))

	# seq info from fai file
	mySeqInfo <- Seqinfo( seqnames = fai_obj$name,
		seqlengths = fai_obj$size,
		isCircular = fai_obj$circ,
		genome = genome_name )

	# GRanges object
	GRanges( seqnames = motif_obj$Contig,
		ranges = IRanges( start = motif_obj$Start, end = motif_obj$End ),
		strand = motif_obj$Strand,
		mcols = motif_obj[ , c( "uid", "Target" ) ],
		seqinfo = mySeqInfo ) %>%
	  return( . )
}

file_exists_not_empty <- function( filename ) {
	return( file.exists( filename ) && file.info( filename )$size > 0 )
}

########################
# data capture objects #
########################
species <- read_tsv( file = here::here( 'changeable_parameters', 'species.txt' ), col_names = 'species' ) %>%
	pull( species )

#################################################
# analyze each species, accumulating group data #
#################################################
for ( species_index in 1:length( species ) )
{
	#####################
	# grab species name #
	#####################
	curr_species <- species[ species_index ]

	message( paste0( "Now processing ", curr_species ) )

	######################
	# generate filenames #
	######################
	cuts_per_gene_filepath <- here::here( 'R_analysis', paste0( curr_species, "_cutsPerGene.csv" ) )
	annotated_cuts_filepath <- here::here( 'R_analysis', paste0( curr_species, "_annotated_Casx_sites.csv" ) )
	sites_filepath <- here::here( 'R_analysis', paste0( curr_species, "_site_count.csv" ) )
	pam_sites_filepath <- here::here( 'R_analysis', paste0( curr_species, "_pam_count.csv" ) )
	genome_size_filepath <- here::here( 'R_analysis', paste0( curr_species, "_genome_size.csv" ) )
	target_list_filepath <- here::here( 'R_analysis', paste0( curr_species, "_targets.csv" ) )

	###############
	# file status #
	###############

	# if these all exist, skip to next
	if ( all(
	file_exists_not_empty( cuts_per_gene_filepath ),
	file_exists_not_empty( annotated_cuts_filepath ),
	file_exists_not_empty( sites_filepath ),
	file_exists_not_empty( pam_sites_filepath ),
	file_exists_not_empty( genome_size_filepath ),
	file_exists_not_empty( target_list_filepath ) ) ) {
		message( paste0( "Skipping all for ", curr_species ) )
		next
	}

	#########################################
	# grab FASTA index info for genome size #
	#########################################
	fai <- read_fasta_index( file = here::here( 'fasta', paste0( curr_species, ".fa.fai" ) ) )

	# Mb genome
	data.frame( MbGenome = sum( as.numeric( fai$size ) ) / 1E6 ) %>%
		write_csv( x = ., path = genome_size_filepath )

	#######################
	# load CasX site list #
	#######################
	# for motif sites, minus strand start and end are actual start and end
	# that makes GRanges freak out because start > end
	# need to flip that here.
	message( paste0( "Loading CasX sites for ", curr_species ) )
	
	all_casx_sites <- read_motif_csv( filename = here::here( 'casx_sites_out', paste0( curr_species, '_casx_sites.csv.gz' ) ) ) %>%
		mutate( .data = ., NewStart = case_when(
			Strand == "+" ~ Start,
			TRUE ~ End
		)) %>%
		mutate( .data = ., End = case_when(
			Strand == "+" ~ End,
			TRUE ~ Start
		)) %>%
		dplyr::select( -Start ) %>%
		rename( x = ., c( 'NewStart' = 'Start' ) ) %>%
		mutate( .data = ., uid = paste( Contig, as.character( Start ), as.character( End ), Strand, sep = "_" ) )

	# get count of site occurences
	site_match_count <- all_casx_sites %>%
		group_by( Target ) %>%
		summarize( ExactMatches = n() ) %>%
		as.data.frame( . )
		
	# write unique sites
	site_match_count %>%
		dplyr::select( Target ) %>%
		write_csv( x = ., path = target_list_filepath )

	# pull list of unique targets
	message( paste0( "Getting unique sites for ", curr_species ) )
	
	unique_target_seqs <- filter( .data = site_match_count, ExactMatches == 1 ) %>%
		pull( .data = ., Target )

	# write count of total and unique sites
	data.frame( total_sites = dim( all_casx_sites )[1], unique_sites = length( unique_target_seqs ) ) %>%
		write_csv( x = ., path = sites_filepath )

	# add to all_casx_sites
	all_casx_sites <- merge( x = all_casx_sites, y = site_match_count, by = 'Target', all = TRUE ) %>%
		mutate( Unique = case_when(
			ExactMatches == 1 ~ TRUE,
			TRUE ~ FALSE )
		)
		
	rm( site_match_count )

	#############
	# PAM usage #
	#############
	all_casx_sites %>%
		group_by( Unique, PAM ) %>%
		summarise(
			AllPam = n()
		) %>%
		as.data.frame( . ) %>%
		write_csv( x = ., path = pam_sites_filepath )

	# this is a long analysis
	# if you had to rerun the pipeline more than once troubleshooting
	# because you maybe ran out of memory or the computer rebooted
	# you don't want to regenerate all files
	# this logic will short circuit and skip files that have already been made
	if ( all(
	file_exists_not_empty( cuts_per_gene_filepath ),
	file_exists_not_empty( annotated_cuts_filepath ) ) ) {
		message( paste0( "Skipping overlaps for ", curr_species ) )
		next
	}

	################
	# import files #
	################
	message( paste0( "Loading gene database from GTF file for ", curr_species ) )
	
	granges_exons <- makeTxDbFromGFF( file = here::here( 'gtfs', paste0( curr_species, ".gtf" ) ), format="gtf" ) %>%
		exonsBy( ., by = 'gene' )
	
	names_of_exons <- names( granges_exons )
	unique_gene_names <- unique( names_of_exons )

	#####################
	# convert to ranges #
	#####################
	message( paste0( "Loading CasX into GRanges for ", curr_species ) )
	granges_casx <- motif_to_GRanges( all_casx_sites, fai, curr_species )

	####################
	# overlaps - exons #
	####################
	message( paste0( "Finding overlaps for ", curr_species ) )
	olap <- findOverlaps( granges_casx, granges_exons, ignore.strand = TRUE )

	olap_exon_table <- mcols( granges_casx )[ queryHits( olap ) , ] %>%
		as.data.frame( . ) %>%
		rename( x = ., c( 'mcols.uid' = 'uid', 'mcols.Target' = 'Target' ) ) %>%
		cbind( ., gene_id = names_of_exons[ subjectHits( olap ) ] )

	# little housekeeping
	# some of these are relatively large memory operations
	# we'll delete some unused objects to free up RAM
	rm( olap )
	rm( granges_casx )
	rm( granges_exons )

	#######################
	# find gRNAs per gene #
	#######################
	# I use the gene_ids and gRNA_Seq for this
	# We're interested in total and unique cutters.
	# That is derived from the grna_seqs and not the PAM site
	message( paste0( "Per gene stats for ", curr_species ) )
	
	total_grna_per_gene <- olap_exon_table %>%
		group_by( gene_id ) %>%
		summarize(
			TotalSites = n()
		)

	unique_grna_per_gene <- olap_exon_table %>%
		filter( .data = ., Target %in% unique_target_seqs ) %>%
		group_by( gene_id ) %>%
		summarize(
			UniqueSites = n()
		)

	full_gene_count_table <- data.frame( gene_id = unique_gene_names ) %>%
		merge( x = ., y = total_grna_per_gene, by = "gene_id", all = TRUE ) %>%
		merge( x = ., y = unique_grna_per_gene, by = "gene_id", all = TRUE ) %>%
		mutate( .data = ., gene_id = as.character( gene_id ) )

	full_gene_count_table$TotalSites[ is.na( full_gene_count_table$TotalSites ) ] = as.integer( 0 )
	full_gene_count_table$UniqueSites[ is.na( full_gene_count_table$UniqueSites ) ] = as.integer( 0 )

	########################
	# quick validity check #
	########################
	if ( any( full_gene_count_table$TotalSites < full_gene_count_table$UniqueSites, na.rm=TRUE ) )
	{
		stop( "More unique gRNAs than total gRNAs for some sites!!!" )
	}

	################################
	# write gene info per organism #
	################################
	# R_analysis/SPECIES_cutsPerGene.csv
	as.data.frame( full_gene_count_table ) %>%
		write_csv( x = ., path = cuts_per_gene_filepath )

	#####################################
	# Add gene annotation to CasX input #
	#####################################
	message( paste0( "Collapse genes for ", curr_species ) )
	
	collapsed_gene_annotations <- dplyr::select( olap_exon_table, uid, gene_id ) %>%
		partition( uid ) %>%
		summarize( exonOverlaps = paste( gene_id, collapse = "|" ) ) %>%
		collect()
	
	message( paste0( "Merging collapsed genes with motif data" ) )
	merge( x = all_casx_sites, y = collapsed_gene_annotations, by = 'uid', all.x = TRUE ) %>%
		mutate( .data = ., exonOverlaps = case_when(
			is.na( exonOverlaps ) ~ "",
			TRUE ~ exonOverlaps ) ) %>%
		dplyr::select( .data = ., uid, Contig, Start, End, Strand, Sequence, PAM, Target, exonOverlaps, ExactMatches, Unique ) %>%
		arrange( .data = ., Contig, Start, End, Strand ) %>%
		as.data.frame( . ) %>%
		write_csv( x = ., path = annotated_cuts_filepath )
}

################
# version info #
################
sessionInfo()
