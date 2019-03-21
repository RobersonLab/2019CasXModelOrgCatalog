library( here )
library( tidyverse )

samples <- read_tsv( file = here::here( "names_experiments_table.txt" ) ) %>%
  filter( mouse == TRUE )

mouse_samples <- pull( .data = samples, filename )	

upset_tibble <- tibble()

for ( idx in 1:length( mouse_samples ) ) {
  curr_name <- samples$filename[ idx ]
  curr_fname <- paste0( curr_name, "_targets.csv.gz" )
  curr_fpath <- here::here( 'R_analysis', curr_fname )
  
  temp <- read_csv( file = curr_fpath, progress = FALSE, col_types = 'c' ) %>%
    mutate( .data = ., species = 1 )
  colnames( temp )[2] = curr_name
  
  if ( sum( dim( upset_tibble ) ) > 0 ) {
    upset_tibble <- merge( upset_tibble, temp, by = "Target", all = TRUE )
  } else {
    upset_tibble = temp
  }
}

col_idx <- which( colnames( upset_tibble ) != "Target" )

for ( idx in 1:length( col_idx ) ) {
	curr_idx <- col_idx[ idx ]
	NA_idx <- which( is.na( upset_tibble[ , curr_idx ] ) )
	
	if ( length( is.na( NA_idx ) ) > 0 ) {
		upset_tibble[ NA_idx, curr_idx ] <- 0
	}
}

write_csv( x = upset_tibble, path = here::here( 'R_analysis', 'mouse_target_upset.csv' ) )
