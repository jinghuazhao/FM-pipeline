#------------------------------------------------#
# Distributed under the Boost Software           #
# License, Version 1.0. (See copy at             #
# http://www.boost.org/LICENSE_1_0.txt)          #
#------------------------------------------------#
#                                                #
# R script to compute Pearson correlations       #
# between variants                               #
#                                                #
# Input: Minimac 1/2 files                       #
#                                                #
# Output:                                        #
# 1) Pairwise correlations are saved in an       #
#    ASCII file                                  #
# 2) Mapping file with variant information       #
#    is used to identify each variant pair       #
#                                                #
# Command-line arguments                         #
# ----------------------                         #
#                                                #
# Arg1  : Filename of the info file (.info)      #
# Arg2  : Filename of the dosage file (.dose.gz) #
# Arg3  : Filename of the file with genomic      #
#         positions and chromosome number        #
#         Can be the same as info file           #
#         Header ("SNP', "Position", "Chr")      #
#         is expected                            #
# Arg4  : String indicating the column name      #
#         in the file specified with Arg3        #
#         that contains the genomic positions    #
# Arg5  : String indicating the column name      #
#         in the file specified with Arg3        #
#         that contains the chromosome number    #
# Arg6  : Chromosome number of the genomic       #
#         for which LD information is to be      #
#         computed                               #
# Arg7  : Start position of the genomic          #
#         region for which LD information        #
#         is to be computed                      #
# Arg8  : End position of the genomic region     #
# Arg9  : Variants with minor allele freq        #
#         above this threshold will be           #
#         included in the computation            #
# Arg10 : Variants with imputation quality       #
#         above this threshold will be           #
#         included in the computation            #
# Arg11 : Basename of the binary output file     #
#         (no file extension)                    #
# Arg12 : Number of threads                       #
#                                                #
# Usage : Rscript \                              #
#         --vanilla \                            #
#         computeCorrelationsMinimac12.r \       #
#         infofile.info \                        #
#         dosagefile.dose.gz \                   #
#         positionsFile \                        #
#         Position \                             #
#         Chromosome \                           #
#         10 \                                   #
#         1939202 \                              #
#         2945274 \                              #
#         0.05 \                                 #
#         0.9 \                                  #
#         outputFile \                           #
#         n                                      #
#                                                #
# Contact : christian.benner@helsinki.fi         #
#------------------------------------------------#

library( parallel )

ptm <- proc.time()

# Parse command-line arguments
args <- commandArgs( trailingOnly = T )
if( length( args ) != 12 ) {
	stop( "Ten arguments must be supplied!", call. = FALSE )
} else {
	for( ii in seq( 7, 8 ) ) {
		number <- as.numeric( args[ ii ] )
		if( is.na( number ) | number <= 0.0 | number != ceiling( number ) ) { 
			stop( sprintf( "Argument %d must be a positive integer value!", ii ), call. = FALSE )
		}
	}
	for( ii in 9 ) {
		number <- as.numeric( args[ ii ] )
		if( is.na( number ) | number < 0.0 | number >= 0.5 ) { 
			stop( sprintf( "Argument %d must be a floating-point value between 0.0 and 0.5!", ii ), call. = FALSE )
		}
	}
	for( ii in 10 ) {
		number <- as.numeric( args[ ii ] )
		if( is.na( number ) | number < 0.0 | number >= 1.0 ) { 
			stop( sprintf( "Argument %d must be a floating-point value between 0.0 and 1.0!", ii ), call. = FALSE )
		}
	}
	if( as.numeric( args[ 7 ] ) >= as.numeric( args[ 8 ] ) ) {
		stop( "Argument 7 (start position) must be smaller than argument 8 (end position)!", call. = FALSE )
	}
	n_thread <- as.numeric( args[ 12 ] )
	if( is.na( n_thread ) | n_thread <= 0 | n_thread != ceiling( n_thread ) ) {
		stop( "Argument 12 (number of threads) must be a positive integer value!", call. = FALSE )
	}
}

#-----------#
# Filenames #
#-----------#

# .info file
info_filename <- args[ 1 ]
if( !file.exists( info_filename ) ) {
	stop( sprintf( "File %s does not exist!", info_filename ), call. = FALSE )
}

# .dose file
dose_filename <- args[ 2 ]
if( !file.exists( dose_filename ) ) {
	stop( sprintf( "File %s does not exist!", dose_filename ), call. = FALSE )
}

# file with variant IDs, genomic positions and chromosome number
positions_filename <- args[ 3 ]
if( !file.exists( positions_filename ) ) {
	stop( sprintf( "File %s does not exist!", positions_filename ), call. = FALSE )
}

# Mapping output file
map_filename <- sprintf( "%s.map", args[ 11 ] )

# Binary file with Pearson correlations
ld_filename <- sprintf( "%s.ld.gz", args[ 11 ] )

#-------------------#
# Columns           #
# Chromosome        #
# Genomic positions #
# Filter            #
#-------------------#

column_positions <- args[ 4 ]
column_chromosome <- args[ 5 ]
chromosome <- as.numeric( args[ 6 ] )
start_position <- as.numeric( args[ 7 ] )
end_position <- as.numeric( args[ 8 ] )
maf_filter <- as.numeric( args[ 9 ] )
imputation_filter <- as.numeric( args[ 10 ] )

#-----------------#
# Read .info file #
#-----------------#

info_data <- read.table( info_filename, header = TRUE, as.is = TRUE )

# Check that column names
# are as expected
if( any( is.na( match( c( "SNP", "Al1", "Al2", "MAF", "Rsq" ), colnames( info_data ) ) ) ) ) {
	stop( sprintf( "Expecting column names 'SNP', 'Al1', 'Al2', 'MAF' and 'Rsq' in %s!", info_filename ), call. = FALSE )
}

#---------------#
# Read position #
# file          #
#---------------#

if( info_filename != positions_filename ) {
	positions_data <- read.table( positions_filename, header = TRUE, as.is = TRUE )
	if( any( is.na( match( c( "SNP", "Position", "Chr" ), colnames( positions_data ) ) ) ) ) {
		stop( sprintf( "Expecting column names 'SNP', 'Position' and 'Chr' in %s!", positions_filename ), call. = FALSE )
	}
} else {
	positions_data <- info_data
}

# Check that there are columns with genomic positions
# and chromosome numbers as specified with command-line
# arguments
if( !( column_positions %in% colnames( positions_data ) ) ) {
	stop( sprintf( "Cannot find column %s in %s!", column_positions, positions_filename ), call. = FALSE )
} else if( !( column_chromosome %in% colnames( positions_data ) ) ) {
	stop( sprintf( "Cannot find column %s in %s!", column_chromosome, positions_filename ), call. = FALSE )
}

# Check that there are genomic position for as
# many variants as there are in the .info file
if( nrow( positions_data ) != nrow( info_data ) ) {
	stop( 
		sprintf( 
			"The number of variants in '%s' and '%s' are not the same!", 
			positions_filename,
			info_filename 
		),
		call. = FALSE 
	)
}

# Check that all genomic positions
# are positive integer values
if( any( sapply( positions_data[, column_positions ], function( x ) {
	number <- as.numeric( x )
	is.na( number ) | number <= 0.0 | number != ceiling( number ) 
} ) ) ) {
	stop( "The genomic positions must be positive integer values!", call. = FALSE )
}

# Check that all chromosome numbers
# are positive integer values
if( any( sapply( positions_data[, column_chromosome ], function( x ) {
	number <- as.numeric( x )
	is.na( number ) | number <= 0.0 | number != ceiling( number ) 
} ) ) ) {
	stop( "The chromsome numbers must be positive integer values!", call. = FALSE )
}

# Check that genomic positions
# are sorted in increasing order
if( 
	!identical( 
		1 : nrow( positions_data ), 
		order( positions_data[, column_chromosome ], positions_data[, column_positions ] ) 
	) 
) {
	stop( "The genomic positions must be sorted in increasing order!", call. = FALSE )
}

# Get indexes of variants within the
# specified genomic range
indexes_variants_to_include <- which( 
	positions_data[, column_positions ] >= start_position & 
	positions_data[, column_positions ] <= end_position & 
	positions_data[, column_chromosome ] == chromosome &
	info_data[, 'MAF' ] >= maf_filter &
	info_data[, 'Rsq' ] >= imputation_filter
)
if( length( indexes_variants_to_include ) == 0 ) {
	stop( "No variants were selected!", call. = FALSE )
}

#-----------------#
# Read .dose file #
#-----------------#

# Read first line only
file_handle <- gzfile( dose_filename, 'r' )
first_line <- strsplit( readLines( file_handle, n = 1 ), split = " |\t" )[[ 1 ]]
close( file_handle )

# Check that the second value is
# equal to "DOSE"
if( grep( "DOSE", first_line[ 2 ] ) != 1) {
	stop( sprintf( "The second column in '%s' must contain 'DOSE'!", dose_filename ), call. = FALSE )
}

# Check that there are dosages for as
# many variants as there are in the
# .info file
if( ( length( first_line ) - 2 ) != nrow( info_data ) ) {
	stop( 
		sprintf( 
			"The number of variants in '%s' and '%s' are not the same!", 
			dose_filename,
			info_filename 
		),
		call. = FALSE 
	)
}

# Specify type of columns to
# be read with scan()
what <- vector( "list", length( first_line ) )
for( ii in ( 2 + indexes_variants_to_include ) ) {
	what[[ ii ]] <- numeric()
}

# Read dosage of variants within
# the specified genomic range
file_handle <- gzfile( dose_filename, 'r' )
dosages <- do.call( "cbind", scan( file_handle, what, quiet = TRUE ) )
close( file_handle )

missdosage <- is.na( dosages )
if( sum( missdosage ) > 0 ) {
	cat( sprintf( "There are %d missing genotypes in '%s' and they will be mean-imputed!\n", sum( missdosage ), dose_filename ) )
}
alt_af <- colMeans( dosages, na.rm = TRUE ) / 2
for( ii in which( colSums( missdosage ) > 0 ) ){ dosages[ missdosage[ , ii ], ii ] <- alt_af[ ii ] * 2 }

# Save number of variants
m_variants <- ncol( dosages )

# Save number of individuals
n_samples <- nrow( dosages )

# Save number of variant pairs
m_pairs <- 0.5 * ( m_variants + m_variants ^ 2 ) - m_variants

#--------------------#
# Write mapping file #
#--------------------#

# Prepare mapping to identify
# each variant pair
map_data <- cbind( 
	SNP = info_data[ indexes_variants_to_include, "SNP" ], 
	Chr = positions_data[ indexes_variants_to_include, column_chromosome ], 
	Position = positions_data[ indexes_variants_to_include, column_positions ], 
	info_data[ indexes_variants_to_include, c( "Al1", "Al2", "MAF", "Rsq" ) ],
	stringsAsFactors = FALSE
)

# Write mapping file
file_handle <- file( map_filename, 'wt' )
writeLines( sprintf( "## Info               : %s", info_filename ), file_handle )
writeLines( sprintf( "## Dose               : %s", dose_filename ), file_handle )
writeLines( sprintf( "## Positions          : %s", positions_filename ), file_handle )
writeLines( sprintf( "## Correlations       : %s", ld_filename ), file_handle )
writeLines( sprintf( "## Chromosome         : %d", chromosome ), file_handle )
writeLines( sprintf( "## Genomic range      : %d-%d", start_position, end_position ), file_handle )
writeLines( sprintf( "## Minor allele freq  : %.2f", maf_filter ), file_handle )
writeLines( sprintf( "## Imputation quality : %.2f", imputation_filter ), file_handle )
writeLines( sprintf( "## Number of variants : %d", m_variants ), file_handle )
writeLines( sprintf( "## Missing variants   : %d", sum( missdosage ) ), file_handle )
writeLines( sprintf( "## Number of pairs    : %s", m_pairs ), file_handle )
writeLines( sprintf( "## Number of samples  : %d", n_samples ), file_handle )
write.table( 
	map_data, 
	file_handle, 
	row.names = FALSE, 
	quote = FALSE, 
	sep = "\t",
)
close( file_handle )

cl <- makeCluster( getOption( "cl.cores", n_thread ) )
file_handle <- gzfile( ld_filename, 'w', compression = 6 )
chunks <- unique( c( seq(1, m_variants + 1, 10), m_variants + 1 ) )
for ( i in 2:length( chunks ) ){
# Compute Pearson correlations
	chunk_len <- chunks[ i ] - chunks[ i - 1 ]
	R <- matrix( NA, chunk_len, m_variants )
	R[ , 1:( chunks[ i ] - 1 ) ] <- round( parCapply( cl, dosages[ , 1:( chunks[ i ] - 1 ) ], get( "cor" ), dosages[ , chunks[ i - 1 ]:( chunks[ i ] - 1 ) ] ), 8 )
	R[ which( abs(R)<0.001, arr.ind=T ) ] <- 0
	naR <- cbind( matrix( FALSE, chunk_len, chunks[ i - 1 ] -1 ), upper.tri( matrix( NA, chunk_len, chunk_len ), diag = TRUE ), matrix( FALSE, chunk_len, m_variants +1 - chunks[ i ] ) )
	R[ naR ] <- NA
#-----------------#
# Write LD matrix #
# as ASCII file   #
#-----------------#
	write.table( R, file = file_handle, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
close( file_handle )
stopCluster( cl )

proc.time() - ptm
