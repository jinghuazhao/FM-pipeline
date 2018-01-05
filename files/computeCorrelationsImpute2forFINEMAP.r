#---------------------------------------------#
# Distributed under the Boost Software        #
# License, Version 1.0. (See copy at          #
# http://www.boost.org/LICENSE_1_0.txt)       #
#---------------------------------------------#
#                                             #
# R script to compute Pearson correlations    #
# between variants                            #
#                                             #
# Input: Impute2 files                        #
#                                             #
# Output:                                     #
# 1) Pairwise correlations are saved in an    #
#    ASCII file                               #
# 2) Mapping file with variant information    #
#    is used to identify each variant pair    #
#                                             #
# Command-line arguments                      #
# ----------------------                      #
#                                             #
# Arg1  : Filename of the info file (_info)   #
# Arg2  : Filename of the gen file (.gen.gz)  #
# Arg3  : Chromosome number of the genomic    #
#         for which LD information is to be   #
#         computed                            #
# Arg4  : Start position of the genomic       #
#         region for which LD information     #
#         is to be computed                   #
# Arg5  : End position of the genomic region  #
# Arg6  : Variants with minor allele freq     #
#         above this threshold will be        #
#         included in the computation         #
# Arg7 : Variants with imputation quality     #
#         above this threshold will be        #
#         included in the computation         #
# Arg8 : Basename of the binary output file   #
#         (no file extension)                 #
# Arg9 : Number of threads                    #
#                                             #
# Usage : Rscript \                           #
#         --vanilla \                         #
# computeCorrelationsImpute2forFINEMAP.r \    #
#         myfile_info \                       #
#         myfile.gen.gz \                     #
#         10 \                                #
#         1939202 \                           #
#         2945274 \                           #
#         0.05 \                              #
#         0.9 \                               #
#         outputFile                          #
#         n                                   #
#                                             #
#---------------------------------------------#

library( parallel )

ptm <- proc.time()

# Parse command-line arguments
args <- commandArgs( trailingOnly = T )
if( length( args ) != 9 ) {
	stop( "Nine arguments must be supplied!", call. = FALSE )
} else {
	for( ii in c( 4, 5 ) ) {
		number <- as.numeric( args[ ii ] )
		if( is.na( number ) | number <= 0.0 | number != ceiling( number ) ) { 
			stop( sprintf( "Argument %d must be a positive integer value!", ii ), call. = FALSE )
		}
	}
	for( ii in 6 ) {
		number <- as.numeric( args[ ii ] )
		if( is.na( number ) | number < 0.0 | number >= 0.5 ) { 
			stop( sprintf( "Argument %d must be a floating-point value between 0.0 and 0.5!", ii ), call. = FALSE )
		}
	}
	for( ii in 7 ) {
		number <- as.numeric( args[ ii ] )
		if( is.na( number ) | number < 0.0 | number >= 1.0 ) { 
			stop( sprintf( "Argument %d must be a floating-point value between 0.0 and 1.0!", ii ), call. = FALSE )
		}
	}
	if( as.numeric( args[ 4 ] ) >= as.numeric( args[ 5 ] ) ) {
		stop( "Argument 7 (start position) must be smaller than argument 8 (end position)!", call. = FALSE )
	}
	n_thread <- as.numeric( args[ 9 ] )
	if( is.na( n_thread ) | n_thread <= 0 | n_thread != ceiling( n_thread ) ) {
		stop( "Argument 9 (number of threads) must be a positive integer value!", call. = FALSE )
	}
}

#-----------#
# Filenames #
#-----------#

# .info file
info_filename <- args[ 1 ]
if( !file.exists( info_filename ) ) {
	stop( sprintf( "File %s not exist!", info_filename ), call. = FALSE )
}

# .gen.gz file
gen_filename <- args[ 2 ]
if( !file.exists( gen_filename ) ) {
	stop( sprintf( "File %s not exist!", gen_filename ), call. = FALSE )
}

# Mapping output file
map_filename <- sprintf( "%s.map", args[ 8 ] )

# Binary file with Pearson correlations
ld_filename <- sprintf( "%s.ld.gz", args[ 8 ] )

#-------------------#
# Columns           #
# Chromosome        #
# Genomic positions #
# Filter            #
#-------------------#

chromosome <- args[ 3 ]
start_position <- as.numeric( args[ 4 ] )
end_position <- as.numeric( args[ 5 ] )
maf_filter <- as.numeric( args[ 6 ] )
imputation_filter <- as.numeric( args[ 7 ] )

#-----------------#
# Read .info file #
#-----------------#

info_data <- read.table( info_filename, header = TRUE, as.is = TRUE )

# Check that column names
# are as expected
if( any( is.na( match( c( "snp_id", "position", "info" ), colnames( info_data ) ) ) ) ) {
	stop( sprintf( "Expecting column names 'snp_id', 'position', and 'info' in %s!", info_filename ), call. = FALSE )
}

# Check that all genomic positions
# are positive integer values
if( any( sapply( info_data[, "position" ], function( x ) {
	number <- as.numeric( x )
	is.na( number ) | number <= 0.0 | number != ceiling( number ) 
} ) ) ) {
	stop( "The genomic positions must be positive integer values!", call. = FALSE )
}

# Check that all chromosome numbers
if( any( sapply( info_data[, "snp_id" ], function( x ) {
	! x %in% c( "---", chromosome, paste( 0, chromosome, sep = '') ) 
} ) ) ) {
	stop( "The chromsome numbers must be correct in snp_id of info file!", call. = FALSE )
}

# Check that genomic positions
# are sorted in increasing order
if( 
	!identical( 
		1 : nrow( info_data ), 
		order( info_data[, "position" ] ) 
	) 
) {
	stop( "The genomic positions must be sorted in increasing order!", call. = FALSE )
}

# Get indexes of variants within the
# specified genomic range
indexes_variants_to_include <- c( 
	info_data$position >= start_position & 
	info_data$position <= end_position & 
	info_data$info >= imputation_filter
)
m0_variants <- sum( indexes_variants_to_include )
if( m0_variants == 0 ) {
	stop( "No variants were selected!", call. = FALSE )
}
info_data <- info_data[ indexes_variants_to_include, ]

#-------------------#
# Read .gen.gz file #
#-------------------#

# Read first line only
file_handle <- gzfile( gen_filename, 'r' )
first_line <- strsplit( readLines( file_handle, n = 1 ), split = " " )[[ 1 ]]
close( file_handle )

# Save number of individuals
n_samples <- ( length( first_line ) - 5 ) / 3
# n_samples <- ( length( first_line ) - 6 ) / 3#For files from qctool
if ( n_samples != ceiling( n_samples ) ) {
	stop( "The gen.gz file is not correct for sample size!", call. = FALSE )
}

what <- vector( "list", length( first_line ) )
what[[ 3 ]] <- numeric()
what[[ 4 ]] <- what[[ 5 ]] <- character()
for ( ii in 6:length( first_line ) ){
	what[[ ii ]] <- numeric()
}

# Read gen of variants within
# the specified genomic range
start_position_index <- min( which( indexes_variants_to_include ) )
m1_variants <- max( which( indexes_variants_to_include ) ) - start_position_index + 1
file_handle <- gzfile( gen_filename, 'r' )
gen_data <- scan( file_handle, what, skip = start_position_index - 1, nlines = m1_variants, quiet = TRUE )
close( file_handle )

# Check that there are gens for as
# many variants as there are in the
# .info file
indexes_variants_to_include <- which( indexes_variants_to_include ) - start_position_index + 1
if( any( gen_data[[ 3 ]][ indexes_variants_to_include ] != info_data$position ) ) {
	stop( 
		sprintf( 
			"The positions of variants in '%s' and '%s' are not the same!", 
			gen_filename,
			info_filename 
		),
		call. = FALSE 
	)
}

ref <- gen_data[[ 4 ]][ indexes_variants_to_include ]
alt <- gen_data[[ 5 ]][ indexes_variants_to_include ]
gen_data[ 3 ] <- gen_data[ 4 ] <- gen_data[ 5 ] <- NULL
gen_data <- do.call( cbind,  gen_data )
dosages <- ( gen_data[ , ( 1:n_samples ) * 3 - 1 ] + 2 * gen_data[ , ( 1:n_samples ) * 3 ] ) / ( gen_data[ , ( 1:n_samples ) * 3 - 2 ] + gen_data[ , ( 1:n_samples ) * 3 - 1 ] + gen_data[ , ( 1:n_samples ) * 3 ] )
rm( gen_data )
dosages <- t( dosages[ indexes_variants_to_include, ] )

# Check that there are no missing dosages
missdosage <- is.na( dosages )
if( sum( missdosage ) > 0 ) {
	cat( sprintf( "There are %d missing genotypes in '%s' and they will be mean-imputed!\n", sum( missdosage ), gen_filename ) )
}

maf <- alt_af <- colMeans( dosages, na.rm = TRUE ) / 2
for( ii in which( colSums( missdosage ) > 0 ) ){ dosages[ missdosage[ , ii ], ii ] <- alt_af[ ii ] * 2 }

maf[ maf > 0.5 ] <- 1 - maf[ maf > 0.5 ]
maf_keep <- maf >= maf_filter
dosages <- dosages[ , maf_keep ]

# Save number of variants
m_variants <- ncol( dosages )

# Save number of variant pairs
m_pairs <- 0.5 * ( m_variants + m_variants ^ 2 ) - m_variants

#--------------------#
# Write mapping file #
#--------------------#

# Prepare mapping to identify
# each variant pair
map_data <- cbind( 
	Chr = chromosome, 
	Position = info_data$position[ maf_keep ],
	Ref = ref[ maf_keep ],
	Alt = alt[ maf_keep ],
	Alt_af = round( alt_af[ maf_keep ], 4 ),
	MAF = round( maf[ maf_keep ], 4 ),
	Info = round( info_data$info[ maf_keep ], 3 )
)

# Write mapping file
file_handle <- file( map_filename, 'wt' )
writeLines( sprintf( "## Info               : %s", info_filename ), file_handle )
writeLines( sprintf( "## Gen                : %s", gen_filename ), file_handle )
writeLines( sprintf( "## Correlations       : %s", ld_filename ), file_handle )
writeLines( sprintf( "## Chromosome         : %s", chromosome ), file_handle )
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
file_handle <- gzfile( ld_filename, 'w' )
chunks <- unique( c( seq(1, m_variants + 1, 1000), m_variants + 1 ) )
#parSapply( cl, 2:length( chunks ), mycor)
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
