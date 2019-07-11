#==================================================
#==================================================
#==================================================
# Libraries
library( bio3d )
library( stringr )
library( Interpol )
library( zoo )

#library( h2o )



#==================================================
#==================================================
#==================================================
# Functions

get_disorder_features <- function(df,type)
{
  df$V2 <- as.character(as.vector(df$V2))
  df$V3 <- as.character(as.vector(df$V3))
  
  #Fraction of residues involved in disorder region in a sequence
  #======================================================================
  fraction_seq_disordered <- sum(df$V3==type)/nrow(df);
  
  #======================================================================
  #Get frequency of each amino acid for all disordered regions in sequence
  list_amino <- c("A",
                  "R",
                  "N",
                  "D",
                  "C",
                  "E",
                  "Q",
                  "G",
                  "H",
                  "I",
                  "L",
                  "K",
                  "M",
                  "F",
                  "P",
                  "S",
                  "T",
                  "W",
                  "Y",
                  "V" )
  disorder_mono <- rep(0,length(list_amino))
  names(disorder_mono) <- list_amino
  
  disorder_info_from_seq <- table(df[df[,]$V3==type,]$V2)/nrow(df[df[,]$V3==type,])
  for (i in 1:length(disorder_info_from_seq))
  {
    amino_id <- names(disorder_info_from_seq[i])
    disorder_mono[amino_id] <- disorder_info_from_seq[amino_id]
  }
  
  #==============================================================
  #Get no of disordered regions in a sequence
  no_disordered_regions <- 0;
  gg <- which(df$V3==type)
  gg_list <- list()
  
  if (sum(gg)==0){
    no_disordered_regions <- 0;
  } else if (length(gg)==1)
  {
    no_disordered_regions <- 1
    gg_list[[1]] <- gg;
  } else 
  {
    count <- 1;
    prev_count <- 1;
    gg_list[[prev_count]] <- gg[1];
    for (j in 1:(length(gg)-1))
    {
      if ((gg[j+1]-gg[j])>1)
      {
        count <- count+1
      } 
      if (count==prev_count) {
        gg_list[[prev_count]] <- c(gg_list[[prev_count]],gg[j+1]);
      } else {
        gg_list[[prev_count+1]] <- gg[j+1];
      }
      prev_count=count;
    }
    no_disordered_regions <- count;
  }
  
  #==============================================================
  #Get binned information for no of disordered regions
  no_disorder_less_than_5 <- 0
  no_disorder_between_5_and_10 <- 0
  no_disorder_greater_10 <- 0;
  
  if (length(gg_list)>0)
  {
    for (i in 1:length(gg_list))
    {
      list_i <- gg_list[[i]];
      if (length(list_i)<5)
      {
        no_disorder_less_than_5 <- no_disorder_less_than_5+1;
      } else if (length(list_i)>=5 && length(list_i) <=10)
      {
        no_disorder_between_5_and_10 <- no_disorder_between_5_and_10+1;
      } else {
        no_disorder_greater_10 <- no_disorder_greater_10 + 1;
      }
    }
  }
  disorder_output <- list(as.numeric(disorder_mono),fraction_seq_disordered,no_disordered_regions,
                          no_disorder_less_than_5,no_disorder_between_5_and_10,
                          no_disorder_greater_10);
  return(disorder_output)
}


find_overlaps <- function(p,s) {
  gg <- gregexpr(paste0("(?=",p,")"),s,perl=TRUE)[[1]]
  if (length(gg)==1 && gg==-1){
    return (0)
  }
  else{
    return (length(gg))
  }
}

get_change_counts <- function(string_vec,s){
  pattern <- "(?="
  for (i in 1:length(string_vec))
  {
    pattern <- paste0(pattern,string_vec[i],"+")
  }
  pattern <- paste0(pattern,")");
  gg <- gregexpr(pattern,s,perl=TRUE)[[1]]
  if (length(gg)==1 && gg==-1){
    return (0)
  }
  else if (length(gg)==1 && gg!=-1){
    return(1)
  }
  else{
    count <- 1;
    for (j in 1:(length(gg)-1))
    {
      if ((gg[j+1]-gg[j])>1)
      {
        count <- count+1
      }
    }
    return(count)
  }
}
#==================================================
# Calculate features for test sequence
PaRSnIP.calc.features.RSAhydro.test <- function( vec.seq,
                                                 SCRATCH.path,
                                                 DISORDER.path,
                                                 output_prefix,
                                                 AA = unlist( strsplit("ACDEFGHIKLMNPQRSTVWY",split = "" ) ),
                                                 SS.3 = unlist( strsplit("CEH",split = "" ) ),
                                                 SS.8 = unlist( strsplit("BCEGHIST",split = "" ) ),
                                                 n.cores = 1 )
{
  #==================================================
  # Preprocess sequence
  # Step 1: Remove all inserts ("-")
  vec.seq <- vec.seq[ vec.seq != "-" ]
  # Step 2: Convert all non-standard amino acids to X
  vec.seq[ !( vec.seq %in% AA ) ] <- "X"
  
  #==================================================
  # Sequence length
  p <- length( vec.seq )
  var.log.seq.len <- log( p )
  
  #==================================================
  # Calculate molecular weight
  # df.mw <- data.frame( read.csv( "DT_MW.csv" ) )
  df.mw <- data.frame( cbind( c( "A",
                                 "R",
                                 "N",
                                 "D",
                                 "C",
                                 "E",
                                 "Q",
                                 "G",
                                 "H",
                                 "I",
                                 "L",
                                 "K",
                                 "M",
                                 "F",
                                 "P",
                                 "S",
                                 "T",
                                 "W",
                                 "Y",
                                 "V" ),
                              c( 89.1,
                                 174.2,
                                 132.1,
                                 133.1,
                                 121.2,
                                 147.1,
                                 146.2,
                                 75.1,
                                 155.2,
                                 131.2,
                                 131.2,
                                 146.2,
                                 149.2,
                                 165.2,
                                 115.1,
                                 105.1,
                                 119.1,
                                 204.2,
                                 181.2,
                                 117.1 ) ) )
  colnames( df.mw ) <- c( "AA",
                          "MW" )
  df.mw$AA <- as.vector( df.mw$AA )
  df.mw$MW <- as.numeric( as.vector( df.mw$MW ) )
  
                               
  vec.mw <- NULL
  for( i in 1:length( vec.seq ) )
  {
    if( nrow( df.mw[ df.mw$AA == vec.seq[ i ], ] ) > 0 )
    {
      vec.mw <- c( vec.mw,
                   df.mw[ df.mw$AA == vec.seq[ i ], ]$MW )
    }
  }
  # var.mw <- sum( vec.mw )
  var.mw <- log( sum( vec.mw ) )
  
  #==================================================
  # Frequency turn-forming residues
  vec.tfr <- 0
  for( i in 1:length( vec.seq ) )
  {
    if( vec.seq[ i ] %in% c( "N", "G", "P", "S" ) )
    {
      vec.tfr <- vec.tfr + 1  
    }
  }
  var.tfr <- vec.tfr / length( vec.seq )
  
  #==================================================
  # Calculate GRAVY index
  var.seq <- paste( vec.seq,
                    collapse = "" )
  var.gravy <- sum( unlist( Interpol::AAdescriptor( var.seq ) ) ) / length( vec.seq )
  
  #==================================================
  # Alipathic index
  vec.ali <- rep( 0,
                  4 )
  vec.ali[ 1 ] <- sum( vec.seq == "A" )
  vec.ali[ 2 ] <- sum( vec.seq == "V" )
  vec.ali[ 3 ] <- sum( vec.seq == "I" )
  vec.ali[ 4 ] <- sum( vec.seq == "L" )
  
  var.ali <- ( vec.ali[ 1 ] + 2.9 * vec.ali[ 2 ] + 3.9 * vec.ali[ 3 ] + 3.9 * vec.ali[ 4 ] ) / p
  
  #==================================================
  # Absolute charge
  vec.ch <- rep( 0,
                 4 )
  vec.ch[ 1 ] <- sum( vec.seq == "R" )
  vec.ch[ 2 ] <- sum( vec.seq == "K" )
  vec.ch[ 3 ] <- sum( vec.seq == "D" )
  vec.ch[ 4 ] <- sum( vec.seq == "E" )
  
  var.ch <- abs( ( ( vec.ch[ 1 ] + vec.ch[ 2 ] - vec.ch[ 3 ] - vec.ch[ 4 ] ) / p ) - 0.03 )
  
  #==================================================
  # Amino acid frequencies
  vec.AA.freq <- rep( NA,
                      length( AA ) )
  for( i in 1:length( AA ) )
  {
    vec.AA.freq[ i ] <- sum( vec.seq == AA[ i ] ) / p
  }

  #==================================================
  # Dipeptide frequencies
  df.dipep <- data.frame( expand.grid( AA,
                                       AA ) )
  vec.dipep <- apply( df.dipep,
                      1,
                      function( vec )
                      {
                        return( paste( as.vector( vec ),
                                       collapse = "" ) )
                      } )
  vec.dipep.freq <- rep( NA,
                         length( vec.dipep ) )
  for( i in 1:length( vec.dipep.freq ) )
  {
    vec.dipep.freq[ i ] <- find_overlaps(vec.dipep[ i ],
                                         paste(vec.seq,
                                               collapse=""))/(p-1)
  }

  #==================================================
  # Tripeptide frequencies
  df.tripep <- data.frame( expand.grid( AA,
                                        AA,
                                        AA ) )
  vec.tripep <- apply( df.tripep,
                       1,
                       function( vec )
                       {
                         return( paste( as.vector( vec ),
                                        collapse = "" ) )
                       } )
  vec.tripep.freq <- rep( NA,
                          length( vec.tripep ) )
  for( i in 1:length( vec.tripep.freq ) )
  {
    vec.tripep.freq[ i ] <- find_overlaps(vec.tripep[ i ], 
                                          paste( vec.seq, 
                                                 collapse = ""))/(p-2)
  }
  # 
  #==================================================
  #==================================================
  #==================================================
  # SCRATCH features

  #==================================================
  # Run SCRATCH
  run.SCRATCH( vec.seq,
               SCRATCH.path,
               output_prefix,
               n.cores )


  #==================================================
  # Run DISORDER
  run.DISORDER( vec.seq,
                DISORDER.path,
                output_prefix )


  #==================================================
  # Disorder information
  df_diso <- read.table(paste0(output_prefix,".diso"),header=FALSE)
  df_pbdat <- read.table(paste0(output_prefix,".pbdat"),header=FALSE)

  disorder_info <- get_disorder_features(df_diso,"*")
  pb_disorder_info <- get_disorder_features(df_pbdat,"^")

  #==================================================
  # 3-state secondary structure classification
  file.ss <- paste( output_prefix,
                    ".ss",
                    sep = "" )
  vec.ss <- as.vector( read.fasta( file.ss )$ali )
  vec.ss.freq <- NULL
  for( i in 1:length( SS.3 ) )
  {
    vec.ss.freq <- c( vec.ss.freq,
                      sum( vec.ss == SS.3[ i ] ) )
  }
  vec.ss.freq <- vec.ss.freq / p

  #================================================
  # Di 3-state change information
   df.ss.dipep <- data.frame( expand.grid( SS.3,
                                          SS.3)
                                        )
   df.ss.dipep$Var1 <- as.character(as.vector(df.ss.dipep$Var1));
   df.ss.dipep$Var2 <- as.character(as.vector(df.ss.dipep$Var2));
   vec.ss.dipep.freq <- NULL
   for (i in 1:nrow(df.ss.dipep))
   {
    #if (df.ss.dipep$Var1[i]!=df.ss.dipep$Var2[i])
    {
      string <- as.character(as.vector(df.ss.dipep[i,]));
      vec.ss.dipep.freq <- c(vec.ss.dipep.freq,
                             get_change_counts(string,paste( as.vector( vec.ss ),
                                                             collapse = "" )))
    }
   }
   vec.ss.dipep.freq <- vec.ss.dipep.freq/sum(vec.ss.dipep.freq);

  #===================================================================
  # Tri 3-state change information
  df.ss.tripep <- data.frame( expand.grid( SS.3,
                                           SS.3,
                                           SS.3)
  )
  df.ss.tripep$Var1 <- as.character(as.vector(df.ss.tripep$Var1));
  df.ss.tripep$Var2 <- as.character(as.vector(df.ss.tripep$Var2));
  df.ss.tripep$Var3 <- as.character(as.vector(df.ss.tripep$Var3));
  vec.ss.tripep.freq <- NULL
  for (i in 1:nrow(df.ss.tripep))
  {
    #if (df.ss.tripep$Var1[i]!=df.ss.tripep$Var2[i] && df.ss.tripep$Var2[i]!=df.ss.tripep$Var3[i])
    {
      string <- as.character(as.vector(df.ss.tripep[i,]));
      vec.ss.tripep.freq <- c(vec.ss.tripep.freq,
                             get_change_counts(string,paste( as.vector( vec.ss ),
                                                             collapse = "" )))
    }
  }
  vec.ss.tripep.freq <- vec.ss.tripep.freq/sum(vec.ss.tripep.freq);

  #==================================================
  # 8-state secondary structure classification
  file.ss8 <- paste( output_prefix,
                     ".ss8",
                     sep = "" )
  vec.ss8 <- as.vector( read.fasta( file.ss8 )$ali )
  vec.ss8.freq <- NULL
  for( i in 1:length( SS.8 ) )
  {
    vec.ss8.freq <- c( vec.ss8.freq,
                       sum( vec.ss8 == SS.8[ i ] ) )
  }
  vec.ss8.freq <- vec.ss8.freq / p

  #================================================
  # Di 8-state change information
  df.ss8.dipep <- data.frame( expand.grid( SS.8,
                                          SS.8)
  )
  df.ss8.dipep$Var1 <- as.character(as.vector(df.ss8.dipep$Var1));
  df.ss8.dipep$Var2 <- as.character(as.vector(df.ss8.dipep$Var2));
  vec.ss8.dipep.freq <- NULL
  for (i in 1:nrow(df.ss8.dipep))
  {
    #if (df.ss8.dipep$Var1[i]!=df.ss8.dipep$Var2[i])
    {
      #print(df.ss8.dipep[i,])
      string <- as.character(as.vector(df.ss8.dipep[i,]));
      vec.ss8.dipep.freq <- c(vec.ss8.dipep.freq,
                             get_change_counts(string,paste( as.vector( vec.ss8 ),
                                                             collapse = "" )))
    }
  }
  vec.ss8.dipep.freq <- vec.ss8.dipep.freq/sum(vec.ss8.dipep.freq);

  #===================================================================
  # Tri 3-state change information
  df.ss8.tripep <- data.frame( expand.grid( SS.8,
                                           SS.8,
                                           SS.8)
  )
  df.ss8.tripep$Var1 <- as.character(as.vector(df.ss8.tripep$Var1));
  df.ss8.tripep$Var2 <- as.character(as.vector(df.ss8.tripep$Var2));
  df.ss8.tripep$Var3 <- as.character(as.vector(df.ss8.tripep$Var3));
  vec.ss8.tripep.freq <- NULL
  for (i in 1:nrow(df.ss8.tripep))
  {
    #if (df.ss8.tripep$Var1[i]!=df.ss8.tripep$Var2[i] && df.ss8.tripep$Var2[i]!=df.ss8.tripep$Var3[i])
    {
      string <- as.character(as.vector(df.ss8.tripep[i,]));
      vec.ss8.tripep.freq <- c(vec.ss8.tripep.freq,
                              get_change_counts(string,paste( as.vector( vec.ss8 ),
                                                              collapse = "" )))
    }
  }
  vec.ss8.tripep.freq <- vec.ss8.tripep.freq/sum(vec.ss8.tripep.freq);

  #==================================================
  # Solvent accessibility prediction at 0%-95% thresholds
  file.acc.20 <- paste( output_prefix,
                        ".acc20",
                        sep = "" )
  vec.acc.20.raw <- scan( file.acc.20,
                          what = "character" )
  vec.acc.20 <- as.numeric( vec.acc.20.raw[ 2:length( vec.acc.20.raw ) ] )
  vec.thresh <- seq( 0, 95, 5 )
  ## vec.acc.20.final <- NULL
  ## for( i in 1:length( vec.thresh ) )
  ## {
  ##   vec.acc.20.final <- c( vec.acc.20.final,
  ##                          sum( vec.acc.20 == vec.thresh[ i ] ) / p )
  ## }
  vec.acc.20.final <- sum( vec.acc.20 == vec.thresh[ 1 ] ) / p
  for( i in 2:length( vec.thresh ) )
  {
    vec.acc.20.final <- c( vec.acc.20.final,
                           sum( vec.acc.20 >= vec.thresh[ i ] ) / p )
  }
  #==================================================
  # Solvent accessibility prediction at 0%-95% thresholds coupled with average hydrophobicity
  ind.thresh <- which( vec.acc.20 == vec.thresh[ 1 ] )
  if( length( ind.thresh ) == 0 )
  {
    vec.rsa.hydro <- 0
  } else
  {
    vec.rsa.hydro <- ( sum( vec.acc.20 == vec.thresh[ 1 ] ) / p ) *
      mean( unlist( Interpol::AAdescriptor( vec.seq[ ind.thresh ] ) ) )
  }
  for( i in 2:length( vec.thresh ) )
  {
    ind.thresh <- which( vec.acc.20 >= vec.thresh[ i ] )
    if( length( ind.thresh ) == 0 )
    {
      vec.rsa.hydro <- c( vec.rsa.hydro,
                          0 )
    } else
    {
      vec.rsa.hydro <- c( vec.rsa.hydro,
                          ( sum( vec.acc.20 >= vec.thresh[ i ] ) / p ) *
                            mean( unlist( Interpol::AAdescriptor( vec.seq[ ind.thresh ] ) ) ) )
    }
  }
  
  
  #==================================================
  #==================================================
  #==================================================
  # Return feature vector
  vec.features <- c( var.log.seq.len,
                     var.mw,
                     var.tfr,
                     var.gravy,
                     var.ali,
                     var.ch,
                     vec.ss.freq,
                     vec.ss.dipep.freq,
                     vec.ss.tripep.freq,
                     vec.ss8.freq,
                     vec.ss8.dipep.freq,
                     vec.ss8.tripep.freq,
                     vec.acc.20.final,
                     vec.rsa.hydro,
                     disorder_info[[1]],
                     disorder_info[[2]],
                     disorder_info[[3]],
                     disorder_info[[4]],
                     disorder_info[[5]],
                     disorder_info[[6]],
                     pb_disorder_info[[1]],
                     pb_disorder_info[[2]],
                     pb_disorder_info[[3]],
                     pb_disorder_info[[4]],
                     pb_disorder_info[[5]],
                     pb_disorder_info[[6]],
                     vec.AA.freq,
                     vec.dipep.freq,
                     vec.tripep.freq
                     )
  
  return( vec.features )
}



#==================================================
# Run SCRATCH
run.SCRATCH <- function( vec.seq,
                         SCRATCH.path,
                         output_prefix,
                         n.cores = 1 )
{
  # Sequence file (FASTA format)
  file.fasta <- paste( output_prefix,
                       ".fasta",
                       sep = "" )
  # Save sequence in a tmp fasta file
  write.fasta( ids = output_prefix,
               seqs = vec.seq,
               file = file.fasta )
  
  # Start running SCRATCH
  print( system2( SCRATCH.path,
                  args = c( file.fasta,
                            output_prefix,
                            n.cores ) ) )
}


#==================================================
# Run DISORDER
run.DISORDER <- function( vec.seq,
                         DISORDER.path,
                         output_prefix
                        )
{
  # Sequence file (FASTA format)
  file.fasta <- paste( output_prefix,
                       ".fasta",
                       sep = "" )
  # Save sequence in a tmp fasta file
  write.fasta( ids = output_prefix,
               seqs = vec.seq,
               file = file.fasta )
  
  # Start running DISORDER
  print( system2( DISORDER.path,
                  args = c(file.fasta
                           ) ) ) 
}

#==================================================
# Convert GBM output into probabilities
get.probabilities <- function( vec.raw,
                               threshold )
{
  # get_class_one_predictions <- as.numeric(as.vector(y.predict.raw$p1))
  
  # Convert the raw score into a probability
  y.predict.prob <- rep( 0,
                         length( vec.raw ) )
  for (i in 1:length( vec.raw ) )
  {
    if (vec.raw[ i ] >= threshold )
    {
      y.predict.prob[ i ] <- 0.5 + 0.5 * ( ( vec.raw[ i ] - threshold ) / ( 1 - threshold ) )
    }
    else
    {
      y.predict.prob[ i ] <- 0.5 - 0.5 * ( ( threshold - vec.raw[ i ] ) / threshold )
    }
  }
  return( y.predict.prob )
}


#==================================================

#==================================================
#==================================================
#==================================================
