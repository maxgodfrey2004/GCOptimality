args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  cat("ERROR: usage is <R binary> gen_aadata.r <output_filename>\n")
}
stopifnot(length(args) == 1)

if (! requireNamespace('seqinr', quietly=TRUE)) {
  install.packages('seqinr')
}

data(aaindex, package='seqinr')

frame <- data.frame(list(
  aaindex[[2]]$I,    # Hydrophobicity
  aaindex[[8]]$I,    # Flexibility
  aaindex[[30]]$I,   # Charge Transfer Capability
  aaindex[[63]]$I,   # Size
  aaindex[[111]]$I,  # Polarity
  aaindex[[177]]$I,  # Refractivity
  aaindex[[213]]$I,  # Average non-bonded energy per atom
  aaindex[[319]]$I   # Accessible surface area
))

colnames(frame) <- c('P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8')
frame.names <- names(aaindex[[2]]$I)
frame <- scale(frame)

# write.csv(frame, args[[1]])
write.table(frame, args[[1]], append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

cat('Successfully wrote data to:', args[[1]], '\n')
cat('Column means:\n')
print(colMeans(frame))