# A row filter for SV tables: rows of the given samples (if any), and with
# either breakpoint in the region (if given); chromosome names are compared
# without a chr prefix
.sv_row_filter = function(sample_ids = NULL, chromosome = NULL, qstart = NULL, qend = NULL){
  chrom = if(is.null(chromosome)) NULL else gsub("^chr", "", chromosome)
  function(x){
    if(!is.null(sample_ids)) x = x[x$tumour_sample_id %in% sample_ids, , drop = FALSE]
    if(!is.null(chrom)){
      a = gsub("^chr", "", x$CHROM_A) == chrom & x$START_A >= qstart & x$START_A <= qend
      b = gsub("^chr", "", x$CHROM_B) == chrom & x$START_B >= qstart & x$START_B <= qend
      x = x[!is.na(a | b) & (a | b), , drop = FALSE]
    }
    x
  }
}
