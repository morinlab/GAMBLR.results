#' @title Read Merge Manta With Liftover.
#'
#' @description Takes a path to bedpe and runs liftover ([GAMBLR.results::liftover_bedpe]) based on the original genome build of the bedpe.
#'
#' @details This is a helper function that is not meant to be used routinely.
#'
#' @param bedpe_paths path to bedpe
#' @param pattern pattern
#' @param out_dir output directory
#'
#' @import dplyr readr
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' manta_bedpe = read_merge_manta_with_liftover(bedpe_paths = "some_path.bedpe",
#'                                              out_dir = "../")
#' }
#'
read_merge_manta_with_liftover = function(bedpe_paths = c(),
                                          pattern = "--matched",
                                          out_dir){

  to_merge = list()
  print(head(bedpe_paths))
  for(thispath in bedpe_paths){
    sample_paths = dir(thispath, pattern = pattern) #DEBUGGING
    print(sample_paths)
    #sample_paths = head(sample_paths,15) #for debugging
    for(sp in sample_paths){
      full_path = paste0(thispath, sp, "/somaticSV.bedpe")
      print(paste("working on HERE:", full_path))
      if(grepl("hg38", full_path)){
        print("using liftOver")
        svbed = liftover_bedpe(full_path) #load and convert to grch37 coordinates
      }else{
        svbed = suppressMessages(read_tsv(full_path, comment = "##", col_types = "cddcddccccccccccccccccc"))
      }

      this_patient = colnames(svbed)[23]
      this_normal = colnames(svbed)[22]
      out_file = paste0(out_dir, "/", this_patient, "--", this_normal, "--hg38Togrch37_sv.tsv")
      print(paste("writing output to", out_file))
      infos = pull(svbed, this_patient)
      infos_n = pull(svbed, this_normal)
      colnames(svbed)[c(1:6)] = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B")
      #all_vafs = get.sv.vaf(infos)
      #svbed$VAF = as.numeric(all_vafs)
      svbed$VAF_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
      svbed$DP_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")),2)[1])})
      svbed$VAF_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
      svbed$DP_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 2)[1])})
      svbed$SOMATIC_SCORE = sapply(svbed$INFO_A, function(x){as.numeric(tail(unlist(strsplit(x, "=")), 1))})

      #filter on PASS, score, VAF
      #svbed_filt = svbed %>%
      #  filter(SCORE > minScore & FILTER == "PASS") %>%
      #  dplyr::select(c(chrom1, start1, end1, chrom2, start2, end2))

      svbed$tumour_sample_id = this_patient
      svbed$normal_sample_id = this_normal
      if(grepl("--unmatched", sp)){
        svbed$pair_status = "unmatched"
      }else{
        svbed$pair_status = "matched"
      }
      print(head(svbed))
      svbed$NAME = "."

      svbed = svbed %>%
        dplyr::select(CHROM_A, START_A, END_A, CHROM_B, START_B, END_B, NAME, SOMATIC_SCORE, STRAND_A, STRAND_B, TYPE, FILTER, VAF_tumour, VAF_normal, DP_tumour, DP_normal, tumour_sample_id, normal_sample_id, pair_status)
      #remove chr prefix from both chromosome names
      svbed = svbed %>%
        mutate(CHROM_A = gsub("chr", "", CHROM_A)) %>%
        mutate(CHROM_B = gsub("chr", "",CHROM_B))

      write_tsv(svbed, out_file, col_names = FALSE)
      #to_merge[[this_patient]] = svbed
    }
  }
}
