deblacklist_regions = function(regions_bed,projection="grch37"){
  if(projection=="grch37"){
    blacklist_file = "/Users/rmorin/git/LLMPP/resources/reference/encode/hg19-blacklist.v2.bed"
  }
  blacklist_df = suppressMessages(read_tsv(blacklist_file,col_names = c("chr","start","end","name")))
  blacklist_dt = as.data.table(blacklist_df)
  colnames(regions_bed)[c(1:3)]=c("chr","start","end")
  if(any(!grepl("chr",regions_bed[,1]))){
    message("adding chr prefix")
    regions_bed = mutate(regions_bed,chr = as.character(chr))
    regions_bed = mutate(regions_bed,chr = paste0("chr",chr))
  }
  regions_dt = as.data.table(regions_bed)

 
  #annotate bins that overlap blacklisted regions
  setkey(blacklist_dt, chr,start,end)
  setkey(regions_dt, chr,start,end)

  annotated = foverlaps(blacklist_dt, regions_dt,by.x =c("chr","start","end"),by.y=c("chr","start","end"),mult = "all",type = "any") %>% dplyr::filter(!is.na(start))
  to_keep = dplyr::filter(regions_bed,!start %in% annotated$start)
  return(to_keep)
}
