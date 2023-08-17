load_combine_cached_bins = function(function_params=list(bin_size=15000),
  additional_details=list(foreground="DLBCL_FL",background="CLL_MM")){

      cache_file_name = paste0(check_config_value(config::get("repo_base")),"cached_results/get_ssm_by_region")
      for (param in names(function_params)[order(names(function_params))]){
        cache_file_name = paste0(cache_file_name,"--",param,"-",function_params[[param]])
      }
      cache_file_name = paste0(cache_file_name,"*")
      for (detail in names(additional_details)){
        cache_file_name = paste0(cache_file_name,"--",detail,"-",additional_details[[detail]])
      }
      cache_file_name = paste0(cache_file_name,".tsv")
     
    cached = Sys.glob(paths = cache_file_name )
  if(!length(cached)==48){
    stop(paste("error: cannot find 48 files using this pattern",cache_file_name))
  }
  binned_list = list()
  for(to_add in cached){
    this_df = suppressMessages(read_tsv(to_add,col_types="ciicnnnn"))
    binned_list[[to_add]]=this_df
  }
  return(bind_rows(binned_list))
}
