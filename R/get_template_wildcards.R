get_template_wildcards = function(parent_key,
                                  template_key){

  if(missing(template_key)){
    wildcard_string = config::get(parent_key)
  }else{
    wildcard_string = config::get(paste0(parent_key,"_wildcards"))[template_key]
  }
  wildcards = stringr::str_split(wildcard_string,",")
  return(unlist(wildcards))
}
