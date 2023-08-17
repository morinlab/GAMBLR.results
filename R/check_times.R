check_times = function(relative_paths,
                       archive_mode = FALSE,
                       force_backup = FALSE){

  local_base = check_config_value(config::get("project_base"))
  remote_base = check_config_value(config::get("project_base",config="default"))
  if(archive_mode){
    if(local_base == remote_base){
      message("checking against local archive")
    }else{
      message("Currently, this mode must be run on the GSC (not remotely)")
      return(NULL)
    }
    local_base = check_config_value(config::get("archive"))

  }
  for(rel_f in relative_paths){
    local_f = paste0(local_base,rel_f)
    remote_f = paste0(remote_base,rel_f)
    print(rel_f)
    if(file.exists(local_f)){
      mtime = file.info(local_f)$mtime
      mtime = stringr::str_remove(mtime,"\\s\\d+:\\d+:\\d+")
      #print(mtime)

      remote_session = check_remote_configuration(auto_connect=TRUE)
      #print(remote_f)
      if(remote_session){
        output = ssh::ssh_exec_internal(ssh_session,paste("stat -L ",remote_f,"| grep Modify"))$stdout

        output = rawToChar(output) %>% stringr::str_extract(.,"\\d+-\\d+-\\d+")
      }else{
        output = file.info(remote_f)$mtime %>% stringr::str_remove("\\s\\d+:\\d+:\\d+")
      }
      remote_time = lubridate::as_date(output)
      local_time = lubridate::as_date(mtime)
      agediff = lubridate::time_length(remote_time - local_time,unit="days")
      if(agediff>0){
        print(paste("Warning! Remote version is",agediff,"days newer than the local file. You probably need to update the following file:"))
        print(rel_f)
      }else{
        message("OK")
      }
    }else{
      if(archive_mode){
        message("local backup of this file doesn't seem to exist:")
        message(local_f)

        copy_no_clobber(remote_f,local_f,force_backup)

      }
    }
  }
}
