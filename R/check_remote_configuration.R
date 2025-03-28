#' @title Check Remote Configuration.
#'
#' @description Check for a remote session and automagically confirm setup will work properly.
#'
#' @details This function determines if a user is working in GAMBLR remotely and, if so, will
#' check if their config is loaded properly and ssh_session is available.
#'
#' @param auto_connect Set to TRUE to ensure an ssh_session is created if absent
#'
#' @return TRUE if a remote session is detected, FALSE otherwise.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' check_remote_configuration()
#' }
#' @keywords internal
check_remote_configuration = function(auto_connect = FALSE){

  remote_gamblr = check_host(auto_connect=auto_connect)
  if(remote_gamblr){
    #code is running on a non-GSC computer. Check that the config is set up properly
    if(!Sys.getenv("R_CONFIG_ACTIVE")=="remote"){
      stop('You seem to be running this on a remote computer but have not set R_CONFIG_ACTIVE properly\nTry running Sys.setenv(R_CONFIG_ACTIVE= "remote")')
    }
  }
  return(remote_gamblr)
}
