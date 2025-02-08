#' @title Check Host.
#'
#' @description Check if code is running remotely and react accordingly.
#'
#' @details The function will (optionally) attempt a connection if necessary, and stores it in a global variable (ssh_session).
#'
#' @param auto_connect Set to TRUE if you want the function to create an ssh session (if necessary).
#' @param verbose Set this to TRUE for verbose messages from the function.
#'
#' @return TRUE if a remote session is detected, FALSE otherwise.
#'
#' @export
#'
#' @examples
#' check_host(auto_connect=TRUE)
#'
#' @keywords internal
check_host = function(auto_connect = FALSE,
                      verbose = FALSE){

  hostname = Sys.info()["nodename"]
  if(grepl("bcgsc.ca",hostname)){
    #we are on the GSC network
    return(FALSE)
  }else{
    # we are on some computer not on the GSC network (needs ssh_session)
    if(exists("ssh_session") && class(ssh_session)=="ssh_session"){
      if(verbose){
        message("active ssh session detected")
      }
    }else{
      if(auto_connect){
        session = get_ssh_session()
        assign("ssh_session", session, envir = .GlobalEnv)
      }else{
        if(verbose){
          message("You appear to be using GAMBLR on your local computer. Be sure to set up an ssh session!")
          message("requires an active VPN connection to the GSC")
          message("?GAMBLR.results::get_ssh_session for more info")
        }
      }
    }
  }
  return(TRUE)
}
