

get_mutation_hotspots = function(pipeline = "oncodriveclustl"){
  
  #matched only
  oncoclustl = "all_the_things/oncodriveclustl-1.1/99-outputs/grch37/BL_DLBCL-BL-like_DLBCL_FL_all--2025-01/3041a58d7bcb6208c51e119149b5dd52/cds/"
  
}



library(ggseqlogo)
lg = get_gambl_colours("lymphgen")

cs1 = make_col_scheme(chars=c('G', 'A', 'M', 'B','L','R','O','P','E','N'),
                      groups=c('G', 'A', 'M', 'B','L','R','O','P','E','N'), 
                      cols=c(lg["EZB"], lg["ST2"], lg["MCD"], lg["BN2"],lg["N1"],
                      lg["A53"],'grey','grey','grey','grey'))

ggplot(  ) + geom_logo(c("GAMBLR","GAMBLR","GOPENR","GAMBLR",
                         "GOPENR","GAMBLR","GAMBNR","GAMBLR",
                         "GOPENR","GOPELR","GAMENR"),seq_type = "other",
                       namespace="ABCDEFGHIJKLMNOPQRSTUVWXYZ",method="probability",
                       col_scheme=cs1) + 
  theme_void()



ggplot(  ) + geom_logo(c("GAMBLR","GAMBLR","GAMBLR","GAMBLR","GOPENR","GAMBLR",
                         "GOPENR","GAMBLR","GAMBNR","GAMBLR",
                         "GOPENR","GOPELR","GAMBLR"),seq_type = "other",
                       namespace="ABCDEFGHIJKLMNOPQRSTUVWXYZ",method="probability",
                       col_scheme=cs1) + 
  theme_void()

