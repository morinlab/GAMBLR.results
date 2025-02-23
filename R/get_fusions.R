#' @title Get gene fusions
#'
#' @description Retrieve fusions identified by FusionCatcher
#'
#' @param projection Which genome build 
#' @param verbose How chatty you want this experience
#' @param keep_genes Vector of genes to restrict to
#' @param drop_genes Vector of genes to drop (i.e. bad/artifacts)
#' @param remove_adjacent_pairs Set to TRUE to drop genes that are already in close proximity
#' @param these_samples_metadata The metadata to link to your fusion table
#' @param join_with_metadata Set to TRUE to get fusions along with the sample metadata in one data frame
#' @param harmonize_immunoglobulin_partners Attempt to clean up IG names so they are all consistently one of IGH, IGK or IGL
#'
#' @return A data frame in a bedpe-like format 
#'
#' @import dplyr readr glue
#' @export
#'
#' @examples
#' 
#' all_fusions = get_gene_fusions()
#'
get_gene_fusions = function(projection = "grch37",
                       verbose=F,
                       keep_genes,
                       drop_genes,
                       remove_adjacent_pairs=FALSE,
                       these_samples_metadata,
                       join_with_metadata=FALSE,
                       harmonize_immunoglobulin_partners=FALSE){
    blacklist = c("EEF1A1","SMG1","CTSS","CTSB","IKZF3","CTSA","TIMM23B","VHL","SNX29","SERINC5","APOL4","TMPO","HMGB2","HMGB1","FRG1","CD44","LINC02328","STK4",
                  "TPM4","TCL6","ERVK9-11","ZNF600","ZNF433-AS1","YARS2","TNFRSF14-AS1","WDR12","YY1AP1","CABYR","TTLL3",
                  "CCDC88C","AC073320.1","ACTB","SWAP70","PTMA","LYZ","CD74","NIBAN3","MLXIPL","SLC7A5","CNR2","HSPA8","AC020656.1","SPN","AC243829.1")
    
    harmonize = function(bedpe){
      bedpe = mutate(bedpe,gene1=case_when(gene1=="IGH@" ~ "IGH",
                                           grepl("^IGHG",gene1) ~ "IGH",
                                           gene1=="IGK@" ~ "IGK",
                                           grepl("^IGK",gene1) ~ "IGK",
                                           grepl("^IGL",gene1) ~ "IGL",
                                           TRUE ~ gene1)) %>%
        mutate(gene2=case_when(gene2=="IGH@" ~ "IGH",
                                   grepl("^IGHG",gene2) ~ "IGH",
                                   gene2=="IGK@" ~ "IGK",
                                   grepl("^IGK",gene2) ~ "IGK",
                                   grepl("^IGL",gene2) ~ "IGL",
                                   TRUE ~ gene2))
      return(bedpe)
    }
    if(!missing(drop_genes)){
      blacklist = unique(c(blacklist,drop_genes))
    }
    output_base = GAMBLR.helpers::check_config_value(config::get("project_base"))
    
    genome_build= "hg38"
    if(projection == "grch37"){
      genome_build= "hg19"
    }
    fusion_path = paste0(output_base,"all_the_things/fusioncatcher-1.0/level_3/merges/all_fusions.{genome_build}.bedpe")
    fusion_path = glue::glue(fusion_path)
    if(verbose){
      print(fusion_path)
    }
    all_fusions = suppressMessages(read_tsv(fusion_path,
                                            progress = FALSE)) %>%
      separate(FUSION,into=c("gene1","gene2"),sep=":")
    if(!missing(keep_genes)){
      all_fusions = dplyr::filter(all_fusions,gene1 %in% keep_genes | gene2 %in% keep_genes)
    }
    all_fusions = dplyr::filter(all_fusions,!gene1 %in% blacklist) %>%
      dplyr::filter(!gene2 %in% blacklist)
    if(remove_adjacent_pairs){
      all_fusions = dplyr::filter(all_fusions,!grepl("10K",FLAGS)) 
    }
    if(!missing(these_samples_metadata)){
      if(join_with_metadata){
        all_fusions = left_join(all_fusions,these_samples_metadata)
      }
    }
    if(harmonize_immunoglobulin_partners){
      all_fusions = harmonize(all_fusions)
      
    }
    return(all_fusions)
}
