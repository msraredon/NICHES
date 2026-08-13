#' Load OmniPath as ground truth
#'
#' @param species Species reference in the Ominipath database
#'
#' @export

LoadOmniPath <- function(species){

# Setup species call
  if (species == 'human'){
    organism = 9606
    }else if (species == 'mouse'){
      organism = 10090
      }else if(species == 'rat'){
        organism = 10116
        }else{
          stop("\nPlease select species for OmniPath mapping. Allows 'human','mouse',or 'rat' ")
        }

#Ligand-Receptor Network ----
# NOTE (NICHES issue #70): newer OmnipathR (>= ~3.17, including 4.1.0) resolves the
# `organism` argument through OmnipathR::ncbi_taxid(), which builds an organism-name
# table by downloading and coalescing species lists from Ensembl, OMA and UniProt.
# When one of those ancillary services is unreachable (e.g. omabrowser.org returning
# HTTP 502), that step aborts with "Can't combine `..1` and `..3`" - even though the
# OmniPath interaction server itself is fine. Since NICHES has already mapped the
# species to its NCBI taxon id above, we do not actually need that translation. The
# helper below first tries the normal OmnipathR path (preserving its caching when the
# services are healthy) and, if that fails, falls back to querying the OmniPath REST
# API directly with the known taxon id, which bypasses the organism-name lookup.
lr_raw <- .niches_fetch_omnipath_ligrec(organism)
omnipath_source <- attr(lr_raw, "niches_omnipath_source") # provenance, retained on the output
lr_Interactions_Omnipath <- lr_raw %>%
  dplyr::select(source_genesymbol,target_genesymbol) %>%
  dplyr::distinct()

#Tag with mechanism name
lr_Interactions_Omnipath$mechanism <- paste(lr_Interactions_Omnipath$source_genesymbol,lr_Interactions_Omnipath$target_genesymbol,sep = '-')

# Identify max number of ligand subunits and max number of receptor subunits (based on "_" as a separator, used in current OmniPath iteration as of 2021-06-07)
source_sub_max <- max(sapply(lr_Interactions_Omnipath$source_genesymbol,function(x) length(strsplit(x,split="_")[[1]])))
target_sub_max <- max(sapply(lr_Interactions_Omnipath$target_genesymbol,function(x) length(strsplit(x,split="_")[[1]])))

# Initialize column names based on how many subunits are in initial database
source_col_names <- paste0("source_",c(1:source_sub_max))
target_col_names <- paste0("target_",c(1:target_sub_max))

# Split into individual columns
temp <- tidyr::separate(data = lr_Interactions_Omnipath,
                        col = source_genesymbol, # Split Source genes
                        into = source_col_names, # Uses initialized column names
                        sep = '_',
                        remove = F)
temp <- tidyr::separate(data = temp,
                        col = target_genesymbol, # Split Target genes
                        into = target_col_names, # Uses initialized column names
                        sep = '_',
                        remove = F)

# Export subunit dataframe
source.subunits <- as.matrix(temp[,source_col_names]) #allows duplicate rownames
rownames(source.subunits) <- temp$source_genesymbol
target.subunits <- as.matrix(temp[,target_col_names]) #allows duplicate rownames
rownames(target.subunits) <- temp$target_genesymbol

ground.truth <- list('source.subunits' = source.subunits,
                     'target.subunits' = target.subunits)
attr(ground.truth, 'omnipath_source') <- omnipath_source # keep data provenance on the object
return(ground.truth)
}

#' Resiliently fetch OmniPath ligand-receptor interactions
#'
#' Internal helper for \code{\link{LoadOmniPath}} that works around a known
#' upstream OmnipathR failure (NICHES issue #70). Newer OmnipathR versions resolve
#' the \code{organism} argument via an organism-name table built from Ensembl, OMA
#' and UniProt downloads; when one of those services is down the call aborts before
#' any interactions are retrieved. This helper first attempts the standard
#' OmnipathR path, and on failure falls back to a direct OmniPath REST query using
#' the already-known NCBI taxon id, which does not require the organism-name lookup.
#'
#' @param organism Integer NCBI taxonomy id (9606 human, 10090 mouse, 10116 rat).
#'
#' @return A data.frame of ligand-receptor interactions containing at least the
#'   \code{source_genesymbol} and \code{target_genesymbol} columns.
#'
#' @keywords internal
#' @noRd
.niches_fetch_omnipath_ligrec <- function(organism){

  # 1) Preferred path: standard OmnipathR call (uses its cache when services are up).
  #    The organism-name download noise is muted; if it fails we recover silently.
  result <- tryCatch(
    .niches_quiet_omnipath(OmnipathR::import_ligrecextra_interactions(organism = organism)),
    error = function(e) e
  )
  if (!inherits(result, "error")){
    ver <- tryCatch(as.character(utils::packageVersion("OmnipathR")), error = function(e) NA_character_)
    message(sprintf(
      "Loaded OmniPath ligand-receptor database via OmnipathR %s (retrieved %s): %d interactions.",
      ver, Sys.Date(), nrow(result)))
    attr(result, "niches_omnipath_source") <-
      list(via = "OmnipathR", version = ver, retrieved = as.character(Sys.Date()), n = nrow(result))
    return(result)
  }

  # 2) Fallback: query the OmniPath REST API directly with the known taxon id.
  #    This avoids OmnipathR's organism-name translation (the step that fails).
  fallback <- tryCatch(.niches_omnipath_rest_ligrec(organism), error = function(e) e)
  if (!inherits(fallback, "error") && nrow(fallback) > 0){
    ver <- .niches_omnipath_server_version()
    message(sprintf(
      "Loaded OmniPath ligand-receptor database (OmniPath server %s, retrieved %s): %d interactions.",
      if (is.na(ver)) "REST API" else ver, Sys.Date(), nrow(fallback)))
    attr(fallback, "niches_omnipath_source") <-
      list(via = "OmniPath REST", version = ver, retrieved = as.character(Sys.Date()), n = nrow(fallback))
    return(fallback)
  }

  # 3) Both paths failed: raise an informative error (this stays loud on purpose).
  stop(
    "Unable to load the OmniPath ligand-receptor database.\n",
    "  OmnipathR error : ", conditionMessage(result), "\n",
    "  REST fallback   : ",
    if (inherits(fallback, "error")) conditionMessage(fallback) else "returned no interactions", "\n",
    "The OmniPath interaction server (omnipathdb.org) may be unreachable from this\n",
    "machine (check network/proxy access to https://omnipathdb.org). See NICHES\n",
    "issue #70 for background.",
    call. = FALSE
  )
}

#' Run an OmnipathR expression with its console logging muted
#'
#' OmnipathR reports failed ancillary downloads (e.g. the OMA/Ensembl species
#' lists) through its own \code{logger}-based console appender, which is not
#' captured by \code{suppressWarnings}/\code{suppressMessages}. Since NICHES
#' handles those failures itself (issue #70), this helper temporarily raises the
#' OmnipathR console log threshold so the expected noise is hidden, restoring the
#' previous level afterwards. Genuine errors still propagate as R conditions.
#'
#' @param expr An expression calling OmnipathR (evaluated lazily inside here).
#'
#' @keywords internal
#' @noRd
.niches_quiet_omnipath <- function(expr){
  old <- getOption("omnipathr.console_loglevel", default = "success")
  try(OmnipathR::omnipath_set_console_loglevel("fatal"), silent = TRUE)
  on.exit(try(OmnipathR::omnipath_set_console_loglevel(old), silent = TRUE), add = TRUE)
  suppressWarnings(suppressMessages(expr))
}

#' Best-effort OmniPath web server version
#'
#' Reads the OmniPath service banner (\code{https://omnipathdb.org/about}) and
#' extracts the server version for provenance reporting. Returns \code{NA} if the
#' banner cannot be read or parsed - version reporting must never block a load.
#'
#' @return Character server version (e.g. "0.1.5"), or \code{NA_character_}.
#'
#' @keywords internal
#' @noRd
.niches_omnipath_server_version <- function(){
  tryCatch({
    txt <- paste(readLines("https://omnipathdb.org/about", warn = FALSE, n = 3), collapse = " ")
    m <- regmatches(txt, regexec("omnipath-server[[:space:]]+([0-9]+(\\.[0-9]+)*)", txt))[[1]]
    if (length(m) >= 2 && nzchar(m[2])) m[2] else NA_character_
  }, error = function(e) NA_character_)
}

#' Direct OmniPath REST query for ligrecextra interactions
#'
#' Downloads the \code{ligrecextra} ligand-receptor dataset straight from the
#' OmniPath web service for a given NCBI taxon id, bypassing OmnipathR's
#' organism-name translation. The web service performs ortholog translation
#' server-side, so mouse (10090) and rat (10116) return organism-appropriate
#' gene symbols, matching what OmnipathR would otherwise produce.
#'
#' @param organism Integer NCBI taxonomy id.
#'
#' @return A data.frame with (at least) \code{source_genesymbol} and
#'   \code{target_genesymbol} columns.
#'
#' @keywords internal
#' @noRd
.niches_omnipath_rest_ligrec <- function(organism){

  url <- sprintf(
    "https://omnipathdb.org/interactions?datasets=ligrecextra&organisms=%s&genesymbols=yes",
    as.integer(organism)
  )

  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  utils::download.file(url, destfile = tmp, quiet = TRUE, mode = "wb")

  df <- utils::read.delim(
    tmp,
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c("source_genesymbol", "target_genesymbol") %in% colnames(df))){
    stop("OmniPath REST response did not contain the expected gene symbol columns.")
  }
  # Drop rows without gene symbols on either side (mirrors usable interactions).
  df <- df[nzchar(df$source_genesymbol) & nzchar(df$target_genesymbol) &
             !is.na(df$source_genesymbol) & !is.na(df$target_genesymbol), , drop = FALSE]

  return(df)
}
