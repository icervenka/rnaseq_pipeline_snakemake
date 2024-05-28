suppressMessages(library(magrittr))

# global constants for ease of dir access
CD1UP <- "../" # nolint
CD2UP <- "../../" # nolint
CD3UP <- "../../../" # nolint


modify_tilda <- function(x, add = TRUE) {
  if (add == TRUE) {
    if (!startsWith(x, "~")) {
      x <- c("~ ", x)
    }
  } else {
    x <- gsub("~", "", x) %>% trimws()
  }
  return(x)
}

name_contrast <- function(c1, c2, sep = "_vs_") {
  return(paste0(c1, sep, c2))
}

get_species_info <- function(species_identifier) {
  species_info <- list(
    mouse_info = c(
      "full_name" = "Mus musculus",
      "abbreviation" = "mmusculus",
      "genus" = "Mus",
      "common" = "mouse",
      "kegg" = "mmu",
      "orgdb" = "org.Mm.eg.db",
      "orgdb_short" = "Mm",
      "meshdb" = "MeSH.Mmu.eg.db",
      "symbol" = "mgi",
      "tax" = 10090
    ),
    human_info = c(
      "full_name" = "Homo sapiens",
      "abbreviation" = "hsapiens",
      "genus" = "Homo",
      "common" = "human",
      "kegg" = "hsa",
      "orgdb" = "org.Hs.eg.db",
      "orgdb_short" = "Hs",
      "meshdb" = "MeSH.Hsa.eg.db",
      "symbol" = "hngc",
      "tax" = 9606
    ),
    rat_info = c(
      "full_name" = "Rattus Norvegicus",
      "abbreviation" = "rnorvegicus",
      "genus" = "Rattus",
      "common" = "rat",
      "kegg" = "rno",
      "orgdb" = "org.Rn.eg.db",
      "orgdb_short" = "Rn",
      "meshdb" = "MeSH.Rno.eg.db",
      "symbol" = "rds",
      "tax" = 10116
    )
  )
  for (x in species_info) {
    if (species_identifier %in% x) {
      sp_arr <- x
      break
    }
  }
  return(sp_arr)
}

translate_gene_ids <- function(gene_ids,
                               sp_info,
                               from_type,
                               to_type = c(
                                 "ENTREZID", "ENSEMBL",
                                 "SYMBOL", "GENENAME"
                               ),
                               method = "bitr",
                               drop = TRUE,
                               fill = TRUE) {
  if (method == "bitr") {
    ids <- suppressMessages(clusterProfiler::bitr(gene_ids,
      fromType = from_type,
      toType = to_type,
      OrgDb = sp_info["orgdb"],
      drop = drop
    ))
    if (!drop && fill) {
      ids <- purrr::map_dfc(ids, function(x) {
        dplyr::coalesce(x, ids[[from_type]])
      }) %>% tidyr::drop_na()
    }
    ids <- ids[!duplicated(ids[[from_type]]), ]
  } else if (method == "biomart") {

    # TODO default to_type
    suppressMessages(library(biomaRt))
    mart <- useMart(
      "ensembl",
      paste0(
        sp_info[["abbreviation"]],
        "_gene_ensembl"
      )
    )
    ids <- getBM(
      filters = "ensembl_gene_id",
      attributes = c(
        "entrezgene_id",
        "ensembl_gene_id",
        paste0(
          sp_info[["symbol"]],
          "_symbol"
        )
      ),
      values = gene_ids,
      mart = mart
    ) %>%
      dplyr::rename(c("ENTREZID", "ENSEMBL", "SYMBOL"))

    ids <- ids[!duplicated(ids$ENSEMBL), ]

    if (drop == TRUE) {
      ids <- ids %>% na.omit()
    }
  }
  return(ids)
}

add_idcolnames <- function(d,
                           sp_info,
                           id_column_name = "ENSEMBL",
                           id_type = "ENSEMBL",
                           required_colnames = c(
                             "ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"
                           )) {
  d <- d %>% dplyr::rename(!!id_type := dplyr::all_of(id_column_name)) # nolint
  data_colnames <- names(d)

  columns_in <- required_colnames %in% data_colnames

  if (!any(columns_in)) {
    stop(paste0(
      "No ID types present in the file. Has to be one of:\n",
      paste(required_colnames, collapse = " ")
    ))
  }

  if (!all(columns_in)) {
    id_column <- required_colnames[columns_in][1]
    ids <- translate_gene_ids(d[, id_column],
      sp_info,
      id_column,
      to_type = c(required_colnames)
    )
    # TODO remove duplicates
    merge_ids <- merge(ids, d, by = required_colnames[columns_in]) %>%
      dplyr::select(dplyr::all_of(required_colnames), dplyr::everything())
    return(merge_ids)
  } else {
    return(d %>% dplyr::select(required_colnames, dplyr::everything()))
  }
}

default_dt <- function(x, opts = NULL) {
  def_opts <- list(
    dom = "Blrtip", # specify content (search box, etc)
    deferRender = TRUE,
    scrollY = 600,
    scroller = TRUE,
    buttons = list(
      I("colvis"),
      "csv",
      "excel"
    )
  )

  DT::datatable(x,
    rownames = TRUE,
    escape = FALSE,
    filter = "top",
    extensions = c("Buttons", "Scroller"),
    style = "bootstrap",
    class = "compact stripe",
    width = "100%",
    options = c(def_opts, opts)
  )
}

get_order_func <- function(mode) {
  if (mode == "mean") {
    return(rowMeans)
  } else if (mode == "diff") {
    return(rowMeans)
  } else if (mode == "variance") {
    return(matrixStats::rowVars)
  } else {
    print("Order function for heatmap not recognized,
    ordering by row variance.")
    return(matrixStats::rowVars)
  }
}
