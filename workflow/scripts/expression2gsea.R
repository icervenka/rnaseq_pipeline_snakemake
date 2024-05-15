suppressMessages(library(stringr))
suppressMessages(library(dplyr))

create_gsea_normalized <- function(data,
                                   metadata,
                                   id_colname,
                                   data_description_col = NULL,
                                   sample_colname = sample,
                                   out = NULL) {
  # create string from supplied variables, needed for certain functions
  samples <- metadata[[deparse(substitute((sample_colname)))]]
  data_cols <- names(data)[names(data) %in% samples]

  # create the data frame to output
  # see https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
  out_df <- data %>%
    dplyr::select(NAME = {{ id_colname }}, dplyr::any_of(data_cols))

  # add description column based on user settings
  if (is.null(data_description_col)) {
    out_df <- out_df %>%
      dplyr::mutate(DESCRIPTION = "") %>%
      dplyr::relocate(DESCRIPTION, .after = NAME)
  } else {
    out_df <- out_df %>%
      dplyr::mutate(DESCRIPTION = {{ data_description_col }}) %>%
      dplyr::relocate(DESCRIPTION, .after = NAME)
  }

  # write to file if specified
  if (!is.null(out)) {
    write.table(out_df,
      out,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  return(out_df)
}

create_gsea_cls <- function(data,
                            metadata,
                            out = NULL,
                            sample_colname = sample,
                            group_colname = group) {
  metadata_filtered <- metadata %>%
    dplyr::select({{ sample_colname }}, {{ group_colname }})

  # create string from supplied variables, needed for certain functions
  samples <- metadata_filtered[[deparse(substitute((sample_colname)))]]
  groups <- metadata_filtered[[deparse(substitute((group_colname)))]] %>%
    stringr::str_replace_all("[:space:]+", "")

  # which data columns and sample indices to include
  data_cols <- names(data)[names(data) %in% samples]
  sample_indices <- match(samples, data_cols)

  # compse the text of the final file
  # see https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
  text_lines <- c(
    paste(
      samples %>% unique() %>% length(),
      groups %>% unique() %>% length(),
      "1"
    ),
    paste("#", paste(groups %>% unique(), collapse = " ")),
    paste(groups[order(sample_indices)], collapse = " ")
  )

  if (!is.null(out)) {
    file_conn <- file(out)
    writeLines(text_lines, file_conn)
    close(file_conn)
  }
  return(text_lines)
}