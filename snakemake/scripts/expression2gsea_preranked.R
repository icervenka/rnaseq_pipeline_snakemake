suppressMessages(library(lazyeval))
suppressMessages(library(dplyr))

create_gsea_rank <- function(data,
                             out = NULL,
                             ranking_formula = ~ -log10(pvalue) * sign(log2FoldChange), # nolint
                             id_colname = ensembl_gene_id,
                             .inf = "replace") {
  out_df <- data %>%
    dplyr::mutate(rank = !!lazyeval::f_rhs(ranking_formula)) %>%
    dplyr::select({{ id_colname }}, rank) %>%
    dplyr::arrange(dplyr::desc(rank)) %>%
    tidyr::drop_na()

  # process the +/- Inf values
  if (.inf == "drop") {
    out_df <- out_df %>%
      dplyr::filter(is.finite(rank))
  } else if (.inf == "replace") {
    # find min and max values that are not infinite
    min_val <- out_df %>%
      dplyr::filter(is.finite(rank)) %>%
      dplyr::pull(rank) %>%
      min()
    max_val <- out_df %>%
      dplyr::filter(is.finite(rank)) %>%
      dplyr::pull(rank) %>%
      max()

    # add 1% to the max value and replace the +/- Inf with new value
    out_df$rank[which(out_df$rank == Inf)] <- max_val + 0.01 * max_val
    out_df$rank[which(out_df$rank == -Inf)] <- min_val + 0.01 * min_val
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
