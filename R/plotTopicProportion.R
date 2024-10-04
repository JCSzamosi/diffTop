#' Plot the topic distribution in all specimens.
#'
#' @param theta_aligned List. returned object of thetaAligned().
#' @param col_names_theta_all Character vector. Column names for long formatted theta. An example is c("iteration", "Sample", "Topic", "topic.dis")
#' @param SampleIDVarInphyloseq Character. SampleID variable name in phyloseq object.
#' @param design A factor of interest, in the form of a formula.
#' @inheritParams alignmentMatrix
#'
#' @return A ggplot2 object. Histrogram of median of topic proportion in each sample.
#' @import phyloseq
#' @import ggplot2
#' @importFrom dplyr mutate summarize group_by left_join ungroup
#' @importFrom stringr str_c
#' @importFrom stats median
#' @export
#'

plotTopicProportion <- function(
  ps,
  theta_aligned,
  K,
  col_names_theta_all = c("iteration", "Sample", "Topic", "topic.dis"),
  chain = 4,
  iterUse = 1000,
  SampleIDVarInphyloseq = "unique_names",
  design = ~pna
){
  Sample <- Topic <- topic.dis <- NULL
  median.topic.dis <- median <- NULL

  # extract the terms from design
  tf <- terms.formula(design)
  term <- attr(tf, 'term.labels')
  
  dimnames(theta_aligned)[[2]] <- phyloseq::sample_names(ps)
  dimnames(theta_aligned)[[3]] <- c(paste0("Topic_", seq(1,K)))


  theta_all <- reshape2::melt(theta_aligned)

  colnames(theta_all) <- col_names_theta_all

  theta_all$Chain <- paste0(
    "Chain ",
    rep(seq(1, chain),
        each = iterUse
    )
  )

  sam <- phyloseq::sample_data(ps) %>%
    data.frame()

  theta_all$Sample <- theta_all$Sample %>%
    as.character()

  theta_all <- dplyr::left_join(
    theta_all,
    sam,
    by =c("Sample" = SampleIDVarInphyloseq)
  )

  theta_all$Chain <- factor(theta_all$Chain)
  theta_all$Topic <- factor(theta_all$Topic)
  theta_all$Sample <- factor(theta_all$Sample)
  theta_all[,term] <- factor(theta_all[,term])


  theta_summary <- theta_all %>%
    dplyr::group_by(
      Sample,
      Topic,
      !!sym(term) 
    ) %>%
    dplyr::summarize(
      median.topic.dis = median(topic.dis)
    ) %>%
    dplyr::ungroup() %>%
    mutate(
      Topic = factor(
        Topic,
        levels = rev(stringr::str_c("Topic_",1:K)
        )
      )
    )

  p <- ggplot2::ggplot(
    theta_summary,
    aes(x = !!sym(term),
        y = Topic,
        fill = !!sym(term))
  )
  p <- p+
    ggplot2::geom_tile(
      aes(alpha = median.topic.dis)
    ) +
    ggplot2::facet_grid(
      .~Sample,
      scale = "free"
    ) +
    ggplot2::xlab(term) +
    ggplot2::scale_alpha(
      name = "median topic distribution"
    ) +
    ggplot2::theme_minimal(
      base_size = 20
    ) +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5),
      strip.text.x = element_text(angle = 90)
    )
  p
}
