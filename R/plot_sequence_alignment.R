#' Plot sequence alignment
#'
#' This function generates a sequence alignment plot using ggplot2 based on the input alignment table.
#'
#' @param alignment_tbl An alignment table containing query and subject information for sequence alignment.
#' @return A ggplot object representing the sequence alignment plot.
#'
#' @examples
#' q <- (c("boo", "fizzbuzz"))
#' s <- ("boofizz")
#'
#' alignment <- alignment_table(q, s)
#' plot_sequence_alignment(alignment_tbl = alignment)
#'
#' # Provide names for query and subject to label the y-axis
#' names(query) <- c("Seq1", "Seq2")
#' names(subject) <- "reference"
#'
#' alignment <- alignment_table(query, subject)
#' plot_sequence_alignment(alignment_tbl = alignment)
#'
#'
#'
#' @import ggplot2 dplyr
#'
#' @export
plot_sequence_alignment <- function(alignment_tbl = alignment_table(query, subject)){
  out_left <- ggplot(alignment_tbl[["mismatch"]], aes(y = PatternName)) +
    geom_linerange(aes(xmin = start, xmax = end), stat = "unique") +
    geom_linerange(data = alignment_tbl[["deletions"]],
                 aes(xmin = start-0.5, xmax = end+0.5),
                 size = 2, color = "white", stat = "unique") +
    geom_linerange(data = alignment_tbl[["insertions"]],
                 aes(xmin = start-width/2, xmax = end-width/2, y = as.numeric(PatternName)+0.4),
                 color = "purple", stat = "unique") +
    geom_linerange(data = alignment_tbl[["insertions"]],
                   aes(x = start-0.6, ymin = PatternName,
                       ymax = as.numeric(PatternName)+0.3),
                   linetype = "dashed" , stat = "unique") +
    geom_linerange(data = alignment_tbl[["insertions"]],
                   aes(x = start-0.4, ymin = PatternName,
                       ymax = as.numeric(PatternName)+0.3),
                   linetype = "dashed", stat = "unique") +
    geom_segment(data = alignment_tbl[["insertions"]],
                   aes(x = start-0.6, xend = start-width/2,
                       y = as.numeric(PatternName)+0.3,
                       yend = as.numeric(PatternName)+0.4),
                 linetype = "dashed", stat = "unique") +
    geom_segment(data = alignment_tbl[["insertions"]],
                 aes(x = start-0.4, xend = end-width/2,
                     y = as.numeric(PatternName)+0.3,
                     yend = as.numeric(PatternName)+0.4),
                 linetype = "dashed", stat = "unique") +
    theme_classic() +
    theme(axis.title = element_blank(),
          legend.position = "none") +
    ggtitle(paste(unique(alignment_tbl$clus_name_p)))

  if(any(!is.na(alignment_tbl[["mismatch"]]$PatternSubstring))){
    out_left <- out_left + geom_point(aes(x = SubjectStart,
                                          color = PatternSubstring)) +
      geom_text(aes(x = SubjectStart, label = PatternSubstring,
                    color = PatternSubstring),
                nudge_y = 0.2, hjust = 0.5)
  }
  return(out_left)
}
