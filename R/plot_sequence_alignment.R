#' Plot sequence alignment
#'
#' This function generates a sequence alignment plot using ggplot2 based on the
#' input alignment table.
#'
#' @param alignment_tbl An alignment table containing query and subject
#' information for sequence alignment. Generated with [alignment_table()].
#' @param insertion_color The color to use for indicating insertions in the
#' alignment. Default is '#21918c'. Can be any output of `colors()` or hex code.
#' @param hide_mismatches A logical value indicating whether to hide mismatches
#' in the alignment plot. Default is FALSE.
#' @return A ggplot object of the sequence alignment plot.
#'
#' @examples
#' q <- (c("boo", "fibububuzz", "bozz", "baofuzz"))
#' s <- "boofizz"
#'
#' alignment <- alignment_table(q, s)
#' pl1 <- plot_sequence_alignment(alignment_tbl = alignment)
#' pl1
#'
#' # Provide names for (some) query and subject elements to label the y-axis
#' names(q) <- c("Seq1", NA, "Seq3")
#' names(s) <- "reference"
#' pl2 <- plot_sequence_alignment(alignment_table(q, s))
#' pl2
#'
#' # Compatible with StringSets from Biostrings
#' dna <- readDNAStringSet(system.file("extdata", "dm3_upstream2000.fa.gz",
#'                                     package="Biostrings"))
#' # The entries dna[2:5] are identical
#' q <- dna[2:4]
#' s <- dna[5]
#' pl3 <- plot_sequence_alignment(alignment_table(q, s))
#' pl3
#'
#' # Let's introduce some SNPs, insertions and deletions
#' q[1] <- as(Biostrings::replaceLetterAt(q[[1]], c(5, 200, 400), "AGC"),
#'             "DNAStringSet")
#' q[2] <- as(c(substr(q[[2]], 300, 1500), substr(q[[2]], 1800, 2000)),
#'            "DNAStringSet")
#' q[3] <- as(Biostrings::replaceAt(q[[3]], 1500,
#'              paste(rep("A", 1000), collapse = "")),
#'            "DNAStringSet")
#' names(q) <- c("mismatches", "deletions", "insertion")
#' names(s) <- substr(names(s)[1], 1, 34)
#'
#' pl4 <- plot_sequence_alignment(alignment_table(q, s))
#' pl4
#'
#' # Compatible with ggplot2 theming
#' pl4 +
#'   ylab("Sequence variants") +
#'   xlab("Length in bp") +
#'   scale_color_viridis_d() +
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#'         axis.title = element_text())
#'
#'
#' @import ggplot2 dplyr
#'
#' @export
plot_sequence_alignment <- function(alignment_tbl =
                                      alignment_table(query, subject),
                                    insertion_color = "#21918c",
                                    hide_mismatches = FALSE){
  # line themeing
  lt <- "dashed"
  co <- "grey"
  ins_w <- 2
  u <- "unique"

  # plot lines representing strings
  pl <- ggplot(filter(alignment_tbl, feature == "mismatch"),
               aes(y = PatternName)) +
    geom_linerange(aes(xmin = start, xmax = end), stat = "unique") +
    # plot white lines representing deletions over the original lines,
    # causing the illusion of a gap.
    geom_linerange(data = filter(alignment_tbl, feature == "deletion"),
                 aes(xmin = start, xmax = end),
                 linewidth = ins_w, color = "white", stat = u) +
    # plot insertions centered on the point of insertion with dashed lines
    # indicating the point of insertion
    # This is the actual insertion
    geom_linerange(data = filter(alignment_tbl, feature == "insertion"),
                 aes(xmin = start-width/2, xmax = end-width/2,
                     y = as.numeric(PatternName)+0.4),
                 color = insertion_color, stat = u) +
    # This is the left leg of the indicator of insertion point
    geom_linerange(data = filter(alignment_tbl, feature == "insertion"),
                   aes(x = start-0.1, ymin = PatternName,
                       ymax = as.numeric(PatternName)+0.3),
                   linetype = lt , stat = u, color = co) +
    # This is the right leg of the indicator of insertion point
    geom_linerange(data = filter(alignment_tbl, feature == "insertion"),
                   aes(x = start+0.1, ymin = PatternName,
                       ymax = as.numeric(PatternName)+0.3),
                   linetype = lt , stat = u, color = co) +
    # This is the left arm of the indicator of insertion point
    geom_segment(data = filter(alignment_tbl, feature == "insertion"),
                 aes(x = start-0.1, xend = start-width/2,
                       y = as.numeric(PatternName)+0.3,
                       yend = as.numeric(PatternName)+0.4),
                 linetype = lt , stat = u, color = co) +
    # This is the right arm of the indicator of insertion poin
    geom_segment(data = filter(alignment_tbl, feature == "insertion"),
                 aes(x = start+0.1, xend = end-width/2,
                     y = as.numeric(PatternName)+0.3,
                     yend = as.numeric(PatternName)+0.4),
                 linetype = lt , stat = u, color = co) +
    theme_classic() +
    theme(axis.title = element_blank(),
          legend.position = "none")

  # Add mismatches as point plot if any exist
  if(any(!is.na(alignment_tbl$PatternSubstring)) &
     !hide_mismatches){
    pl <- pl + geom_point(aes(x = SubjectStart,
                                          color = PatternSubstring)) +
      geom_text(aes(x = SubjectStart, label = PatternSubstring,
                    color = PatternSubstring),
                nudge_y = 0.2, hjust = 0.5)
  }

  return(pl)
}
