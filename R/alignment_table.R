#' Alignment Table
#'
#' Generate a table of mismatches and indels between one or many query sequences and a subject sequence.
#'
#' @param query An object of class XStringSet containing the query sequences.
#' @param subject An object of class XStringSet containing the subject sequence.
#'
#' @return A list containing tibbles with information on mismatches and indels.
#'
#' @examples
#' query_seq <- DNAStringSet(c("ACCGTACCTGG", "ACCTTGG"))
#' subject_seq <- DNAStringSet("ACCGTACCGGG")
#' x <- alignment_table(query_seq, subject_seq)
#' x
#'
#' @import dplyr
#' @importFrom Biostrings DNAStringSet AAStringSet pairwiseAlignment mismatchTable width
#'
#' @export
alignment_table <- function(query = XStringSet,
                            subject = XStringSet){
  if(is.null(names(query))){
    names(query) <- paste0("query", 1:length(query))
  }
  if(is.null(names(subject))){
    names(subject) <- "subject"
  }

  all_sequences <- c(query, subject)

  alig <- pairwiseAlignment(all_sequences, subject)

  ranges <- alig@subject@range %>% as.data.frame()

  deletions <- alig@pattern@indel %>% as.data.frame() %>%
    as_tibble() %>%
    left_join(tibble(group = 1:length(names(all_sequences)),
              PatternName = names(all_sequences)), by = "group") %>%
    mutate(PatternName = factor(PatternName,
                                levels = c(names(subject),
                                           names(query))))

  insertions <- alig@subject@indel %>% as.data.frame() %>%
    as_tibble() %>%
    left_join(tibble(group = 1:length(names(all_sequences)),
                     PatternName = names(all_sequences)), by = "group") %>%
    mutate(PatternName = factor(PatternName,
                                levels = c(names(subject),
                                           names(query))))

  mismatches <- alig %>%
    mismatchTable() %>%
    as_tibble()
  if(nrow(mismatches) == 0){
    mismatches <- add_row(mismatches, PatternId = 1:length(names(all_sequences)))
  }
  mismatches <- full_join(mismatches,
            tibble(PatternId = 1:length(names(all_sequences)),
                     PatternName = names(all_sequences),
                     PatternLength = width(all_sequences),
                     ranges), by = join_by(PatternId))

  mismatches_ref <- select(mismatches, contains("Subject")) %>%
    filter(!is.na(SubjectSubstring)) %>%
    dplyr::rename(PatternSubstring = "SubjectSubstring") %>%
    bind_cols(
      select(mismatches, PatternId, PatternName,
             PatternLength, start, end, width) %>%
        filter(PatternName == names(subject))
    )
  mismatches <- bind_rows(mismatches, mismatches_ref) %>%
    mutate(PatternName = factor(PatternName,
                                levels = c(names(subject),
                                           names(query)))) %>%
    mutate(end = ifelse(end < (PatternLength-sum(insertions$width)+start),
                        (PatternLength-sum(insertions$width)+start), end)) %>%
    unique()
  out <- list(mismatches, insertions, deletions)
  names(out) <- c("mismatch", "insertions", "deletions")
  return(out)
}
