#' Alignment Table
#'
#' Generate a table of mismatches and indels between one or many query
#' sequences and a subject sequence.
#'
#' @param query A string or vector of strings or object of class XStringSet
#' containing the query sequences/strings.
#' @param subject A string or object of class XStringSet containing the subject
#' sequence/strin. Must be of length 1.
#' @param ... Any additional parameters are passed on to
#' [pwalign::pairwiseAlignment()]. This allows for adjusting alignment
#' algorithm and parameters.
#'
#' @return A list containing tibbles with information on mismatches and indels.
#'
#' @examples
#' query_seq <- Biostrings::DNAStringSet(c("ACCGTACCTGG", "ACCTTGG"))
#' subject_seq <- Biostrings::DNAStringSet("ACCGTACCGGG")
#' alignment_table(query_seq, subject_seq)
#'
#' # Works with any string
#' query_string <- c("boo", "fizzbuzz")
#' subject_string <- "boofizz"
#' alignment_table(query_string, subject_string)
#'
#' @import dplyr
#' @importFrom pwalign pairwiseAlignment mismatchTable width
#'
#' @export
alignment_table <- function(query = XStringSet,
                            subject = XStringSet,
                            ...){
  if(is.null(names(query))){
    names(query) <- paste0("query", 1:length(query))
  }
  if(any(is.na(names(query)))){
    names(query)[is.na(names(query))] <- paste0("query",
                                                which(is.na(names(query))))
  }
  if(is.null(names(subject))){
    names(subject) <- "subject"
  }

  #Align all seqs to subject, also subject itself
  all_sequences <- c(query, subject)
  alig <- pairwiseAlignment(all_sequences, subject, ...)

  ranges <- alig@subject@range %>% as.data.frame()

  #Extract insertions first to use in "end" adjustment for mismatches
  insertions <- alig@subject@indel %>% as.data.frame() %>%
    as_tibble() %>%
    left_join(tibble(group = 1:length(names(all_sequences)),
                     PatternName = names(all_sequences)), by = "group") %>%
    mutate(PatternName = factor(PatternName,
                                levels = c(names(subject),
                                           names(query))))
  #Extract mismatches, add empty rows if none detected
  mismatches <- alig %>%
    mismatchTable() %>%
    as_tibble()
  if(nrow(mismatches) == 0){
    mismatches <- add_row(mismatches,
                          PatternId = 1:length(names(all_sequences)))
  }
  mismatches <- full_join(mismatches,
            tibble(PatternId = 1:length(names(all_sequences)),
                     PatternName = names(all_sequences),
                     PatternLength = width(all_sequences),
                     ranges), by = join_by(PatternId)) %>%
    # start and end need to be adjusted so that each character is of width 1
    # centered on its number, e.g. character 4 is from 3.5-4.5
    mutate(start = start - 0.5,
           end = end + 0.5)
  # include the reverse mismatches against the subject to show the original
  # characters on the subject
  mismatches_ref <- dplyr::select(mismatches, contains("Subject")) %>%
    filter(!is.na(SubjectSubstring)) %>%
    dplyr::rename(PatternSubstring = "SubjectSubstring") %>%
    bind_cols(
      dplyr::select(mismatches, PatternId, PatternName,
             PatternLength, start, end, width) %>%
        filter(PatternName == names(subject))
    )
  #combine with mismatches
  mismatches <- bind_rows(mismatches, mismatches_ref) %>%
    mutate(PatternName = factor(PatternName,
                                levels = c(names(subject),
                                           names(query)))) %>%
    mutate(end = ifelse(end < (PatternLength-sum(insertions$width)+start),
                        (PatternLength-sum(insertions$width)+start), end)) %>%
    unique() %>%
    mutate(feature = "mismatch")

  # Extract deletions
  deletions <- alig@pattern@indel %>% as.data.frame() %>%
    as_tibble() %>%
    left_join(tibble(group = 1:length(names(all_sequences)),
                     PatternName = names(all_sequences)), by = "group") %>%
    mutate(PatternName = factor(PatternName,
                                levels = c(names(subject),
                                           names(query)))) %>%
    # start and end need to be adjusted so that each character deletion is of
    # width 1 centered on its number, e.g. character deletion 4 is from 3.5-4.5
    mutate(start = start - 0.5,
           end = end + 0.5,
           feature = "deletion")


    # Careful, subject@indel is relative to the query, unlike pattern@indel,
    # Adjust by merging with mismatch tibble and adding start.
    insertions <- rename(insertions, start_ins = "start") %>%
    left_join(dplyr::select(mismatches, start, PatternName),
              by = join_by(PatternName)) %>%
    mutate(end = end + start,
           start = start_ins-1+start,
           feature = "insertion")

    #Combine mismatches, deletions and insertions into one tbl
    out <- bind_rows(mismatches,
              dplyr::select(insertions, start, end, width,
                            PatternName, feature),
              dplyr::select(deletions, start, end, width,
                            PatternName, feature))
  return(out)
}
