library(Biostrings)

test_that("Check if basic alignment works", {
    query <- DNAStringSet("ACGT")
    subject <- DNAStringSet("ACGT")
    result <- alignment_table(query, subject)
    expect_equal(dim(result), c(2, 15))
    expect_true(all(is.na(result$PatternSubstring)))
})

test_that("Check if mismatch handling works", {
    query <- DNAStringSet("ACGT")
    subject <- DNAStringSet("AAAA")
    result <- alignment_table(query, subject)
    expect_true("mismatch" %in% unique(result$feature))
})

