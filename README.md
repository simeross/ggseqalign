<img src="inst/hexlogo/hexlogo.png" alt="ggseqalign hex sticker" width="400"/>

## Intro
This is an R package to perform pairwise alignments of strings and plot them in a minimal style that is suitable for strings/sequences of any length. It is compatible with ggplot2 and DNA/AA sequence objects from Biostrings.

I am currently working on a Bioconductor submission for this package, but it should be usable in its current state.

## Installation
### R version 4.4 and above
`ggseqalign` can be installed from it's original source on GitLab (requires `devtools`)
```
devtools::install_git("https://gitlab.com/nmbu.no/ipv/lim-rossmann/ggseqalign.git")
```
### earlier R versions
If you use any R version older than `4.4.0`, install version 0.1.0 of this package variants
```
devtools::install_git("https://gitlab.com/nmbu.no/ipv/lim-rossmann/ggseqalign@e1a2e99ee7d079d2e4331bc101382cfebe125fa8")
```
### Bioconductor
This will be updated for installation instructions from Bioconductor if the submission is successful.

## Quick start
All you need is this package and some strings to align and you are ready to go.
```
library(ggseqalign)

query_strings <- (c("boo", "fibububuzz", "bozz", "baofuzz"))
subject_string <- "boofizz"

alignment <- alignment_table(query_strings, subject_string)

plot_sequence_alignment(alignment)
```

To align DNA or AA sequences from a fasta file, read them in with `Biostrings`
```
library(ggseqalign)
library(Biostrings)

query_sequences <- Biostrings::readDNAStringSet("my_multi_sequence_fasta.fa")
subject_sequence <- Biostrings::readDNAStringSet("my_reference_fasta.fa")

alignment <- alignment_table(query_sequences, subject_sequence)

pl <- plot_sequence_alignment(alignment)
pl
```

To style the plot generated above to your own taste, use `ggplot2`, for example:
```
library(ggplot2)
 
pl +
 ylab("Sequence variants") +
 xlab("Length in bp") +
 scale_color_viridis_d() +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_text())
```

## Bug reports
If you come across bugs, please [submit an issue](https://gitlab.com/nmbu.no/ipv/lim-rossmann/ggseqalign/-/issues)

### License
This package is licenced under the Artistic License 2.0.