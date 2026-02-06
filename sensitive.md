### Sensitive mode

Sensitive mode selects representative genomes based on connectivity, defined as the highest sum of weighted degrees in the ANI graph. It applies a greedy maximal weighted dominating set algorithm, where weighted degree of a bin is calculated as the sum of positive differences between the pairwise ANI and the ANI threshold.

We compared N50 and gene counts between representative sets from sensitive and no-reassembly modes to assess whether representatives in sensitive mode exhibit systematic bias toward shorter sequences with inflated ANI to neighboring bins.

Our results showed that there is no sytematic bias between two modes.

Sensitive mode is suitable if you want to generate custom reference databases for taxonomic assignment.
