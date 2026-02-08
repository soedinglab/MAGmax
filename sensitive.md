### Sensitive mode

Sensitive mode selects representative genomes based on connectivity, defined as the highest sum of weighted degrees in the ANI graph. It applies a greedy maximum weighted dominating set algorithm, where weighted degree of a bin is calculated as the sum of positive differences between the pairwise ANI and the ANI threshold.

We compared N50 and gene counts between representative sets from sensitive and no-reassembly modes to assess whether representatives in sensitive mode exhibit systematic bias toward shorter sequences with inflated ANI to neighboring bins. For this, we used 20K genome set that was used in benchmarking magmax and galah.

Our results showed that there is no sytematic bias between two modes (upper and lower triangle values in figures show the percentage of representatives with higher N50 or gene counts relative to those selected by the other mode).

<img width="5319" height="2297" alt="n50_genecount_twomodecomp" src="https://github.com/user-attachments/assets/fe07ab10-d185-447e-b024-33c07a9d8e47" />

Taxonomic labels of the representatives and member genomes are consistent between both modes. The number of final replicated genomes, however, differs between the two modes, with 12976 genomes in `--sensitive` and 12935 in `--no-reassembly`. This difference is expected due to the use of different algorithms.
