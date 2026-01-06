### Comparison of MAGmax and Galah for genome dereplication

Galah is a fast dereplication tool and provides option to use CheckM2 quality scores similar to MAGmax. Here, we compare galah and magmax to evalute their performance.

Briefly, Galah’s clustering algorithm proceeds as follows:

1) genomes are ordered by quality scores (highest is first with index 0)

2) preclusters (or connected components) are obtained using single-linkage clustering

3) within each precluster, the first index genome (highest quality genomes in the preculster) is selected as a representative. Genomes are itereted in index order. If any genome in the iteration has ANI < threshold with the repesentative genome, it becomes a new cluster representative. Likewise, all cluster representatives are found based on the condition that ANI < threshold with existing representatives.

4) non-representative bins are assigned to the closest representative based on their ANI values.

In this way, galah ensures that (i) each representative is <99% ANI to each other representative, and (ii) all members are >=99% ANI to the representative.


### How does MAGmax differ from Galah?

Within each precluster component, MAGmax finds clusters using maximial cliques identification algorithm, i.e., genomes in the clusters share ANI > threshold with all cluster members.

A genome can be part of multiple clusters such that any pair that shares ANI > threshold should be in its cluster and any pair with ANI < threshold should be in different clusters.

If reassembly mode is enabled, MAGmax merges bins within the cluster and reassembles them to improve completeness and reduce contamination. Otherwise, it selects the best quality genome as representative.


### Demonstration using a simple real dataset


#### **Input:** 9 real high-quality (>90% completeness, <5% contamination) bacterial genomes. Indexed 0-8.

        Index   Genomefilename  Completeness    Contamination   quality_score (completeness - 5*contamination)
        0   GCA_948538755.1_MGBC165981  94.21   0.09    93.76
        1   GCA_948610485.1_MGBC109371  90.77   0.17    89.92
        2   GCA_948908505.1_MGBC155662  96.27   0.09    95.82
        3   GCA_948919045.1_MGBC146531  90.67   0.05    90.42
        4   GCA_948681465.1_MGBC155251  97.83   0.03    97.68
        5   GCMeta_00504745  90.97	0.11    90.42
        6   GCMeta_00505367  97.19	0.07    96.84
        7   GCMeta_01449380  96.43	0.04    96.23
        8   GCMeta_01454054  90.24  0.05    89.99

#### ANI threshold for depreplication: **99.9%**

Skani pairs with ANI >= 99.9%:

        0 - 1, 2, 3, 5, 6, 7, 8

        1 - 0, 2, 3, 4, 5, 6, 7, 8

        2 - 0, 1, 3, 4, 5, 6, 7, 8

        3 - 0, 1, 2, 4, 5, 6, 7, 8

        4 - 1, 2, 3, 5, 6, 7, 8

        5 - 0, 1, 2, 3, 4, 6, 7, 8

        6 - 0, 1, 2, 3, 4, 5, 7

        7 - 0, 1, 2, 3, 4, 5, 6, 8

        8 - 0, 1, 2, 3, 4, 5, 7

Only genomic pairs *0-4* and *6-8* share ANI < 99.9%, while all other pairs share ANI >= 99.9%.

**Galah result:** 2 genomes {indexed 0, 4}

**MAGmax result:** 1 genome {indexed 4} *(without reassembly mode)*


#### Why do the numbers of dereplicated bins differ?

1. Both tools, form the same precluster with all 9 genomes in one connected component. 

2. Galah sorts the genomes by quality score and selects the best genome (index 4) as representative. It iterates through the sorted genome list in the component and finds that genome at index 0 shares ANI < threshold with already selected representative. Genome at index 0 becomes a new representative, even though there is another genome at index 6 with >99.9% ANI to the new representative has better quality score.

3. MAGmax identifies 4 clusters using maximal cliques algorithm which keeps genomes 0-4 and 6-8 in different clusters ({0, 1, 2, 3, 5, 6, 7}, {0, 1, 2, 3, 5, 7, 8}, {1, 2, 3, 4, 5, 6, 7} & {1, 2, 3, 4, 5, 7, 8}).

4. MAGmax considers **all pairs with 99.9% ANI** and chooses the best quality genome between the pairs, whereas Galah ingores many such pairs.


<span style="color:blue"> **Overall, MAGmax algorithmically stricter than Galah and ensures that always best quality genome in the genomic pair is selected as representative.**</span>


### Computational efficiency

Runtime and peak memory usage of magmax and galah were measured on three datasets with different numbers of input genomes (9, 500, and 20K).
For each dataset, the best performance of the two tools was compared based on three replicate runs under the same CPU and memory limits.

Galah runs faster on the small dataset, whereas no significant differences in runtime or peak memory usage were observed between the two tools on larger datasets.

MAGmax also supports using a precomputed ANI file as direct input, reducing the runtime for the 20K dataset from 4,522 seconds to 357 seconds (~12× speedup). This enables iterative dereplication across different cutoffs at a much faster pace.

<img width="7119" height="2930" alt="computational_efficiency" src="https://github.com/user-attachments/assets/cc9c8b46-fb24-442a-9ab8-4c38ff9759c1" />

