### Generate custom species-level database

Genome dereplication is not always perfect due to inherent limitation of hierarchical clustering algorithms used in derpelication tools (dRep and Galah). Alternatively, taxonomic classification using GTDBtk followed by grouping genomes by taxonomy assignment is other option for dereplication, but it has limitations too. 1) ANI radius of under represented species may be inaccurate causing wrong taxonomy labeling. 2) Novel species canot be assigned. Combining dereplication and taxonomic classification can enhance the discovery of novel species with improved accuracy.

