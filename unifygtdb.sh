#!/bin/bash
# Usage: bash unify.sh gtdb_taxonomy.tsv gtdbtk_species_representatives.tsv memberships.tsv > unified.tsv

GTDB_TAX="$1"
GTDBTK="$2"
MEMBERS="$3"

awk '
BEGIN { OFS="\t"; print "#user_representative\tgtdb_representative" }

# File 1: gtdb_taxonomy.tsv -> species -> GTDB reference genome
FILENAME==ARGV[1] {
    if ($0 ~ /^#/) next
    genome = $1
    split($0, a, "s__")
    species = a[length(a)]
    if (species != "" && !(species in gtdb_ref))
        gtdb_ref[species] = genome
    gtdb_species[species] = 1
    next
}

# File 2: gtdbtk_species_representatives.tsv -> your_rep -> species
FILENAME==ARGV[2] {
    if ($0 ~ /^#/) next
    your_rep = $1
    split($0, a, "s__")
    species = a[length(a)]
    bin_species[your_rep] = species
    next
}

# File 3: memberships.tsv -> process each representative
FILENAME==ARGV[3] {
    if ($0 ~ /^#/) next
    your_rep = $1
    species = (your_rep in bin_species) ? bin_species[your_rep] : ""
    if (species != "" && species in gtdb_ref) {
        print your_rep, gtdb_ref[species]   # user rep -> matched GTDB ref
        covered[species] = 1
    } else {
        print your_rep, "unknown"           # no GTDB species match
    }
    next
}

END {
    # GTDB species not covered by any user representative
    for (species in gtdb_ref) {
        if (!(species in covered)) {
            ref = gtdb_ref[species]
            print ref, ref
        }
    }
}
' "$GTDB_TAX" "$GTDBTK" "$MEMBERS"
