#!/bin/bash
# Download the IMGT GENE-DB reference sequences (all functionalities: F, ORF, and all P), and keep the IG V/D/J-REGION entries for each species we use.
# These are only used for their functionality labels (see data/germlines/write-functionalities.py): the default germline sets themselves come from other places (see data/germlines/README.md).
# NOTE that data/germlines/imgt-download/human/ is a different, smaller IMGT download ("F+ORF+in-frame P"), which is missing most pseudogenes.
set -eu
outdir=$(dirname $(realpath $0))
url='https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP'
tmpfile=$(mktemp)
curl -sS -o $tmpfile "$url"
for spstr in 'human:Homo sapiens' 'macaque:Macaca mulatta'; do
    species=${spstr%%:*}
    imgt_species=${spstr#*:}
    # header: >accession|gene|species|functionality|region|...  NOTE the species field can have a suffix for the individual/strain (e.g. 'Macaca mulatta_17573'), so we match the start of it
    awk -v sp="$imgt_species" 'BEGIN{FS="|"} /^>/{keep = (($3 == sp || index($3, sp "_") == 1) && $2 ~ /^IG[HKL][VDJ]/ && $5 ~ /^[VDJ]-REGION$/)} keep' $tmpfile > $outdir/$species.fasta
    echo "$species: $(grep -c '>' $outdir/$species.fasta) sequences"
done
rm $tmpfile
echo "downloaded $(date +%Y-%m-%d) from $url" > $outdir/download-info.txt
