#!/bin/sh
# Install dante into $PREFIX. noarch package; copies the source tree
# (including the bundled tool-data/protein_domains database) into
# $PREFIX/share/dante and symlinks the executables into $PREFIX/bin.
set -x -e

DANTE_DIR="${PREFIX}/share/dante"

mkdir -p "${PREFIX}/bin"
mkdir -p "${DANTE_DIR}"
cp -r . "${DANTE_DIR}"

# Executables
ln -s "${DANTE_DIR}/dante"                          "${PREFIX}/bin/dante"
ln -s "${DANTE_DIR}/dante_gff_output_filtering.py"  "${PREFIX}/bin/dante_gff_output_filtering.py"
ln -s "${DANTE_DIR}/coverage2gff.py"                "${PREFIX}/bin/coverage2gff.py"
ln -s "${DANTE_DIR}/dante_gff_to_dna.py"            "${PREFIX}/bin/dante_gff_to_dna.py"
ln -s "${DANTE_DIR}/fasta2database.py"              "${PREFIX}/bin/fasta2database.py"
ln -s "${DANTE_DIR}/parse_aln.py"                   "${PREFIX}/bin/parse_aln.py"
ln -s "${DANTE_DIR}/summarize_gff.R"                "${PREFIX}/bin/summarize_gff.R"
