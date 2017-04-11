#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $DIR
test_data="$DIR/test_data"
classification_tbl='!!!! set up the path !!!!!'
pdb='!!! set up the path'

######## Protein Domais Finder
## single_seq, for/rev strand of mapping
$DIR/protein_domains_pd.py -q $test_data/GEPY_test_long_1 -pdb $pdb -cs $classification_tbl -dir $PWD/tmp/single_fasta/
## multifasta
$DIR/protein_domains_pd.py -q $test_data/vyber-Ty1_01.fasta -pdb $pdb -cs $classification_tbl -dir $PWD/tmp/multifasta/
## multifasta_win
$DIR/protein_domains_pd.py -q $test_data/vyber-Ty1_01.fasta -pdb $pdb -cs $classification_tbl -wd 3100 -od 1500 -dir $PWD/tmp/multifasta_win

## testing if outputs are the same in case of using sliding window and not
if [[ $(diff $PWD/tmp/multifasta/output_domains.gff $PWD/tmp/multifasta_win/output_domains.gff) -eq 0 ]];then
	echo "Testing output of sliding window comparing to no window accomplished sucessfuly"
else
	echo "WARNING! There is difference between outputs of sliding window and no window used"
fi
######## Protein Domains Filter
## default params
$DIR/domains_filtering.py -dom_gff $PWD/tmp/single_fasta/output_domains.gff 
if [[ -e $PWD/tmp/single_fasta/domains_filtered.gff ]] && [[ -e $PWD/tmp/single_fasta/dom_prot_seq.txt ]] ; then
	echo -e "Filtered file and protein seqs file for default parameters exists"
else
	echo -e "Filtered outputs for default parameters are missing"
fi
if [[ $(cat $PWD/tmp/single_fasta/domains_filtered.gff | wc -l) -gt 1 ]];then
	echo "File was correctly filtered using default parameters"
fi
## Ty1-RT filtering
$DIR/domains_filtering.py -dom_gff $PWD/tmp/multifasta/output_domains.gff -sd Ty1-RT
if [[ -e $PWD/tmp/multifasta/domains_filtered.gff ]] && [[ -e $PWD/tmp/multifasta/dom_prot_seq.txt ]]; then
	echo -e "Filtered file and protein seqs file of Ty1-RT domains exist"
else
	echo -e "Filtered outputs of Ty1-RT domains are missing"
fi
if [[ $(cat $PWD/tmp/multifasta/domains_filtered.gff | wc -l) -gt 1 ]];then
	echo "File was correctly filtered for Ty1-RT domains"
fi
