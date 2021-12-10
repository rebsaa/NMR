#!/bin/bash

SPECIES=H.glaber
mkdir -p $SPECIES\/INPUT
for FILE in $(find ../Results/Trinity/$SPECIES -wholename '*/Trinity.fasta'); do

	SAMPLE=$(echo $FILE | cut -d '/' -f 5)
	PEP=../Results/Trinity/$SPECIES/$SAMPLE/TransDecoder/Trinity.fasta.transdecoder.pep
	CDS=../Results/Trinity/$SPECIES/$SAMPLE/TransDecoder/Trinity.fasta.transdecoder.cds

	##Extract the id of the isoform with the longest cds per gene
	> $SPECIES\/$SAMPLE\_longest_prediction.txt
	for GENE in $(grep '>' $CDS | cut -d ' ' -f1 | sed 's/_i[0-9].*//g;'| sort -u); do
	grep '>' $CDS | grep $GENE |cut -d ' ' -f 1,2 | sort -V -r -k2 | head -n 1 | cut -d ' ' -f1 >> $SPECIES\/$SAMPLE\_longest_prediction.txt
	done

	cat $FILE | tr '\n' '\t' | sed 's/>/\n>/g' | grep -wf <(sed 's/\.p[0-9]*//g' $SPECIES\/$SAMPLE\_longest_prediction.txt) | tr '\t' '\n' | sed "/^\s*$/d;s/^>/>${SAMPLE}\_/g" | cut -d ' ' -f 1 > $SPECIES\/$SAMPLE.fasta
	cat $PEP | tr '\n' '\t' | sed 's/>/\n>/g' | grep -wf $SPECIES\/$SAMPLE\_longest_prediction.txt | tr '\t' '\n' | sed "/^\s*$/d;s/^>/>${SAMPLE}\_/g" | cut -d ' ' -f 1 | sed 's/\.p[0-9]*$//g'> $SPECIES\/INPUT/$SAMPLE.pep
	cat $CDS | tr '\n' '\t' | sed 's/>/\n>/g' | grep -wf $SPECIES\/$SAMPLE\_longest_prediction.txt | tr '\t' '\n' | sed "/^\s*$/d;s/^>/>${SAMPLE}\_/g" | cut -d ' ' -f 1 | sed 's/\.p[0-9]*$//g'> $SPECIES\/$SAMPLE.cds

done
