#!/bin/bash
#Determine the longest peptide sequence per orthogroup, run in OrthoFinder/$SPECIES/INPUT/OrthoFinder/Results*/Orthogroup_Sequences

SPECIES=H.glaber
>$SPECIES/$SPECIES\_longest_peptide_per_orthogroup.pep

for FILE in $(find . -wholename './H.glaber/INPUT/OrthoFinder/*/Orthogroup_Sequences/*fa'); do 
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' $FILE | paste - - | awk '{print $0 "\t" length($2)}'| sort -Vrk 3 | head -n 1 | cut -f1 >> $SPECIES\/longest_peptide_per_orthogroup.pep
 
done

for SAMPLE in $(cut -d '_' -f 1-3 $SPECIES\/$SPECIES\_longest_peptide_per_orthogroup.pep | sed 's/>//g' | sort -u); do
grep -wf $SPECIES\/$SPECIES\_longest_peptide_per_orthogroup.pep <(cat $SPECIES\/$SAMPLE\.fasta | paste - -) | tr '\t' '\n' > $SPECIES\/$SAMPLE\_longest_transcript_per_orthogroup.fasta
done

cat $SPECIES\/*longest_transcript_per_orthogroup.fasta > $SPECIES\/$SPECIES\_reference_transcript.fasta
