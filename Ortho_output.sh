#Determine the longest peptide sequence per orthogroup, run in OrthoFinder/$SPECIES/INPUT/OrthoFinder/Results*/Orthogroup_Sequences

for FILE in *fa; do 
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' $FILE | paste - - | awk '{print $0 "\t" length($2)}'| sort -Vrk 3 | head -n 1 | cut -f1 >> longest_peptides.txt
done
