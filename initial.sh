mkdir ../BUSCO
cd ../BUSCO

wget https://busco-data.ezlab.org/v4/data/lineages/mammalia_odb10.2020-09-10.tar.gz
tar -xf mammalia_odb10.2020-09-10.tar.gz 
rm mammalia_odb10.2020-09-10.tar.gz

mkdir ../SILVA
mkdir ../SILVA/SILVA_rRNA
cd ../SILVA

#Download fasta files for lsu and ssu rRNA
wget https://ftp.arb-silva.de/current/Exports/SILVA_138.1_LSUParc_tax_silva.fasta.gz
wget https://ftp.arb-silva.de/current/Exports/SILVA_138.1_SSUParc_tax_silva.fasta.gz

#Add origin to sequence identifier
zcat SILVA_138.1_LSUParc_tax_silva.fasta.gz | sed '/^>/s/$/ lsu/g' |  gzip > lsu.tmp
zcat SILVA_138.1_SSUParc_tax_silva.fasta.gz | sed '/^>/s/$/ ssu/g' | gzip > ssu.tmp

#Concatenate both files and remove blanks, in lines containing rRNA substitute U with T
zcat ssu.tmp lsu.tmp | sed 's/ /_/g' | sed '/^>/!s/U/T/g' | gzip > silva.fasta.gz

#Combine sequence and identifier in one line, sort for unique sequences, return to original format
zcat silva.fasta.gz | sed -e '/^>/s/$/@/' -e 's/^>/§/' | tr -d '\n' | tr "§" "\n" | tr "@" "\t" | sort -u -f -k2 | sed -e 's/^/>/' -e 's/\t/\n/' | gzip > unique_silva.fasta.gz

rm *tmp
rm SILVA_138.1*

#build bowtie2 index
bowtie2-build --large-index --threads 6 unique_silva.fasta.gz SILVA_rRNA/SILVA_rRNA

#prepare Diamond database

mkdir Diamond
cd Diamond

for ENS in $(cat diamond_ref_pep.txt); do
echo $ENS
wget $ENS
FILE=$(echo $ENS | cut -d '/' -f 9)
echo $FILE
NAME=$(echo $FILE | cut -d '.' -f 1)
diamond makedb --in $FILE -d $NAME --ignore-warnings
rm $FILE
done
