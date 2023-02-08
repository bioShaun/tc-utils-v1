fasta=$1

scriptDir="/public/scripts/pyTools/fasta"

seqkit fx2tab -g -l ${fasta} > ${fasta}.gc.txt

awk '$3 <= 120' ${fasta}.gc.txt |cut -f1,2,5 > part1.seq.gc

awk '$3 > 120 && $3 <= 200' ${fasta}.gc.txt |cut -f1 > part2.seq.txt

seqkit fx2tab ${fasta}| grep -f part2.seq.txt | seqkit tab2fx | seqkit sliding -s 1 -W 120 - | seqkit fx2tab -g -l - > part2.seq.raw.gc

python "${scriptDir}/selectGC50.py" part2.seq.raw.gc part2.seq.gc   

awk '$3 > 200' ${fasta}.gc.txt |cut -f1 > part3.seq.txt

seqkit fx2tab ${fasta}| grep -f part3.seq.txt | seqkit tab2fx | seqkit sliding -s 1 -W 120 - | seqkit fx2tab -g -l - > part3.seq.raw.gc

python "${scriptDir}/selectStartEnd.py" part3.seq.raw.gc part3.seq.gc

cat *seq.gc | sort > final.${fasta}.gc.txt
