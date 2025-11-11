hap_count=0
for f in simdat/denisovan_early/hmmix/decoded/decoded.NAMH_1.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_1.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_2.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_2.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_3.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_3.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_4.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_4.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_5.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_5.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_6.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_6.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_7.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_7.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_8.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_8.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_9.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_9.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_10.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_10.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_11.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_11.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_12.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_12.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_13.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_13.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_14.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_14.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_15.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_15.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_16.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_16.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_17.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_17.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_18.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_18.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_19.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_19.hap2.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_20.hap1.txt simdat/denisovan_early/hmmix/decoded/decoded.NAMH_20.hap2.txt;
do
echo $f
echo $hap_count
id=$(basename "$f" | cut -d. -f2) 
echo $id
hap=$(basename "$f" | cut -d. -f3 | sed 's/\.txt//')  # hap2
echo $hap
grep 'Archaic' "$f" | awk '($6>=0.8)' | awk -v id="$id" -v hap="$hap" -v count="$hap_count" 'BEGIN {OFS="\t"} {print count, $0, id, hap}' >> simdat/denisovan_early/fragments/denisovan_early_hmmix_all.bed
hap_count=$((hap_count+1))
done

for i in {1..20};
do
echo $i
out="simdat/denisovan_early/fragments/denisovan_early_hmmix_$i.bed"
echo $out
awk -v i="$i" '$2 == i {print}' simdat/denisovan_early/fragments/denisovan_early_hmmix_all.bed > "$out"
done