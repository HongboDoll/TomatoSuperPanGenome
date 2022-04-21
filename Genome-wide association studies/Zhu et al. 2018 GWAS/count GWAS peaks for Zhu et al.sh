#!/bin/bash

rm -rf results; mkdir results
head -1 362_annotated_metabolite_321_sample_phenotype.xls | cut -f '2-9999' | sed 's/\t/\n/g' | grep -v '^$' | while read i
do 
	name=`awk '$1=="'"$i"'"' 362_annotated_metabolite_ID_name_corres.xls | awk -F '\t' '{print $2}'`
	emmax -v -d 10 -t 321_tomato_SV -p phenotype/${i}_phenotype.xls_order.xls -k PAV.aBN.kinf -c PCA_old.txt -o results/$i
	cat 321_tomato_SV.map | paste - results/${i}.ps | awk '{print $2,$1,$4,$8}' | sed '1i\SNP CHR BP P' | sed 's/ /\t/g' > results/${i}.321_tomato_SV.out.txt
	./manhattan_qq.R results/${i}.321_tomato_SV.out.txt results/${i}.321_tomato_SV.out.txt.png 1.53E-5 $name
done


ls results/*txt | while read i
do
	name=`echo $i | awk -F '/' '{print $2}' | awk -F '.' '{print $1}'`
	awk '$4<1.53E-5' $i | ./output_peak_lead_SV.py 1000000 results/${name}.PAV.1MbNonredundant.xls results/${name}.PAV.1MbWindow.count
done

cat <(echo -e "CHR\tPOS") <(awk '{print $1,$2}' results/SlFM2000.PAV.1MbWindow.count) > coordinates
ls results/*count | while read i
do
	name=`echo $i | awk -F '/' '{print $2}' | awk -F '.' '{print $1}'`
	cat <(echo -e "${name}") <(awk '{print $3}' $i) > ${i}.Count
done

paste coordinates results/*Count > 362_annotated_metabolite_1MbWindow_count.xls

./count_8classes_metabolites_peak.py 362_annotated_metabolite_1MbWindow_count.xls ./362_annotated_metabolite_ID_name_corres.xls_4_plot > 362_annotated_metabolite_8classes_1MbWindow_count.xls

for i in `seq 1 12`
do
	awk '$1=="'"$i"'"' 362_annotated_metabolite_8classes_1MbWindow_count.xls > 362_annotated_metabolite_8classes_1MbWindow_count_chr${i}.xls 
	./plot_each_chromosome.R 362_annotated_metabolite_8classes_1MbWindow_count_chr${i}.xls t${i}.pdf ${i}
	awk '{print $1,$2,$3,$4+$5+$6+$7+$8+$9+$10+$11}' 362_annotated_metabolite_8classes_1MbWindow_count_chr${i}.xls > 362_annotated_metabolite_8classes_1MbWindow_count_ALL_chr${i}.xls
done

cat *_count_ALL_* > 362_annotated_metabolite_8classes_1MbWindow_count_ALL.xls

