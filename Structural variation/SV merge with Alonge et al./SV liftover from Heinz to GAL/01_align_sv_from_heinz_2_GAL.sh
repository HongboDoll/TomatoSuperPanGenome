#!/bin/bash

rm -rf shell; mkdir shell
grep -v '#' merged.ont.v1.0.vcf | awk '$5!~/\[/&&$5!~/\]/' | while read i
do
chr=`echo $i | awk '{print $1}'`
pos=`echo $i | awk '{print $2}'`
start=`expr $pos - 200`
end=`expr $pos + 200`

echo -e """#!/bin/bash
samtools faidx heinz_ref/${chr}.fa ${chr}:${start}-${end} > query_seq/${chr}_${pos}_query.fa
nucmer -t 1 gal_ref/${chr}.fa query_seq/${chr}_${pos}_query.fa --prefix mummer_results/${chr}_${pos}
delta-filter -1 mummer_results/${chr}_${pos}.delta > mummer_results/${chr}_${pos}.f.delta
show-coords -THrcl mummer_results/${chr}_${pos}.f.delta > mummer_results/${chr}_${pos}.filter.1coords

""" > shell/${chr}_${pos}.sh1 && chmod 755 shell/${chr}_${pos}.sh1
done

rm command
rm -rf mummer_results && mkdir mummer_results
ls shell | while read i
do
	echo -e "shell/${i}" >> command
done

split -a 3 -d -l 580 command command_ && chmod 755 command*

ls command_* | while read i
do
	cat <(echo -e "#!/bin/bash") $i > t && mv t $i
	sbatch -J ${i} -p queue1 --qos=queue1 -N 1 --ntasks-per-node=1 -e %x.err -o %x.out "./$i" 
done


