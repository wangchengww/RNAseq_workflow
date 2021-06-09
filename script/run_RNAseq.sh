#  run_RNAseq.sh
#  
#  Copyright 2021 WangPF <wangpf0608@126.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
#!/bin/bash
# Program:
# 	RNAseq workflow
# History:
# 	20210118	First release
# Author:
# 	WangPF

## 加载配置文件
. .conf

#########################################################################
## 准备
cd ${work_dir}/db
hisat2-build -p ${thread} ${genome} ${index}
gffread ${gff} -T -o ${gtf}
convert2bed --input=gff --output=bed < ${gff} > ${bed}
#########################################################################


IFS_OLD=$IFS
IFS=$'\n'

for i in $(cat ${sample})
do

IFS=$'\t'
i=($i)
IFS=$IFS_OLD

# 过滤、质控
mkdir -p ${work_dir}/00.data/01.clean_data/QC ${work_dir}/00.data/QC
cd ${work_dir}/00.data/01.clean_data

fastp -i ${i[2]} -o ./${i[1]}_1.clean.fastq.gz \
      -I ${i[3]} -O ./${i[1]}_2.clean.fastq.gz \
      --json=./${i[1]}.json --html=${i[1]}.html --report_title="${i[1]} fastp report" \
      --thread=8

# 比对
mkdir -p ${work_dir}/01.Mapping
cd ${work_dir}/01.Mapping

if [ $strandSpecific = 0 ]; then
	hisat2 --new-summary -p ${thread} \
           -x ${index} \
           -1 ../00.data/01.clean_data/${i[1]}_1.clean.fastq.gz \
           -2 ../00.data/01.clean_data/${i[1]}_2.clean.fastq.gz \
           -S ${i[1]}.sam \
           1> ${i[1]}.log 2>&1
elif [ $strandSpecific = 1 ]; then
	hisat2 --new-summary -p ${thread} \
           -x ${index} --rna-strandness RF \
           -1 ../00.data/01.clean_data/${i[1]}_1.clean.fastq.gz \
           -2 ../00.data/01.clean_data/${i[1]}_2.clean.fastq.gz \
           -S ${i[1]}.sam \
           1> ${i[1]}.log 2>&1
elif [ $strandSpecific = 2 ]; then
	hisat2 --new-summary -p ${thread} \
           -x ${index} --rna-strandness FR \
           -1 ../00.data/01.clean_data/${i[1]}_1.clean.fastq.gz \
           -2 ../00.data/01.clean_data/${i[1]}_2.clean.fastq.gz \
           -S ${i[1]}.sam \
           1> ${i[1]}.log 2>&1
fi

samtools sort -@ ${thread} -o ${i[1]}.sort.bam ${i[1]}.sam
samtools index ${i[1]}.sort.bam
rm ${i[1]}.sam

# 定量
mkdir -p  ${work_dir}/02.Quantification
cd ${work_dir}/02.Quantification
Rscript run-featurecounts.R \
	-b ../01.Mapping/${i[1]}.sort.bam \
	-g ${gtf} -o ${i[1]} \
	--nthread ${thread} \
	--featureType ${featureType} \
	--attrType ${attrType} \
	--strandSpecific ${strandSpecific}

done

cd ${work_dir}/00.data/00.raw_data
fastqc -o ./QC --nogroup --threads ${thread} *[fastq\|fq].gz
cd ${work_dir}/00.data/01.clean_data
fastqc -o ./QC --nogroup --threads ${thread} *clean.fastq.gz
Rscript stat.R

cd ${work_dir}/01.Mapping
echo -e "Sample,Total read,Mapping read,Mapping rate,Unique mapping read,Unique mapping rate" > align_stat.csv
for i in $(cut -f2 ${sample}); do perl alignStat.pl $i; done >> align_stat.csv

cd ${work_dir}/03.Merge_result

cut -f2 ${sample} | sed 's/^/..\/02.Quantification\//' | sed 's/$/.count/' > genes.quant_files.txt
perl script/abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files genes.quant_files.txt --out_prefix genes

# 差异表达
cd ${work_dir}/04.DE_analysis
perl ${trinity}/Analysis/DifferentialExpression/run_DE_analysis.pl \
    --matrix ../03.Merge_result/genes.counts.matrix \
    --method ${de_method} \
    --samples_file ../00.data/samples.txt \
    --contrasts contrasts.txt 

# 富集分析
cd ${work_dir}/05.GO_KEGG
for de_result in ${work_dir}/04.DE_analysis/DESeq2.*/*DE_results
do
Rscript enrich.R --de_result ${de_result} \
	--de_log2FoldChange ${de_log2FoldChange} \
	--de_padj ${de_padj} \
	--enrich_pvalue ${enrich_pvalue} \
	--enrich_qvalue ${enrich_qvalue} \
	--orgdb ${orgdb} \
	--species ${species}
done

cd ${work_dir}/06.GSEA
for de_result in ${work_dir}/04.DE_analysis/DESeq2.*/*DE_results
do
Rscript GSEA.R --de_result ${de_result} \
	--enrich_pvalue ${enrich_pvalue} \
	--orgdb ${orgdb} \
	--draePdf ${pdf} \
	--species ${species}
done
