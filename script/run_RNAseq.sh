#!/bin/bash
# Program:
# 	RNAseq workflow
# History:
# 	2020/09/22	First release
# Author:
# 	WangPF

##
. .conf

#########################################################################
## 准备
cd ${work_dir}/refseq
hisat2-build -p ${thread} ${genome} ${index}
gffread ${gff} -T -o ${gtf}
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

if [ strandSpecific = 0 ]; then
	hisat2 --new-summary -p ${thread} \
           -x ${index} \
           -1 ../00.data/01.clean_data/${i[1]}_1.clean.fastq.gz \
           -2 ../00.data/01.clean_data/${i[1]}_2.clean.fastq.gz \
           -S ${i[1]}.sam \
           1> ${i[1]}.log 2>&1
elif [ strandSpecific = 1 ]; then
	hisat2 --new-summary -p ${thread} \
           -x ${index} --rna-strandness RF \
           -1 ../00.data/01.clean_data/${i[1]}_1.clean.fastq.gz \
           -2 ../00.data/01.clean_data/${i[1]}_2.clean.fastq.gz \
           -S ${i[1]}.sam \
           1> ${i[1]}.log 2>&1
elif [ strandSpecific = 2 ]; then
	hisat2 --new-summary -p ${thread} \
           -x ${index} --rna-strandness FR \
           -1 ../00.data/01.clean_data/${i[1]}_1.clean.fastq.gz \
           -2 ../00.data/01.clean_data/${i[1]}_2.clean.fastq.gz \
           -S ${i[1]}.sam \
           1> ${i[1]}.log 2>&1
fi

samtools sort -@ ${thread} -o ${i[1]}.sort.bam ${i[1]}.sam
samtools index ${i[1]}.sort.bam

# 定量
mkdir -p  ${work_dir}/02.Quantification
cd ${work_dir}/02.Quantification
Rscript run-featurecounts.R \
	-b ../01.Mapping/${i[1]}.sort.bam \
	-g ${gtf} -o ${i[1]} \
	--nthread ${thread} \
	--attrType ${attrType} \
	--strandSpecific ${strandSpecific}

done

cd ${work_dir}/03.Merge_result

ls ../02.Quantification/*.count >genes.quant_files.txt
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
