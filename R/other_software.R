
micro_works={list(
  "fastp"=paste0('
  mkdir -p data/fastp
  ~/miniconda3/envs/waste/bin/fastp -w 8 -i data/rawdata/${sample}_f1.fastq.gz -o data/fastp/${sample}_1.fq.gz \\
    -I data/rawdata/${sample}_r2.fastq.gz -O data/fastp/${sample}_2.fq.gz -j data/fastp/${i}.json
'),
  "rm_human"=paste0('
  mkdir -p data/rm_human
  ~/miniconda3/envs/waste/bin/bowtie2 -p 8 -x ~/db/humangenome/hg38 -1 data/fastp/${sample}_1.fq.gz \\
    -2 data/fastp/${sample}_2.fq.gz -S data/rm_human/${sample}.sam \\
    --un-conc data/rm_human/${sample}.fq --very-sensitive
  rm data/rm_human/${sample}.sam
'),
  "kraken"=paste0('
  mkdir -p result/kraken
  python /share/home/jianglab/shared/krakenDB/K2ols/kraken2M.py -t 16 \\
      -i data/rm_human \\
      -c 0.05 \\
      -s .1.fq,.2.fq \\
      -o result/kraken \\
      -d /share/home/jianglab/shared/krakenDB/mydb2 \\
      -k ~/miniconda3/envs/waste/bin/kraken2 \\
      -kt /share/home/jianglab/shared/krakenDB/K2ols/KrakenTools
  mkdir -p result/kraken/kreport
  mv result/kraken/*_kreport.txt result/kraken/kreport/
  python ~/db/script/format_kreports.py -i result/kraken/kreport \\
      -kt /share/home/jianglab/shared/krakenDB/K2ols/KrakenTools --save-name2taxid
'),
  "kraken2"=paste0('
  mkdir -p result/kraken2
  ~/miniconda3/envs/waste/bin/kraken2 --threads 8 \\
    --db ~/db/kraken2/stand_krakenDB \\
    --confidence 0.05 \\
    --output result/kraken2/${sample}.output \\
    --report result/kraken2/${sample}.kreport \\
    --paired \\
    data/rm_human/${sample}.1.fq \\
    data/rm_human/${sample}.2.fq
'),
  "megahit"=paste0('
  #single sample
  mkdir -p result/megahit
  mkdir -p result/megahit/contigs
  ~/miniconda3/envs/waste/bin/megahit -t 8 -1 data/rm_human/${sample}.1.fq \\
    -2 data/rm_human/${sample}.2.fq \\
    -o result/megahit/${sample} --out-prefix ${sample}
  sed -i "/>/s/>/>${sample}_/" result/megahit/${sample}/${sample}.contigs.fa
  mv result/megahit/${sample}/${sample}.contigs.fa result/megahit/contigs/
'),
  "megahit2"=paste0('
  #multi-sample\u6df7\u62fc
  #\u9700\u8981\u5927\u5185\u5b58\uff0c\u5efa\u8bae3\u500dfq.gz\u5185\u5b58
  mkdir -p data/com_read
  mkdir -p result/megahit2

  for i in `cat samplelist`
  do
      echo ${i}
      cat data/rm_human/${i}.1.fq >> data/com_read/R1.fq
      cat data/rm_human/${i}.2.fq >> data/com_read/R2.fq
  done

  ~/miniconda3/envs/waste/bin/megahit -t 8 -1 data/com_read/R1.fq \\
    -2 data/com_read/R2.fq \\
    -o result/megahit2/ --out-prefix multi_sample
'),
  "prodigal"=paste0("
  mkdir -p result/prodigal
  mkdir -p result/prodigal/gene_fa
  mkdir -p result/prodigal/gene_gff
  ~/miniconda3/envs/waste/bin/prodigal -i result/megahit/contigs/${sample}.contigs.fa \\
      -d result/prodigal/gene_fa/${sample}.gene.fa \\
      -o result/prodigal/gene_gff/${sample}.gene.gff \\
      -p meta -f gff

  mkdir -p result/prodigal/fullgene_fa
  grep 'partial=00' result/prodigal/gene_fa/${sample}.gene.fa | cut -f1 -d ' '| sed 's/>//' > result/prodigal/gene_fa/${sample}.fullid
  seqkit grep -f result/prodigal/gene_fa/${sample}.fullid result/prodigal/gene_fa/${sample}.gene.fa > result/prodigal/fullgene_fa/${sample}.gene.fa

"),
  "cluster"=paste0('
  echo "start merge"
  cat result/prodigal/gene_fa/*.gene.fa > result/prodigal/all.gene.fa
  cat result/prodigal/fullgene_fa/*.gene.fa > result/prodigal/all.fullgene.fa
  echo "done"

  mkdir -p result/NR
  mmseqs easy-linclust result/prodigal/all.gene.fa result/NR/nucleotide mmseqs_tmp \\
    --min-seq-id 0.9 -c 0.9 --cov-mode 1  --threads 8

  # ~/miniconda3/envs/waste/bin/cd-hit-est -i result/prodigal/all.gene.fa \\
  #   -o result/NR/nucleotide.fa \\
  #   -aS 0.9 -c 0.9 -G 0 -g 0 -T 0 -M 0

  mv result/NR/nucleotide_rep_seq.fasta result/NR/nucleotide.fa
  rm result/NR/nucleotide_all_seqs.fasta
  ~/miniconda3/envs/waste/bin/seqkit translate result/NR/nucleotide.fa > result/NR/protein.fa
  sed \'s/\\*//g\' result/NR/protein.fa > result/NR/protein_no_end.fa
  rm result/NR/protein.fa
'),
  "seq_stat"=paste0('
  test_file=`head -n1 $samplelist`
  if [ -f result/megahit/contigs/${test_file}.contigs.fa ]; then
    ~/miniconda3/envs/waste/bin/seqkit stats result/megahit/contigs/*.contigs.fa >result/megahit/contig_stats
  fi
  if [ -f result/prodigal/gene_fa/${test_file}.gene.fa ]; then
    ~/miniconda3/envs/waste/bin/seqkit stats result/prodigal/gene_fa/*.gene.fa >result/prodigal/gene_fa_stats
  fi
  if [ -f result/prodigal/fullgene_fa/${test_file}.gene.fa ]; then
    ~/miniconda3/envs/waste/bin/seqkit stats result/prodigal/fullgene_fa/*.gene.fa >result/prodigal/fullgene_fa_stats
  fi
  if [ -f result/NR/nucleotide.fa ]; then
    ~/miniconda3/envs/waste/bin/seqkit stats result/NR/nucleotide.fa >result/NR/nucleotide_stat
  fi
'),
  "salmon-index"=paste0('
  ## \u5efa\u7d22\u5f15, -t\u5e8f\u5217, -i \u7d22\u5f15
  # \u5927\u70b9\u5185\u5b58
  mkdir -p result/salmon
  ~/miniconda3/envs/waste/share/salmon/bin/salmon index \\
    -t result/NR/nucleotide.fa \\
    -p 4 \\
    -i result/salmon/index
'),
  "salmon-quant"=paste0('
  mkdir -p result/salmon/quant
  ~/miniconda3/envs/waste/share/salmon/bin/salmon quant \\
      --validateMappings \\
      -i result/salmon/index -l A -p 4 --meta \\
      -1 data/rm_human/${sample}.1.fq \\
      -2 data/rm_human/${sample}.2.fq \\
      -o result/salmon/quant/${sample}.quant
'),
  "salmon-merge"=paste0("
  ls result/salmon/quant/*.quant/quant.sf |awk -F'/' '{print $(NF-1)}' |sed 's/.quant//' >tmp_finished
  diff_rows -f samplelist -s tmp_finished >tmp_diff
  # \u8ba1\u7b97\u7ed3\u679c\u7684\u884c\u6570
  line_count=$( cat tmp_diff| wc -l)

  # \u68c0\u67e5\u884c\u6570\u662f\u5426\u5927\u4e8e1\uff0c\u5982\u679c\u662f\u5219\u7ec8\u6b62\u811a\u672c
  if [ \"$line_count\" -gt 1 ]; then
    echo 'sf\u6587\u4ef6\u548csamplelist\u6570\u91cf\u4e0d\u4e00\u81f4'
    exit 1
  fi

  ##mapping rates
  cat result/salmon/quant/*.quant/logs/salmon_quant.log |grep 'Mapping rate = '|awk '{print $NF}'> tmp
  paste samplelist tmp > result/salmon/mapping_rates

  ## combine
  ~/miniconda3/envs/waste/share/salmon/bin/salmon quantmerge \\
      --quants result/salmon/quant/*.quant \\
      -o result/salmon/gene.TPM
  ~/miniconda3/envs/waste/share/salmon/bin/salmon quantmerge \\
      --quants result/salmon/quant/*.quant \\
      --column NumReads -o result/salmon/gene.count
  sed -i '1 s/.quant//g' result/salmon/gene.*
"),
  "eggnog"=paste0("
  ## \u5927\u70b9\u5185\u5b58, \u6570\u636e\u5e93\u670950G\u5de6\u53f3

  conda activate func
  ## diamond\u6bd4\u5bf9\u57fa\u56e0\u81f3eggNOG 5.0\u6570\u636e\u5e93, 1~9h\uff0c\u9ed8\u8ba4diamond 1e-3
  mkdir -p result/eggnog
  emapper.py --no_annot --no_file_comments --override \\
    --data_dir ~/db/eggnog5 \\
    -i result/NR/protein_no_end.fa \\
    --cpu 8 -m diamond \\
    -o result/eggnog/protein
  ## \u6bd4\u5bf9\u7ed3\u679c\u529f\u80fd\u6ce8\u91ca, 1h
  emapper.py \\
    --annotate_hits_table result/eggnog/protein.emapper.seed_orthologs \\
    --data_dir ~/db/eggnog5 \\
    --cpu 8 --no_file_comments --override \\
    -o result/eggnog/output

  ## \u6dfb\u8868\u5934, 1\u5217\u4e3aID\uff0c9\u5217KO\uff0c16\u5217CAZy\uff0c21\u5217COG\uff0c22\u5217\u63cf\u8ff0
  sed '1 i Name\\tortholog\\tevalue\\tscore\\ttaxonomic\\tprotein\\tGO\\tEC\\tKO\\tPathway\\tModule\\tReaction\\trclass\\tBRITE\\tTC\\tCAZy\\tBiGG\\ttax_scope\\tOG\\tbestOG\\tCOG\\tdescription' \\
    result/eggnog/output.emapper.annotations \\
    > result/eggnog/eggnog_anno_output
"),
  "cazy"=paste0("
  mkdir -p result/dbcan2
  diamond blastp   \\
  	--db ~/db/dbcan2/CAZyDB.07312020  \\
  	--query result/NR/protein_no_end.fa \\
  	--threads 8 -e 1e-5 --outfmt 6 \\
  	--max-target-seqs 1 --quiet \\
  	--out result/dbcan2/gene_diamond.f6
"),
  "rgi"=paste0("
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate rgi
  mkdir -p result/card
  ~/miniconda3/envs/rgi/bin/rgi main \\
    --input_sequence result/NR/protein_no_end.fa \\
    --output_file result/card/protein \\
    --input_type protein --num_threads 8 \\
    --clean --alignment_tool DIAMOND # --low_quality #partial genes
"),
  "vfdb"=paste0("
  mkdir -p result/vfdb
  diamond blastp   \\
  	--db ~/db/VFDB/VFDB_setB_pro \\
  	--query result/NR/protein_no_end.fa \\
  	--threads 8 -e 1e-5 --outfmt 6 \\
  	--max-target-seqs 1 --quiet \\
  	--out result/vfdb/gene_diamond.f6
"),
  "summary"=paste0("
  mkdir -p result/summ_table
  if [ -f result/eggnog/eggnog_anno_output ]; then
  # \u6c47\u603b\uff0c9\u5217KO\uff0c16\u5217CAZy\u6309\u9017\u53f7\u5206\u9694\uff0c21\u5217COG\u6309\u5b57\u6bcd\u5206\u9694\uff0c\u539f\u59cb\u503c\u7d2f\u52a0
  ~/db/script/summarizeAbundance.py \\
    -i result/salmon/gene.count \\
    -m result/eggnog/eggnog_anno_output \\
    -c '9,16,21' -s ',+,+*' -n raw \\
    -o result/summ_table/eggnog
  sed -i 's/^ko://' summ_table/eggnog.KO.raw.txt
  fi

  if [ -f result/card/protein.txt ]; then
  awk 'BEGIN {FS = \"\\t\"; OFS = \"\\t\"} {split($1, a, \" \"); $1 = a[1]} 1' result/card/protein.txt >result/card/protein_format_id.txt
  ~/db/script/summarizeAbundance.py \\
    -i result/salmon/gene.count \\
    -m result/card/protein_format_id.txt \\
    -c '11' -s ';' -n raw \\
    -o result/summ_table/card
  fi

  if [ -f result/dbcan2/gene_diamond.f6 ]; then
  # \u63d0\u53d6\u57fa\u56e0\u4e0edbcan\u5206\u7c7b\u5bf9\u5e94\u8868
  perl ~/db/script/format_dbcan2list.pl \
    -i result/dbcan2/gene_diamond.f6 \
    -o result/dbcan2/gene.list
  # \u6309\u5bf9\u5e94\u8868\u7d2f\u8ba1\u4e30\u5ea6
  ~/db/script/summarizeAbundance.py \\
    -i result/salmon/gene.count \\
    -m result/dbcan2/gene.list \\
    -c '2' -s ',' -n raw \\
    -o result/dbcan2/cazy
  fi

  if [ -f result/vfdb/gene_diamond.f6 ]; then
  sed -i '1 i Name\tvf\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' \\
    result/vfdb/gene_diamond.f6
  ~/db/script/summarizeAbundance.py \\
    -i result/salmon/gene.count \\
    -m result/vfdb/gene_diamond.f6 \\
    -c '2' -s ';' -n raw \\
    -o result/summ_table/card
  fi
")
)}

#' Microbiome sbatch
#'
#' @param work_dir work_dir
#' @param step "fastp","rm_human","megahit","prodigal","salmon-quant",...
#' @param all_sample_num all sample number
#' @param array array number
#' @param partition partition
#' @param mem_per_cpu mem_per_cpu, "2G"
#' @param cpus_per_task cpus_per_task
#'
#' @export
micro_sbatch=function(work_dir="/share/home/jianglab/pengchen/work/asthma/",
                      step="fastp",all_sample_num=40,array=1,
                      partition="cpu",cpus_per_task=1,mem_per_cpu="2G"){
  num_each_array=ceiling(all_sample_num/array)
  if(!endsWith(work_dir,"/"))work_dir=paste0(work_dir,"/")
  array=ifelse(array>1,paste0("1-",array),"1")

  header={paste0("#!/bin/bash
#SBATCH --job-name=",step,"
#SBATCH --output=",work_dir,"log/%x_%A_%a.out
#SBATCH --error=",work_dir,"log/%x_%A_%a.err
#SBATCH --array=",array,"
#SBATCH --partition=",partition,"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=",cpus_per_task,"
#SBATCH --mem-per-cpu=",mem_per_cpu,"
")}

  set={paste0('samplelist=',work_dir,'samplelist

echo start: `date +"%Y-%m-%d %T"`
start=`date +%s`

####################
')}

  set2={paste0('echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
START=$SLURM_ARRAY_TASK_ID

NUMLINES=',num_each_array,' #how many sample in one array

STOP=$((SLURM_ARRAY_TASK_ID*NUMLINES))
START="$(($STOP - $(($NUMLINES - 1))))"

#set the min and max
if [ $START -lt 1 ]
then
  START=1
fi
if [ $STOP -gt ',all_sample_num,' ]
then
  STOP=',all_sample_num,'
fi

echo "START=$START"
echo "STOP=$STOP"

####################
')}

  loop1={paste0('for (( N = $START; N <= $STOP; N++ ))
do
  sample=$(head -n "$N" $samplelist | tail -n 1)
  echo $sample')}

  if(step%in%names(micro_works))work={micro_works[[step]]}
  else work="
  do something"

  loop2="done
"

  end={paste0('
##############
echo end: `date +"%Y-%m-%d %T"`
end=`date +%s`
echo TIME:`expr $end - $start`s')}

  if(step=="kraken"){
    res_text=paste0("#!/bin/bash
#SBATCH --job-name=kraken2M
#SBATCH --output=",work_dir,"log/%x_%A_%a.out
#SBATCH --error=",work_dir,"log/%x_%A_%a.err
#SBATCH --time=14-00:00:00
#SBATCH --partition=mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=100G

fqp=",work_dir,"rm_human
python /share/home/jianglab/shared/krakenDB/K2ols/kraken2M.py -t 32 \\
    -i ${fqp} \\
    -c 0.05 \\
    -s _1.fastq,_2.fastq \\
    -o ",work_dir,"result/kraken/ \\
    -d /share/home/jianglab/shared/krakenDB/mydb2 \\
    -k ~/miniconda3/envs/waste/bin/kraken2 \\
    -kt /share/home/jianglab/shared/krakenDB/K2ols/KrakenTools")
  }

  if(step%in%c("fastp","rm_human","kraken2","megahit","prodigal","salmon-quant"))
    res_text=paste0(header,set,set2,loop1,work,loop2,end)
  else res_text=paste0(header,set,work,end)

  lib_ps("clipr", library = F)
  clipr::write_clip(res_text)
  message(res_text)
}

#' Prepare the result from fastp (.json file)
#'
#' @param jsonfiles the directory contains .json file
#' @param prefix default c("Raw","Clean"), for the before filtering and after filtering.
#'
#' @return data.frame
#' @export
#'
pre_fastp=function(jsonfiles,prefix=c("Raw","Clean")){
  lib_ps("jsonlite","dplyr",library = FALSE)
  result_list <- jsonfiles
  result_dict <- list()
  merge_result_dict <- list()
  if(length(prefix)<2)prefix=c("Raw","Clean")
  for (i in result_list) {
    result_dict[[i]] <- jsonlite::fromJSON(i)
  }

  for (k in names(result_dict)) {
    v <- result_dict[[k]]
    key <- strsplit(basename(k), "\\.")[[1]][1]
    merge_result_dict[[key]] <- data.frame(
      Raw_reads = v$summary$before_filtering$total_reads,
      Raw_bases = v$summary$before_filtering$total_bases,
      Raw_Q20 = v$summary$before_filtering$q20_rate * 100,
      Raw_Q30 = v$summary$before_filtering$q30_rate * 100,
      Raw_GC = v$summary$before_filtering$gc_content * 100,
      Clean_reads = v$summary$after_filtering$total_reads,
      Clean_bases = v$summary$after_filtering$total_bases,
      Clean_Q20 = v$summary$after_filtering$q20_rate * 100,
      Clean_Q30 = v$summary$after_filtering$q30_rate * 100,
      Clean_GC = v$summary$after_filtering$gc_content * 100,
      Duplication = v$duplication$rate * 100,
      Insert_size = v$insert_size$peak,
      Clean_reads_Raw_reads = v$filtering_result$passed_filter_reads / v$summary$before_filtering$total_reads * 100
    )
  }

  df<- do.call(rbind, merge_result_dict)
  colnames(df)=c(paste0(prefix[1],"_",c("reads","bases","Q20","Q30","GC")),
                 paste0(prefix[2],"_",c("reads","bases","Q20","Q30","GC")),
                 "Duplication","Insert_size",
                 paste0(prefix[2],"/",prefix[1]))
  df
}
