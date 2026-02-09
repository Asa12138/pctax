# This is a R script to run other software for microbiome analysis.
micro_works <- {
  list(
    "fastp" = paste0("
  mkdir -p data/fastp
  ~/miniconda3/envs/waste/bin/fastp -w 8 -i data/rawdata/${sample}_f1.fastq.gz -o data/fastp/${sample}_1.fq.gz \\
    -I data/rawdata/${sample}_r2.fastq.gz -O data/fastp/${sample}_2.fq.gz -j data/fastp/${sample}.json
"),
    "rm_human" = paste0("
  mkdir -p data/rm_human
  ~/miniconda3/envs/waste/bin/bowtie2 -p 8 -x ~/db/humangenome/hg38 -1 data/fastp/${sample}_1.fq.gz \\
    -2 data/fastp/${sample}_2.fq.gz -S data/rm_human/${sample}.sam \\
    --un-conc data/rm_human/${sample}.fq --very-sensitive
  rm data/rm_human/${sample}.sam
"),
    "kraken" = paste0("
  mkdir -p result/kraken
  python /share/home/jianglab/shared/krakenDB/K2ols/kraken2M.py -t 16 \\
      -i data/rm_human \\
      -m paired-end \\
      -c 0.05 \\
      -s .1.fq,.2.fq \\
      -o result/kraken \\
      -d /share/home/jianglab/shared/krakenDB/mydb2 \\
      -k ~/miniconda3/envs/waste/bin/kraken2 \\
      -kt /share/home/jianglab/shared/krakenDB/K2ols/KrakenTools
  mkdir -p result/kraken/kreport
  mv result/kraken/*_kreport.txt result/kraken/kreport/
  python /data/home/jianglab/share/pc_DB/format_kreports-PC.py -i result/kraken/kreport \\
      -kt /share/home/jianglab/shared/krakenDB/K2ols/KrakenTools --save-name2taxid
  mv format_report result/kraken/
"),
    "kraken2" = paste0("
  mkdir -p result/kraken2
  ~/miniconda3/envs/waste/bin/kraken2 --threads 8 \\
    --db ~/db/stand_krakenDB \\
    --confidence 0.05 \\
    --output result/kraken2/${sample}.output \\
    --report result/kraken2/${sample}.kreport \\
    --paired \\
    data/rm_human/${sample}.1.fq \\
    data/rm_human/${sample}.2.fq
"),
    "humann" = paste0("
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate biobakery3

  mkdir -p data/paired
  cat data/rm_human/${sample}.1.fq data/rm_human/${sample}.2.fq > data/paired/${sample}_paired.fq

  mkdir -p result/humann
  ~/miniconda3/envs/biobakery3/bin/humann --input data/paired/${sample}_paired.fq \\
    --output result/humann/ --threads 24

  ln result/humann/${sample}_paired_humann_temp/${sample}_paired_metaphlan_bugs_list.tsv result/humann/
  rm -rf result/humann/${sample}_paired_humann_temp
"),
    "megahit" = paste0('
  #single sample
  mkdir -p result/megahit
  mkdir -p result/megahit/contigs
  ~/miniconda3/envs/waste/bin/megahit -t 8 -1 data/rm_human/${sample}.1.fq \\
    -2 data/rm_human/${sample}.2.fq \\
    -o result/megahit/${sample} --out-prefix ${sample}
  sed -i "/>/s/>/>${sample}_/" result/megahit/${sample}/${sample}.contigs.fa
  mv result/megahit/${sample}/${sample}.contigs.fa result/megahit/contigs/
'),
    "metaspades" = paste0('
  #single sample
  mkdir -p result/metaspades
  mkdir -p result/metaspades/contigs
  ~/miniconda3/envs/metawrap/bin/metaspades.py -t 8 -k 21,33,55,77,99,127 --careful \
    -1 data/rm_human/${sample}.1.fq \
    -2 data/rm_human/${sample}.2.fq \
    -o result/metaspades/${sample}

  sed -i "/>/s/>/>${sample}_/" result/metaspades/${sample}/scaffolds.fasta
  mv result/metaspades/${sample}/scaffolds.fasta result/metaspades/contigs/${sample}_scaffolds.fasta
'),
    "megahit2" = paste0("
  #multi-sample\u6df7\u62fc
  #\u9700\u8981\u5927\u5185\u5b58\uff0c\u5efa\u8bae3\u500dfq.gz\u5185\u5b58
  mkdir -p data/com_read
  rm -r result/megahit2

  for i in `cat samplelist`
  do
      echo ${i}
      cat data/rm_human/${i}.1.fq >> data/com_read/R1.fq
      cat data/rm_human/${i}.2.fq >> data/com_read/R2.fq
  done

  ~/miniconda3/envs/waste/bin/megahit -t 8 -1 data/com_read/R1.fq \\
    -2 data/com_read/R2.fq \\
    -o result/megahit2/ --out-prefix multi_sample
"),
    "prodigal" = paste0("
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
    "metawrap" = paste0("
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate metawrap

  mkdir -p result/metawrap
  mkdir -p result/metawrap/gene_fa
  ~/miniconda3/envs/waste/bin/prodigal -i result/megahit/contigs/${sample}.contigs.fa \\
      -d result/prodigal/gene_fa/${sample}.gene.fa \\
      -o result/prodigal/gene_gff/${sample}.gene.gff \\
      -p meta -f gff

  mkdir -p result/prodigal/fullgene_fa
  grep 'partial=00' result/prodigal/gene_fa/${sample}.gene.fa | cut -f1 -d ' '| sed 's/>//' > result/prodigal/gene_fa/${sample}.fullid
  seqkit grep -f result/prodigal/gene_fa/${sample}.fullid result/prodigal/gene_fa/${sample}.gene.fa > result/prodigal/fullgene_fa/${sample}.gene.fa

"),
    "cluster" = paste0('
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

  mv result/NR/nucleotide_rep_seq.fasta result/NR/nucleotide.fna
  rm result/NR/nucleotide_all_seqs.fasta
  ~/miniconda3/envs/waste/bin/seqkit translate result/NR/nucleotide.fna > result/NR/protein.faa
  sed \'s/\\*//g\' result/NR/protein.faa > result/NR/protein_no_end.faa
  rm result/NR/protein.faa
'),
    "seq_stat" = paste0("
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
  if [ -f result/NR/nucleotide.fna ]; then
    ~/miniconda3/envs/waste/bin/seqkit stats result/NR/nucleotide.fna >result/NR/nucleotide_stat
  fi
"),
    "salmon-index" = paste0("
  ## \u5efa\u7d22\u5f15, -t\u5e8f\u5217, -i \u7d22\u5f15
  # \u5927\u70b9\u5185\u5b58\uff0c30\u500d\uff1f
  mkdir -p result/salmon
  ~/miniconda3/envs/waste/share/salmon/bin/salmon index \\
    -t result/NR/nucleotide.fna \\
    -p 4 \\
    -i result/salmon/index
"),
    "salmon-quant" = paste0("
  mkdir -p result/salmon/quant
  ~/miniconda3/envs/waste/share/salmon/bin/salmon quant \\
      --validateMappings \\
      -i result/salmon/index -l A -p 4 --meta \\
      -1 data/rm_human/${sample}.1.fq \\
      -2 data/rm_human/${sample}.2.fq \\
      -o result/salmon/quant/${sample}.quant
"),
    "salmon-merge" = paste0("
  ls result/salmon/quant/*.quant/quant.sf |awk -F'/' '{print $(NF-1)}' |sed 's/.quant//' >tmp_finished
  diff_rows -f $samplelist -s tmp_finished >tmp_diff
  # \u8ba1\u7b97\u7ed3\u679c\u7684\u884c\u6570
  line_count=$( cat tmp_diff| wc -l)

  # \u68c0\u67e5\u884c\u6570\u662f\u5426\u5927\u4e8e1\uff0c\u5982\u679c\u662f\u5219\u7ec8\u6b62\u811a\u672c
  if [ \"$line_count\" -gt 1 ]; then
    echo 'sf\u6587\u4ef6\u548csamplelist\u6570\u91cf\u4e0d\u4e00\u81f4'
    exit 1
  fi

  ##mapping rates
  cat result/salmon/quant/*.quant/logs/salmon_quant.log |grep 'Mapping rate = '|awk '{print $NF}'> tmp
  paste $samplelist tmp > result/salmon/mapping_rates

  ## combine
  ~/miniconda3/envs/waste/share/salmon/bin/salmon quantmerge \\
      --quants result/salmon/quant/*.quant \\
      -o result/salmon/gene.TPM
  ~/miniconda3/envs/waste/share/salmon/bin/salmon quantmerge \\
      --quants result/salmon/quant/*.quant \\
      --column NumReads -o result/salmon/gene.count
  sed -i '1 s/.quant//g' result/salmon/gene.*
"),
    "eggnog" = paste0("
  ## \u5927\u70b9\u5185\u5b58, \u6570\u636e\u5e93\u670950G\u5de6\u53f3
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate func
  ## diamond\u6bd4\u5bf9\u57fa\u56e0\u81f3eggNOG 5.0\u6570\u636e\u5e93, 1~9h\uff0c\u9ed8\u8ba4diamond 1e-3
  mkdir -p result/eggnog
  emapper.py --no_annot --no_file_comments --override \\
    --data_dir ~/db/eggnog5 \\
    -i result/NR/protein_no_end.faa \\
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
    "cazy" = paste0("
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate func
  mkdir -p result/dbcan2
  diamond blastp   \\
  	--db ~/db/dbcan2/CAZyDB.07312020  \\
  	--query result/NR/protein_no_end.faa \\
  	--threads 8 -e 1e-5 --outfmt 6 \\
  	--max-target-seqs 1 --quiet \\
  	--out result/dbcan2/gene_diamond.f6
"),
    "rgi" = paste0("
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate rgi
  mkdir -p result/card
  ~/miniconda3/envs/rgi/bin/rgi main \\
    --input_sequence result/NR/protein_no_end.faa \\
    --output_file result/card/protein \\
    --input_type protein --num_threads 8 \\
    --clean --alignment_tool DIAMOND # --low_quality #partial genes
"),
    "vfdb" = paste0("
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate func
  mkdir -p result/vfdb
  diamond blastp   \\
  	--db ~/db/VFDB/VFDB_setB_pro \\
  	--query result/NR/protein_no_end.faa \\
  	--threads 8 -e 1e-5 --outfmt 6 \\
  	--max-target-seqs 1 --quiet \\
  	--out result/vfdb/gene_diamond.f6
"),
    "ice" = paste0("
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate func
  mkdir -p result/ice
  diamond blastp   \\
  	--db ~/db/ARG_VF_ICE \\
  	--query result/NR/protein_no_end.faa \\
  	--threads 8 -e 1e-5 --outfmt 6 \\
  	--max-target-seqs 1 --quiet \\
  	--out result/ice/gene_diamond.f6
"),
    "summary" = paste0("
  mkdir -p result/summ_table
  if [ -f result/eggnog/eggnog_anno_output ]; then
  # \u6c47\u603b\uff0c9\u5217KO\uff0c16\u5217CAZy\u6309\u9017\u53f7\u5206\u9694\uff0c21\u5217COG\u6309\u5b57\u6bcd\u5206\u9694\uff0c\u539f\u59cb\u503c\u7d2f\u52a0
  ~/script/others/summarizeAbundance.py \\
    -i result/salmon/gene.count \\
    -m result/eggnog/eggnog_anno_output \\
    -c '9,16,21' -s ',+,+*' -n raw \\
    -o result/summ_table/eggnog
  sed -i 's/^ko://' summ_table/eggnog.KO.raw.txt
  fi

  if [ -f result/card/protein.txt ]; then
  awk 'BEGIN {FS = \"\\t\"; OFS = \"\\t\"} {split($1, a, \" \"); $1 = a[1]} 1' result/card/protein.txt >result/card/protein_format_id.txt
  ~/script/others/summarizeAbundance.py \\
    -i result/salmon/gene.count \\
    -m result/card/protein_format_id.txt \\
    -c '11' -s ';' -n raw \\
    -o result/summ_table/card
  fi

  if [ -f result/dbcan2/gene_diamond.f6 ]; then
  # \u63d0\u53d6\u57fa\u56e0\u4e0edbcan\u5206\u7c7b\u5bf9\u5e94\u8868
  perl ~/script/format_dbcan2list.pl \
    -i result/dbcan2/gene_diamond.f6 \
    -o result/dbcan2/gene.list
  # \u6309\u5bf9\u5e94\u8868\u7d2f\u8ba1\u4e30\u5ea6
  ~/script/others/summarizeAbundance.py \\
    -i result/salmon/gene.count \\
    -m result/dbcan2/gene.list \\
    -c '2' -s ',' -n raw \\
    -o result/dbcan2/cazy
  fi

  if [ -f result/vfdb/gene_diamond.f6 ]; then
  sed -i '1 i Name\\tvf\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore' \\
    result/vfdb/gene_diamond.f6
  ~/script/others/summarizeAbundance.py \\
    -i result/salmon/gene.count \\
    -m result/vfdb/gene_diamond.f6 \\
    -c '2' -s ';' -n raw \\
    -o result/summ_table/vfdb
  fi

  if [ -f result/ice/gene_diamond.f6 ]; then
  sed -i '1 i Name\\tice\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore' \\
    result/ice/gene_diamond.f6
  ~/script/others/summarizeAbundance.py \\
    -i result/salmon/gene.count \\
    -m result/ice/gene_diamond.f6 \\
    -c '2' -s ';' -n raw \\
    -o result/summ_table/ice
  fi
")
  )
}

virome_works <- {
  list(
    "genomad" = paste0("
  # ca rgi
  mkdir -p result/genomad_out
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate rgi
  ~/miniconda3/envs/rgi/bin/genomad \\
    end-to-end result/megahit/contigs/${sample}.contigs.fa \\
    result/genomad_out/${sample}_out ~/db/genomad_db/genomad_db \\
    --disable-nn-classification --cleanup -t 4
"),
    "checkv" = paste0("
  # ca rgi
  mkdir result/checkv_out/
  mkdir result/all_virus/
  cat result/genomad_out/*_out/*.contigs_summary/*.contigs_virus.fna > result/all_virus/all_virus.fna
  cat result/genomad_out/*_out/*.contigs_summary/*.contigs_virus_proteins.faa  > result/all_virus/all_virus_proteins.faa
  cat result/genomad_out/*_out/*.contigs_summary/*.contigs_virus_summary.tsv |grep -v 'seq_name' > result/all_virus/all_contigs_virus_summary.tsv
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate rgi
  checkv end_to_end result/all_virus/all_virus.fna result/checkv_out/ \\
    -t 8 -d ~/db/genomad_db/checkv-db-v1.5 --remove_tmp
"),
    "vcontact2" = paste0("
  # ca vContact2
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate vContact2
  vcontact2 --raw-proteins all_virus_gene.faa \\
    --proteins-fp all_virus_genomes_g2g.csv \\
    --db 'ProkaryoticViralRefSeq211-Merged' \\
    --rel-mode 'Diamond' --pcs-mode MCL \\
    --vcs-mode ClusterONE --threads 8 \\
    --output-dir vConTACT2_Results_211
"),
    "genomad" = paste0("
  # ca rgi
  mkdir -p result/genomad_out
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate rgi
  ~/miniconda3/envs/rgi/bin/genomad \\
    end-to-end result/megahit/contigs/${sample}.contigs.fa \\
    result/genomad_out/${sample}_out ~/db/genomad_db/genomad_db \\
    --disable-nn-classification --cleanup -t 4
")
  )
}

micro_works <- append(micro_works, virome_works)

#' Microbiome sbatch
#'
#' @param work_dir work_dir
#' @param step "fastp","rm_human","megahit","prodigal","salmon-quant",...
#' @param all_sample_num all sample number
#' @param array array number
#' @param partition partition
#' @param mem_per_cpu mem_per_cpu, "2G"
#' @param cpus_per_task cpus_per_task
#' @return No value
#' @export
micro_sbatch <- function(work_dir = "/share/home/jianglab/pengchen/work/asthma/",
                         step = "fastp", all_sample_num = 40, array = 1,
                         partition = "cpu", cpus_per_task = 1, mem_per_cpu = "2G") {
  num_each_array <- ceiling(all_sample_num / array)
  if (!endsWith(work_dir, "/")) work_dir <- paste0(work_dir, "/")
  array <- ifelse(array > 1, paste0("1-", array), "1")

  header <- {
    paste0("#!/bin/bash
#SBATCH --job-name=", step, "
#SBATCH --output=", work_dir, "log/%x_%A_%a.out
#SBATCH --error=", work_dir, "log/%x_%A_%a.err
#SBATCH --array=", array, "
#SBATCH --partition=", partition, "
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=", cpus_per_task, "
#SBATCH --mem-per-cpu=", mem_per_cpu, "
")
  }

  set <- {
    paste0("samplelist=", work_dir, 'samplelist

echo start: `date +"%Y-%m-%d %T"`
start=`date +%s`

####################
')
  }

  set2 <- {
    paste0('echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
START=$SLURM_ARRAY_TASK_ID

NUMLINES=', num_each_array, ' #how many sample in one array

STOP=$((SLURM_ARRAY_TASK_ID*NUMLINES))
START="$(($STOP - $(($NUMLINES - 1))))"

#set the min and max
if [ $START -lt 1 ]
then
  START=1
fi
if [ $STOP -gt ', all_sample_num, " ]
then
  STOP=", all_sample_num, '
fi

echo "START=$START"
echo "STOP=$STOP"

####################
')
  }

  loop1 <- {
    paste0('for (( N = $START; N <= $STOP; N++ ))
do
  sample=$(head -n "$N" $samplelist | tail -n 1)
  echo $sample
  start1=`date +%s`')
  }

  if (step %in% names(micro_works)) {
    work <- {
      micro_works[[step]]
    }
  } else {
    work <- "
  do something"
  }

  loop2 <- "  end1=`date +%s`
  echo `expr $end1 - $start1`s
done
"

  end <- {
    paste0('
##############
echo end: `date +"%Y-%m-%d %T"`
end=`date +%s`
echo TIME:`expr $end - $start`s')
  }

  if (step == "kraken") {
    res_text <- paste0("#!/bin/bash
#SBATCH --job-name=kraken2M
#SBATCH --output=", work_dir, "log/%x_%A_%a.out
#SBATCH --error=", work_dir, "log/%x_%A_%a.err
#SBATCH --time=14-00:00:00
#SBATCH --partition=mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=100G

fqp=", work_dir, "rm_human
python /share/home/jianglab/shared/krakenDB/K2ols/kraken2M.py -t 32 \\
    -i ${fqp} \\
    -c 0.05 \\
    -s _1.fastq,_2.fastq \\
    -o ", work_dir, "result/kraken/ \\
    -d /share/home/jianglab/shared/krakenDB/mydb2 \\
    -k ~/miniconda3/envs/waste/bin/kraken2 \\
    -kt /share/home/jianglab/shared/krakenDB/K2ols/KrakenTools")
  }

  if (step %in% c(
    "fastp", "rm_human", "kraken2", "megahit", "prodigal", "metaspades", "salmon-quant", "humann",
    "genomad"
  )) {
    res_text <- paste0(header, set, set2, loop1, work, loop2, end)
  } else {
    res_text <- paste0(header, set, work, end)
  }

  lib_ps("clipr", library = FALSE)
  clipr::write_clip(res_text, allow_non_interactive = TRUE)
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
pre_fastp <- function(jsonfiles, prefix = c("Raw", "Clean")) {
  lib_ps("jsonlite", library = FALSE)
  result_list <- jsonfiles
  result_dict <- list()
  merge_result_dict <- list()
  if (length(prefix) < 2) prefix <- c("Raw", "Clean")
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

  df <- do.call(rbind, merge_result_dict)
  colnames(df) <- c(
    paste0(prefix[1], "_", c("reads", "bases", "Q20", "Q30", "GC")),
    paste0(prefix[2], "_", c("reads", "bases", "Q20", "Q30", "GC")),
    "Duplication", "Insert_size",
    paste0(prefix[2], "/", prefix[1])
  )
  df
}

#' Prepare the result from assembly_stats (.json file)
#'
#' @param jsonfiles the directory contains .json file
#'
#' @return data.frame
#' @export
#'
pre_assembly_stats <- function(jsonfiles) {
  lib_ps("jsonlite", library = FALSE)
  result_list <- jsonfiles
  result_dict <- list()
  merge_result_dict <- list()
  for (i in result_list) {
    result_dict[[i]] <- jsonlite::fromJSON(i)
  }

  for (k in names(result_dict)) {
    v <- result_dict[[k]]
    key <- strsplit(basename(k), "\\.")[[1]][1]
    v$`Contig Stats`$name <- key
    merge_result_dict[[key]] <- as.data.frame(v$`Contig Stats`)
  }

  df <- do.call(rbind, merge_result_dict)
  df <- cbind(df["name"], df[, -ncol(df)])
  df
}

#' Preprocess GTDB-Tk Classification Results
#'
#' This function reads and processes the output files from a GTDB-Tk `classify` workflow.
#' It combines bacterial (bac120) and archaeal (ar53) classification summaries and phylogenetic trees (if available) into a unified format.
#'
#' @param classify_dir A character string specifying the path to the GTDB-Tk `classify` output directory.
#'                     This directory should contain files like `gtdbtk.bac120.summary.tsv`, `gtdbtk.bac120.classify.tree`, etc.
#'
#' @return A list with two components:
#' \describe{
#'   \item{gtdb_res}{A data frame containing the combined taxonomic classification for all genomes.
#'                   The `classification` column is parsed into standard taxonomic ranks (Domain to Species).}
#'   \item{tree}{A phylogenetic tree (phylo object) combining the bacterial and (if present) archaeal trees.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if the provided directory exists and contains the necessary `*.summary.tsv` files.
#'   \item Reads the bacterial backbone tree.
#'   \item If an archaeal tree file exists, it binds it to the bacterial tree.
#'   \item Reads and combines all `*.summary.tsv` files in the directory.
#'   \item Parses the semicolon-separated `classification` string into separate columns for each taxonomic rank.
#'   \item Ensures the resulting taxonomy table has standard ranks (Domain, Phylum, Class, Order, Family, Genus, Species).
#' }
#' @export
pre_gtdb_tk <- function(classify_dir) {
  lib_ps("ape", library = FALSE)
  # 1. Check if the directory exists and contains summary files
  if (!dir.exists(classify_dir)) {
    stop("The provided directory '", classify_dir, "' does not exist.")
  }

  summary_files <- list.files(classify_dir, pattern = ".*\\.summary\\.tsv$", full.names = TRUE)
  if (length(summary_files) == 0) {
    stop("No files matching the pattern '*.summary.tsv' were found in the directory '", classify_dir, "'.")
  }

  # 2. Read and combine the phylogenetic trees
  # Read bacterial tree
  bac_tree_file <- file.path(classify_dir, "gtdbtk.backbone.bac120.classify.tree")
  if (!file.exists(bac_tree_file)) {
    stop("The bacterial tree file 'gtdbtk.backbone.bac120.classify.tree' was not found in '", classify_dir, "'.")
  }
  tree1 <- ape::read.tree(bac_tree_file)

  # Check and read archaeal tree if it exists, then bind it
  ar_tree_file <- file.path(classify_dir, "gtdbtk.ar53.classify.tree")
  if (file.exists(ar_tree_file)) {
    tree2 <- ape::read.tree(ar_tree_file)
    tree1 <- ape::bind.tree(tree1, tree2)
  }

  # 3. Read and combine all summary TSV files
  gtdb_list <- lapply(summary_files, utils::read.table, sep = "\t", header = TRUE)
  gtdb_res <- do.call(rbind, gtdb_list)

  # 4. Parse the classification column
  # Split the classification string by semicolon
  gtdb_tax <- pcutils::strsplit2(gtdb_res$classification, ";", colnames = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

  # 5. (Optional) Apply further taxonomy table preprocessing
  gtdb_tax <- pre_tax_table(gtdb_tax, tax_levels = c("d", "p", "c", "o", "f", "g", "s", "st"), add_prefix = FALSE)

  # 6. Combine user genome ID with the parsed taxonomy table
  # Keep the user_genome column and bind it with the taxonomy data
  final_result <- cbind(gtdb_res["user_genome"], gtdb_tax)

  tree1 <- ape::drop.tip(tree1, setdiff(tree1$tip.label, final_result$user_genome))

  # 7. Return the results as a list
  result_list <- list(
    gtdb_res = final_result,
    tree = tree1
  )

  return(result_list)
}

#' Plot GTDB-Tk Phylogenetic Tree with Taxonomic Coloring
#'
#' This function creates a circular phylogenetic tree visualization from GTDB-Tk results,
#' with tips colored by specified taxonomic level. It can optionally represent only one genome per species.
#'
#' @param gtdb_res A data frame containing GTDB-Tk classification results, typically
#'                 the output from \code{\link{pre_gtdb_tk}} function.
#' @param tree A phylogenetic tree object (phylo class) containing all genomes.
#' @param tax_level Character specifying the taxonomic level for coloring tips.
#'                  Default is "Phylum". Other options include "Class", "Order", "Family", "Genus", "Species".
#' @param represented Logical indicating whether to include only one representative genome
#'                    per species. Default is TRUE.
#' @param layout Character specifying the tree layout. Options include "fan", "circular",
#'               "rectangular", etc. Default is "fan".
#' @param branch_size Numeric value for branch line width. Default is 0.2.
#' @param color_na Character specifying color for tips with missing taxonomic information.
#'                 Default is "black".
#'
#' @return A ggplot object containing the phylogenetic tree visualization.
#'
#' @export
plot_gtdb_tr <- function(gtdb_res, tree, tax_level = "Phylum", represented = TRUE,
                         layout = "fan", branch_size = 0.2, color_na = "black") {
  lib_ps("ggtree", "ggtreeExtra", library = FALSE)
  requireNamespace("ggplot2")
  Species <- `_color` <- NULL
  # Input validation
  if (!inherits(tree, "phylo")) {
    stop("The 'tree' parameter must be a phylogenetic tree object of class 'phylo'")
  }
  if (!is.data.frame(gtdb_res)) {
    stop("The 'gtdb_res' parameter must be a data frame")
  }
  if (!"user_genome" %in% names(gtdb_res)) {
    stop("The 'gtdb_res' parameter must contain a 'user_genome' column with genome identifiers")
  }
  if (!tax_level %in% names(gtdb_res)) {
    stop("Taxonomic level '", tax_level, "' not found in gtdb_res data frame")
  }

  # Filter genomes: either all genomes or one representative per species
  if (represented) {
    genome_info_f <- dplyr::distinct(gtdb_res, Species, .keep_all = TRUE)
  } else {
    genome_info_f <- gtdb_res
  }

  # Prune the tree to include only selected genomes
  tips_to_keep <- genome_info_f$user_genome
  tips_in_tree <- tips_to_keep[tips_to_keep %in% tree$tip.label]

  if (length(tips_in_tree) == 0) {
    stop("No genome names in 'gtdb_res$user_genome' match the tip labels in the phylogenetic tree")
  }

  our_tr <- ape::drop.tip(tree, setdiff(tree$tip.label, tips_in_tree))

  # Prepare tree data for plotting
  our_tr_df <- fortify(our_tr)

  # Merge with taxonomic information
  tax_info <- genome_info_f[, c("user_genome", tax_level)]
  our_tr_df <- dplyr::left_join(our_tr_df, tax_info, by = c("label" = "user_genome"))

  # Clean taxonomic names by removing prefix (e.g., "p__", "c__", etc.)
  prefix <- tolower(substr(tax_level, 1, 1))
  our_tr_df[[tax_level]] <- gsub(paste0(prefix, "__"), "", our_tr_df[[tax_level]])
  our_tr_df[["_color"]] <- our_tr_df[[tax_level]]

  # Generate color palette using your existing get_cols function
  legends <- get_cols(our_tr_df[["_color"]])

  # Create the tree plot
  p <- ggtree::ggtree(our_tr_df,
    layout = layout, mapping = aes(color = `_color`),
    linewidth = branch_size
  ) +
    scale_color_manual(
      values = legends, name = tax_level, na.value = color_na,
      label = c(names(legends), "")
    ) +
    guides(colour = guide_legend(
      order = 1,
      override.aes = list(
        linewidth = 1,
        linetype = setNames(c(rep(1, length(legends)), NA), legends)
      )
    ))

  return(p)
}

#' Visualize CheckM2 Genome Quality Assessment Results
#'
#' This function creates a scatter plot showing genome completeness vs contamination
#' with optional marginal density plots to display distributions.
#'
#' @param checkm2_df Data frame. CheckM2 results containing at least columns:
#'   'Completeness', 'Contamination', and 'Name'.
#' @param add_marginal Logical. Whether to add marginal density plots using ggExtra.
#'   Default is TRUE.
#' @param marginal_type Character. Type of marginal plot: "density", "histogram",
#'   "boxplot", or "violin". Default is "density".
#' @param point_size Numeric. Size of points in scatter plot. Default is 0.6.
#' @param base_size Numeric. Base font size for the plot. Default is 14.
#' @param quality_thresholds List. Custom thresholds for quality classification.
#'   Default list(high_comp = 90, high_contam = 5, med_comp = 70, med_contam = 10)
#' @param filter_data Logical. Whether to filter low-quality genomes. Default is TRUE.
#' @param min_quality_score Numeric. Minimum quality score for filtering. Default is 50.
#' @param min_completeness Numeric. Minimum completeness for filtering. Default is 50.
#' @param max_contamination Numeric. Maximum contamination for filtering. Default is 10.
#'
#' @return A ggplot object or ggExtra plot object if marginal plots are added.
#'
#' @export
plot_checkm2_res <- function(checkm2_df,
                             add_marginal = TRUE,
                             marginal_type = "density",
                             point_size = 0.6,
                             base_size = 14,
                             quality_thresholds = list(
                               high_comp = 90, high_contam = 5,
                               med_comp = 70, med_contam = 10
                             ),
                             filter_data = TRUE,
                             min_quality_score = 50,
                             min_completeness = 50,
                             max_contamination = 10) {
  Completeness <- Contamination <- Group <- Quality_score <- Group1 <- Group2 <- SampleID <- ratio <- NULL
  lib_ps("ggExtra", library = FALSE)
  # Input validation
  required_cols <- c("Completeness", "Contamination", "Name")
  missing_cols <- setdiff(required_cols, colnames(checkm2_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Create a copy of the data
  genome_info <- checkm2_df

  # Quality classification
  genome_info <- genome_info %>%
    mutate(Group = case_when(
      Completeness >= quality_thresholds$high_comp & Contamination <= quality_thresholds$high_contam ~ "High quality",
      Completeness >= quality_thresholds$med_comp & Contamination <= quality_thresholds$med_contam ~ "Medium quality",
      TRUE ~ "Partial"
    ))

  # Calculate quality score
  genome_info <- genome_info %>%
    mutate(Quality_score = Completeness - 5 * Contamination)

  # Filter data if requested
  if (filter_data) {
    genome_info <- genome_info %>%
      filter(
        Quality_score >= min_quality_score,
        Completeness >= min_completeness,
        Contamination <= max_contamination
      )

    if (nrow(genome_info) == 0) {
      warning("No genomes passed the quality filters. Showing all data.")
      genome_info <- checkm2_df %>%
        mutate(
          Group = "All",
          Quality_score = Completeness - 5 * Contamination
        )
    }
  }

  # Calculate group statistics for labels
  genome_count <- genome_info %>%
    count(Group) %>%
    mutate(
      ratio = round(100 * (n / sum(n)), 1),
      Group1 = paste0(Group, " (", n, ", ", ratio, "%)")
    )

  # Create mapping for group labels
  tmp <- setNames(genome_count$Group1, genome_count$Group)
  genome_info$Group2 <- tmp[genome_info$Group]

  # Ensure factor levels are in logical order
  quality_levels <- c("High quality", "Medium quality", "Partial", "All")
  existing_levels <- intersect(quality_levels, names(tmp))
  genome_info$Group2 <- factor(genome_info$Group2,
    levels = unname(tmp[existing_levels])
  )

  # Extract SampleID from Name (assuming format "SampleID_rest")
  genome_info$SampleID <- gsub("_.*", "", genome_info$Name)

  # Create base scatter plot
  p <- ggplot(genome_info) +
    geom_point(aes(x = Completeness, y = Contamination, col = Group2),
      size = point_size, alpha = 0.7
    ) +
    scale_color_manual(
      values = c("#78c679", "#80b1d3", "#fdb462", "#ff6b6b"),
      name = NULL
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = c(0.33, 0.8),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
      axis.text = element_text(colour = "black", size = base_size - 1),
      panel.grid.major = element_line(linewidth = 0.3),
      panel.grid.minor = element_line(linewidth = 0.1)
    ) +
    labs(
      title = "CheckM2 Genome Quality Assessment",
      subtitle = paste("Total genomes:", nrow(genome_info)),
      x = "Completeness (%)",
      y = "Contamination (%)"
    )

  # Add marginal plots if requested
  if (add_marginal) {
    if (!requireNamespace("ggExtra", quietly = TRUE)) {
      warning("ggExtra package not installed. Proceeding without marginal plots.")
      return(p)
    }

    p <- ggExtra::ggMarginal(
      p = p,
      type = marginal_type,
      margins = "both",
      size = 5,
      colour = "black",
      groupFill = TRUE,
      alpha = 0.3
    )
  }

  return(p)
}


#' Preprocess MPA Species Abundance and Taxonomy Data
#'
#' This function reads a species abundance profile from kraken2 format_report output,
#' performs filtering to remove low-abundance and unwanted taxa, and extracts
#' a clean taxonomy table for downstream analysis.
#'
#' @param dir Character. Path to the directory containing the input file
#'   "mpa_profile_species.txt". Default is an empty string (current directory).
#' @param exclude Character. Pattern to exclude specific taxa (e.g., "g__Streptococcus").
#'   Uses grepl() for pattern matching. Default is NULL (no exclusion).
#' @param relative_threshold Numeric. Relative abundance threshold for filtering
#'   low-abundance taxa. Taxa with mean relative abundance below this threshold
#'   will be removed. Default is 1e-4.
#'
#' @return A list with two components:
#'   \item{species}{A matrix of filtered species abundance data}
#'   \item{species_taxonomy}{A data.frame of taxonomy information for each species}
#'
#' @details
#' The function performs the following steps:
#' 1. Reads the Metaphlan species profile table
#' 2. Removes "unclassified" entries and "cellular_organisms" category
#' 3. Filters out taxa matching the \code{exclude} pattern (if provided)
#' 4. Applies relative abundance filtering using \code{pcutils::rm_low}
#' 5. Extracts and formats taxonomy information from the Metaphlan-style names
#' 6. Cleans species names by removing the "s__" prefix
#'
#' @importFrom readr read_table
#' @importFrom tibble column_to_rownames
#' @importFrom pcutils rm_low strsplit2
#' @export
pre_format_report <- function(dir = "", exclude = NULL, relative_threshold = 1e-4) {
  # Construct file path
  file_path <- file.path(dir, "mpa_profile_species.txt")

  # Read and preprocess abundance data
  species_abundance <- readr::read_table(file_path, show_col_types = FALSE) %>%
    tibble::column_to_rownames("clade")

  # Remove unclassified entries and cellular organisms
  species_abundance <- species_abundance[
    rownames(species_abundance) != "unclassified" &
      !grepl("cellular_organisms", rownames(species_abundance)),
  ]

  # Apply exclusion filter if specified
  if (!is.null(exclude)) {
    species_abundance <- species_abundance[
      !grepl(exclude, rownames(species_abundance)),
    ]
  }

  # Filter low abundance taxa
  species_filtered <- pcutils::rm_low(species_abundance,
    relative_threshold = relative_threshold
  )

  # Extract taxonomy information
  taxonomy_parts <- pcutils::strsplit2(rownames(species_filtered), "\\|") %>%
    as.data.frame()

  colnames(taxonomy_parts) <- c(
    "Domain", "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus", "Species"
  )

  # Handle special case for Viruses
  taxonomy_parts[taxonomy_parts$Domain == "d__Viruses", "Kingdom"] <- "k__Viruses"

  # Clean species names and set as rownames
  clean_species_names <- gsub("s__", "", taxonomy_parts$Species)
  rownames(species_filtered) <- clean_species_names
  rownames(taxonomy_parts) <- clean_species_names

  # Return results as a named list
  list(
    species = species_filtered,
    species_taxonomy = taxonomy_parts
  )
}
#' Preprocess geNomad and CheckV Output Results
#'
#' This function automatically processes geNomad output files by detecting
#' sample names from the directory structure and optionally integrates CheckV
#' quality assessment results.
#'
#' @param genomad_out_dir Character. Path to the geNomad output directory.
#'   This directory should contain sample-specific subdirectories with the
#'   pattern "*.contigs_summary".
#' @param checkV_out_dir Character. Optional path to the CheckV output directory.
#'   If provided, quality summary will be integrated. Default is NULL.
#' @param provirus Logical. Whether to identify and separate provirus sequences.
#'   Default is TRUE.
#' @param filter Logical. Whether to apply quality filtering to viral sequences.
#'   Default is TRUE.
#' @param min_length Numeric. Minimum sequence length for filtering. Default is 1000.
#' @param min_completeness Numeric. Minimum completeness score for CheckV filtering.
#'   Default is 50.
#' @param checkV_out_prefix Character. Optional prefix to remove from CheckV contig IDs.
#'
#' @return An object of class "virus_res" containing four components:
#'   \item{sample}{Detected sample name}
#'   \item{virus_summary}{Integrated data frame with geNomad and optional CheckV results}
#'   \item{virus_genes}{Gene-level annotations from geNomad}
#'   \item{valid_virus}{Filtered high-quality viral sequences}
#'
#' @details
#' The function automatically detects sample names by searching for directories
#' with the pattern "*.contigs_summary" within the genomad_out_dir. It then
#' extracts the sample name by removing the ".contigs_summary" suffix.
#'
#' @examples
#' \dontrun{
#' # Basic usage - sample name will be automatically detected
#' virus_results <- pre_genomad(genomad_out_dir = "~/Documents/R/Lung_virome/data/genomad_out2/")
#'
#' # Access the detected sample name
#' sample_name <- virus_results$sample
#' print(paste("Detected sample:", sample_name))
#' }
#'
#' @importFrom dplyr filter pull left_join
#' @importFrom readr read_delim
#' @export
pre_genomad <- function(genomad_out_dir = "", checkV_out_dir = NULL,
                        provirus = TRUE, filter = TRUE, checkV_out_prefix = NULL, min_length = 1000,
                        min_completeness = 50) {
  topology <- virus_score <- n_hallmarks <- marker_enrichment <- contig_length <- completeness <- contig_id <- NULL
  # Automatically detect sample name from directory structure
  contigs_dirs <- list.files(genomad_out_dir,
    pattern = "\\.contigs_summary$",
    full.names = FALSE
  )

  if (length(contigs_dirs) == 0) {
    stop("No '.contigs_summary' directories found in: ", genomad_out_dir)
  }

  # Extract sample name by removing the .contigs_summary suffix
  sample <- gsub("\\.contigs_summary$", "", contigs_dirs[1])

  if (length(contigs_dirs) > 1) {
    warning("Multiple '.contigs_summary' directories found. Using first one: ", sample)
  }

  cat("Detected sample:", sample, "\n")

  # Construct geNomad file paths using detected sample name
  genomad_summary_file <- file.path(
    genomad_out_dir,
    paste0(sample, ".contigs_summary"),
    paste0(sample, ".contigs_virus_summary.tsv")
  )

  genomad_genes_file <- file.path(
    genomad_out_dir,
    paste0(sample, ".contigs_summary"),
    paste0(sample, ".contigs_virus_genes.tsv")
  )

  # Check if geNomad files exist
  if (!file.exists(genomad_summary_file)) {
    stop("geNomad summary file not found: ", genomad_summary_file)
  }
  if (!file.exists(genomad_genes_file)) {
    stop("geNomad genes file not found: ", genomad_genes_file)
  }

  # Read geNomad virus summary
  virus_summary <- readr::read_delim(genomad_summary_file, show_col_types = FALSE)
  colnames(virus_summary)[1] <- "contig_id"

  # Read geNomad virus genes
  virus_genes <- readr::read_delim(genomad_genes_file, show_col_types = FALSE)
  virus_genes$contig_id <- gsub("_\\d+$", "", virus_genes$gene)

  # Integrate CheckV results if provided
  if (!is.null(checkV_out_dir)) {
    checkv_summary_file <- file.path(checkV_out_dir, "quality_summary.tsv")

    if (file.exists(checkv_summary_file)) {
      checkV <- readr::read_delim(checkv_summary_file, show_col_types = FALSE)
      if (!is.null(checkV_out_prefix)) checkV$contig_id <- gsub(checkV_out_prefix, "", checkV$contig_id)
      virus_summary <- dplyr::left_join(virus_summary, checkV, by = "contig_id")
      cat("CheckV results successfully integrated\n")
    } else {
      warning("CheckV summary file not found: ", checkv_summary_file)
    }
  }

  # Handle provirus sequences
  if (provirus && "topology" %in% names(virus_summary)) {
    Provirus <- dplyr::filter(virus_summary, topology == "Provirus")
    cat("Found", nrow(Provirus), "provirus sequences\n")
    virus_summary <- dplyr::filter(virus_summary, topology != "Provirus")
  }

  # Apply quality filtering
  valid_virus <- data.frame()
  if (filter) {
    valid_virus <- rbind(
      dplyr::filter(
        virus_summary,
        (length > 5000 & length <= 10000) &
          virus_score > 0.9 &
          n_hallmarks > 0 &
          marker_enrichment > 2
      ),
      dplyr::filter(
        virus_summary,
        length > 10000 &
          ((virus_score > 0.8 & n_hallmarks > 0) | marker_enrichment > 5)
      )
    )

    # Additional CheckV filtering if available
    if (!is.null(checkV_out_dir) && exists("checkV")) {
      checkV_high_quality <- dplyr::filter(
        checkV,
        contig_length > min_length &
          completeness > min_completeness
      )
      valid_virus <- dplyr::filter(valid_virus, contig_id %in% checkV_high_quality$contig_id)
    }

    message("After quality filtering,", nrow(valid_virus), " high-quality viral sequences retained\n")
  } else {
    valid_virus <- virus_summary
  }

  # Create results object with sample name
  virus_res <- list(
    sample = sample,
    virus_summary = virus_summary,
    virus_genes = virus_genes,
    valid_virus = valid_virus
  )

  class(virus_res) <- "virus_res"
  return(virus_res)
}

#' Print method for virus_res objects
#'
#' @param x An object of class virus_res
#' @param ... Additional arguments (not used)
#'
#' @export
print.virus_res <- function(x, ...) {
  cat("Virus Analysis Results (virus_res object)\n")
  cat("=======================================\n")
  cat("Sample:", x$sample, "\n")
  cat("Total viral sequences:", nrow(x$virus_summary), "\n")
  cat("High-quality sequences:", nrow(x$valid_virus), "\n")
  cat("Genes annotated:", nrow(x$virus_genes), "\n")
}

#' Plot Individual Phage Genome Structure with Annotations
#'
#' This function creates a circular genome map for a single phage, displaying gene annotations,
#' functional categories, and other genomic features. It automatically extracts topology information
#' and provides flexible parameter customization.
#'
#' @param one_phage Character. Phage sequence identifier (e.g., "k141_10408").
#' @param genomad_out_res List. Output object from pre_genomad function containing virus summary and gene data.
#' @param anno_table Data frame. Optional annotation table with gene information for left join.
#'   Must contain a 'gene' column for merging with genomad gene data.
#' @param y_var Character. Variable for y-axis mapping. Default is "strand".
#' @param fill_var Character. Variable for fill color mapping. Default is "COG".
#' @param label_wrap Integer. Width for wrapping gene description labels. Default is 20.
#' @param palette Character. Color palette for fill categories. Default is "Set3".
#' @param label_var Character vector. Columns to use for gene description labels.
#'
#' @return A ggplot object displaying the circular phage genome map.
#'
#' @details
#' The function automatically extracts topology information from the virus summary data
#' and creates a polar coordinate visualization showing:
#' - Gene arrows indicating direction and position
#' - Flexible y-axis and fill color mappings
#' - Genome length and automatically detected topology information
#' - Gene descriptions with intelligent label placement
#'
#' @export
plot_one_phage <- function(one_phage, genomad_out_res, anno_table = NULL,
                           y_var = "strand", fill_var = "COG",
                           label_var = c("annotation_description"),
                           label_wrap = 20, palette = "Set3") {
  lib_ps("gggenes", "ggrepel", "stringr", library = FALSE)
  bp <- contig_id <- start <- end <- NULL
  # Input validation
  if (!one_phage %in% genomad_out_res$virus_summary$contig_id) {
    stop("Phage ID '", one_phage, "' not found in virus summary data")
  }

  # Extract gene data for the specific phage
  tmp_gene <- dplyr::filter(genomad_out_res$virus_genes, contig_id == one_phage)

  if (nrow(tmp_gene) == 0) {
    stop("No gene data found for phage: ", one_phage)
  }

  # Extract topology information automatically
  phage_info <- dplyr::filter(genomad_out_res$virus_summary, contig_id == one_phage)
  topology <- if ("topology" %in% colnames(phage_info)) {
    as.character(phage_info$topology[1])
  } else {
    "Unknown"
  }
  if ("taxonomy" %in% colnames(phage_info)) {
    tmp <- strsplit(phage_info$taxonomy, ";")[[1]]
    taxonomy <- tmp[max(which(tmp != ""))]
  } else {
    taxonomy <- "Unknown"
  }

  # Merge with annotation table if provided
  if (!is.null(anno_table)) {
    # Check if annotation table has the required 'gene' column
    if (!"gene" %in% colnames(anno_table)) {
      warning("Annotation table must contain a 'gene' column for merging. Proceeding without annotations.")
    } else {
      # Merge annotation data using the gene column
      tmp_gene <- dplyr::left_join(tmp_gene, anno_table, by = "gene")
    }
  }

  # Create a description for labeling, 优先取最前面非NA的部分
  if (all(label_var %in% colnames(tmp_gene))) {
    tmp_gene$desc <- apply(tmp_gene[, label_var], 1, function(x) {
      desc <- x[!is.na(x) & x != ""]
      if (length(desc) > 0) {
        return(desc[1])
      } else {
        return(NA)
      }
    })
  } else {
    tmp_gene$desc <- NA
  }

  # Ensure the fill variable exists, create dummy if not
  if (!fill_var %in% colnames(tmp_gene)) {
    tmp_gene[[fill_var]] <- "Unknown"
    warning("Fill variable '", fill_var, "' not found in data. Using default 'Unknown' values.")
  }

  # Ensure the y-axis variable exists
  if (!y_var %in% colnames(tmp_gene)) {
    # If strand information is not available, set all to 1 for single strand display
    tmp_gene[[y_var]] <- 1
    message("Y-axis variable '", y_var, "' not found. Displaying as single strand.")
  }

  # Calculate genome length and adjust x-axis limits based on topology
  bp <- max(c(tmp_gene$start, tmp_gene$end), na.rm = TRUE)

  # Adjust x-axis limits based on topology type
  if (grepl("DTR", topology, ignore.case = TRUE)) {
    x_end <- bp
  } else if (grepl("ITR", topology, ignore.case = TRUE)) {
    x_end <- bp * 1.05
  } else {
    x_end <- bp * 1.1 # Extra space for other topologies
  }

  # Create dynamic aesthetic mappings
  aes_mapping <- aes(
    xmin = start, xmax = end,
    y = .data[[y_var]],
    fill = .data[[fill_var]]
  )

  # Create the base plot
  p <- ggplot(tmp_gene, aes_mapping) +
    # Genome backbone
    geom_segment(aes(x = 0, xend = bp, y = .data[[y_var]], yend = .data[[y_var]]),
      linewidth = 0.5, color = "gray50"
    ) +
    # Gene arrows with dynamic forward direction
    gggenes::geom_gene_arrow(aes(forward = (.data[[y_var]] > 0)),
      arrowhead_height = unit(3, "mm"),
      arrowhead_width = unit(2, "mm"),
      arrow_body_height = unit(2, "mm")
    ) +
    # Convert to circular genome map
    coord_polar() +
    ylim(c(-10, max(tmp_gene[[y_var]], na.rm = TRUE) + 2)) +
    xlim(c(0, x_end)) +
    theme_void() +
    # Color scheme
    scale_fill_brewer(type = "qual", palette = palette, na.value = "grey90") +
    theme(legend.position = "right") +
    # Genome information annotation
    annotate(
      geom = "text", x = 0, y = -10,
      label = stringr::str_glue("{one_phage}\n{bp} bp; {topology}\n{taxonomy}"),
      size = 4, hjust = 0.5
    )

  # Add gene labels if there are not too many genes
  if (nrow(tmp_gene) <= 50) {
    p <- p + ggrepel::geom_label_repel(
      aes(
        x = (start + end) / 2, y = .data[[y_var]],
        label = stringr::str_wrap(desc, label_wrap)
      ),
      box.padding = 0.3, size = 2.5, max.overlaps = 15,
      segment.curvature = 0.01, label.r = 0
    )
  } else {
    message("Too many genes (", nrow(tmp_gene), ") for labels. Labels omitted for clarity.")
  }

  return(p)
}

#' Visualize Contigs Quality Metrics
#'
#' This function creates a scatter plot to visualize the quality metrics of contigs
#' from geNomad and CheckV analysis results.
#'
#' @param genomad_out_res List. Output object from pre_genomad function containing virus summary data.
#' @param show_valid Logical. Whether to display only valid (high-quality) viral sequences.
#'   Default is FALSE (show all sequences).
#' @param point_size Numeric. Size of the points in the plot. Default is 1.4.
#' @param point_alpha Numeric. Transparency of the points (0-1). Default is 0.6.
#' @param base_size Numeric. Base font size for the plot. Default is 15.
#'
#' @return A ggplot object displaying contigs quality metrics.
#'
#' @export
plot_contigs_quality <- function(genomad_out_res, show_valid = FALSE,
                                 point_size = 1.4, point_alpha = 0.6,
                                 base_size = 15) {
  completeness <- checkv_quality <- topology <- NULL
  # Input validation
  if (!inherits(genomad_out_res, "virus_res")) {
    stop("Input must be a virus_res object from pre_genomad function")
  }

  # Select data based on show_valid parameter
  if (show_valid) {
    plot_data <- genomad_out_res$valid_virus
    if (nrow(plot_data) == 0) {
      warning("No valid viral sequences found. Showing all sequences instead.")
      plot_data <- genomad_out_res$virus_summary
    }
  } else {
    plot_data <- genomad_out_res$virus_summary
  }

  if (nrow(plot_data) == 0) {
    stop("No data available for plotting")
  }
  # Create the plot
  p <- ggplot(plot_data) +
    geom_point(
      aes(
        x = length / 1000,
        y = completeness,
        color = checkv_quality,
        shape = topology
      ),
      size = point_size,
      alpha = point_alpha
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 4)),
      shape = guide_legend(override.aes = list(size = 4))
    ) +
    labs(
      x = "Contig length (kb)",
      y = "Completeness (%)",
      color = "CheckV quality",
      shape = "Topology"
    ) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black", size = base_size - 2)
    )

  return(p)
}
