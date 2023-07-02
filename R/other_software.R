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

