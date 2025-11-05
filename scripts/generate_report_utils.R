## ---- I/O ----

read_fastp_duplication_batch <- function(fastp_path, n_cpu){
  read_fastp_duplication <- function(path) {
    j <- fromJSON(path, simplifyVector = TRUE)
    h <- j$duplication$rate
    tibble(duplication_rate = h)
  }
  
  files <- list.files(fastp_path, "*.json", full.names = TRUE)
  pbmclapply(files,
             function(f){
               s <- f %>% basename() %>% str_remove("\\.fastp\\.json")
               read_fastp_duplication(f) %>% 
                 mutate(sample = s)
             },
             mc.cores = n_cpu) %>% 
    list_rbind()
}

read_fastp_insert_size_batch <- function(fastp_path, n_cpu){
  read_fastp_insert_size <- function(path) {
    j <- fromJSON(path, simplifyVector = TRUE)
    h <- j$insert_size$histogram
    tibble(insert_size = seq_along(h) - 1L, count = as.integer(h))
  }
  
  files <- list.files(fastp_path, "*.json", full.names = TRUE)
  pbmclapply(files,
             function(f){
               s <- f %>% basename() %>% str_remove("\\.fastp\\.json")
               read_fastp_insert_size(f) %>% 
                 mutate(sample = s)
             },
             mc.cores = n_cpu) %>% 
    list_rbind()
}

read_fastp_filtered_counts_batch <- function(fastp_path, n_cpu){
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  
  read_fastp_filtered_counts <- function(path) {
    j <- fromJSON(path, simplifyVector = TRUE)
    af <- j$summary$after_filtering %||% j$reads$after_filtering %||% j$after_filtering
    if (is.null(af)) stop("after_filtering section not found in fastp JSON")
    tibble(
      total_reads = as.numeric(af[["total_reads"]] %||% af[["reads"]]),
      total_bases = as.numeric(af[["total_bases"]] %||% af[["bases"]])
    )
  }
  files <- list.files(fastp_path, "*.json", full.names = TRUE)
  pbmclapply(files,
             function(f){
               s <- f %>% basename() %>% str_remove("\\.fastp\\.json")
               read_fastp_filtered_counts(f) %>% 
                 mutate(sample = s)
             },
             mc.cores = n_cpu) %>% 
    list_rbind()
}

read_qualimap_depth_across_reference_batch <- function(qualimap_path, n_cpu){
  read_qualimap_depth_across_reference <- function(path){
    read_tsv(path, col_names = c("pos",	"depth", "depth_std"), comment = "#")
  }
  get_sample_name_qualimap <- function(path, qualimap_path){
    splitted <- path %>% 
      str_remove(qualimap_path) %>% 
      str_split_1(pattern = "/")
    splitted <- splitted[splitted != ""] 
    splitted[1]
  }
  
  files <- list.files(qualimap_path, "coverage_across_reference.txt", 
                      full.names = TRUE, recursive = TRUE)
  pbmclapply(files,
             function(f){
               s <- f  %>% get_sample_name_qualimap(qualimap_path)
               read_qualimap_depth_across_reference(f) %>% 
                 mutate(sample = s)
             },
             mc.cores = n_cpu) %>% 
    list_rbind()
}

read_bedtools_coverage_batch <- function(coverage_path, n_cpu){
  read_bedtools_coverage <- function(path){
    read_tsv(path, 
             col_names = c(
               "seqname", "start", "end",
               "n_reads", "covered_bases",
               "window_size", "coverage"))
  }
  
  files <- list.files(coverage_path, "*_coverage.txt", 
                      full.names = TRUE)
  pbmclapply(files,
             function(f){
               s <- f %>% basename() %>% str_remove("_coverage\\.txt")
               read_bedtools_coverage(f) %>% 
                 mutate(sample = s)
             },
             mc.cores = n_cpu) %>% 
    list_rbind()
}

## ---- Helpers ----

get_q_from_hist <- function(p, x, y){
  y_cum <- cumsum(y)/sum(y)
  idx <- which(y_cum <=p) %>% max()
  x[idx]
}

get_bimodal_means <- function(x){
  x <- x[!is.na(x)]
  x <- x[is.finite(x)]
  x <- x[x >= quantile(x, 0.01) & x <= quantile(x, 0.99)]
  km <- kmeans(x, centers = 2, nstart = 20, iter.max = 100)
  means <- tapply(x, km$cluster, mean)
  sds   <- tapply(x, km$cluster, sd)
  ns    <- as.integer(table(km$cluster))
  ord   <- order(as.numeric(means))  # low -> high
  
  tibble(
    mode       = c("low", "high"),
    mean       =  as.numeric(means)[ord],
    sd         = as.numeric(sds)[ord],
    n          = ns[ord],
    prop       = ns[ord] / sum(ns)
  )
}

get_unimodal_means <- function(x){
  x <- x[!is.na(x)]
  x <- x[is.finite(x)]
  x <- x[x >= quantile(x, 0.01) & x <= quantile(x, 0.99)]
  
  tibble(
    mode       = "mean",
    mean       = mean(x),
    sd         = sd(x),
    n          = length(x),
    prop       = NA
  )
}
