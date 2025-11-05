## Files needed
## data_dir/proj/qc/*_qc_files/fastp/*sample*.json
## data_dir/proj/qc/*_qc_files/qualimap/*sample*/raw_data_qualimapReport/coverage_across_reference.txt
## data_dir/proj/qc/*_qc_files/coverage/*sample*_coverage.txt
## data_dir/proj/qc/*_qc_files/*_count_sites.tsv

rm(list = ls())

library(jsonlite)
library(pbmcapply)
library(ggpubr)
library(tidyverse)

proj <- "dug_test_against_Fcp"
n_cpu <- 4
bimodal <- TRUE # estimate mean depth and coverage using bimodal distribution

data_dir <- "mapping_qc/data/"
source("mapping_qc/generate_report_utils.R")

## ---- # reads, # bases (fastp) ----
filtered_counts <- read_fastp_filtered_counts_batch(
  file.path(data_dir, proj, "qc", paste0(proj, "_qc_files"), "fastp"),
  n_cpu)

max_reads_M     <- max(filtered_counts$total_reads/1e6,  na.rm = TRUE)
max_bases_100M  <- max(filtered_counts$total_bases/1e8,  na.rm = TRUE)
scale_factor    <- ifelse(is.finite(max_reads_M / max_bases_100M) && max_bases_100M > 0,
                          max_reads_M / max_bases_100M, 1)

filtered_counts_plot_df <- filtered_counts %>%
  mutate(reads_M      = total_reads / 1e6,
         bases_100M   = total_bases / 1e8,
         bases_scaled = bases_100M * scale_factor)

p1 <- filtered_counts_plot_df %>% 
  ggplot(aes(x = sample)) +
    geom_col(aes(y = reads_M)) +
    geom_point(aes(y = bases_scaled), 
               size = 2, color = "red") +
    coord_flip() +
    xlab("Sample") +
    scale_y_continuous(
      name = "Total Reads (M)",
      sec.axis = sec_axis(~ . / scale_factor, name = "Total Bases as red dots (100M)")
    ) +
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
    theme_classic()

## ---- duplication rate (fastp) ----
duplication_rates <- read_fastp_duplication_batch(
  file.path(data_dir, proj, "qc", paste0(proj, "_qc_files"), "fastp"),
  n_cpu)
p2 <- duplication_rates %>% 
  ggplot(aes(x = sample)) +
  geom_col(aes(y = duplication_rate)) +
  coord_flip() +
  xlab("Sample") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    name = "Duplication rate",
  ) +
  theme_classic()

## ---- Insert size distribution (fastp) ----
insert_size <- read_fastp_insert_size_batch(
  file.path(data_dir, proj, "qc", paste0(proj, "_qc_files"), "fastp"),
  n_cpu
) %>% 
  filter(count > 0)

p3 <- insert_size %>% 
  group_by(sample) %>% 
  arrange(insert_size) %>% 
  summarize(is25 = get_q_from_hist(0.25, insert_size, count),
            is50 = get_q_from_hist(0.50, insert_size, count),
            is75 = get_q_from_hist(0.75, insert_size, count),
            .groups = "drop") %>% 
  ggplot(aes(x = sample)) +
    geom_point(aes(y = is50)) +
    geom_errorbar(aes(ymin = is25, ymax = is75)) +
    coord_flip() +
    xlab("Sample") + 
    ylab("Insert Size (bp)")


## ---- Depth (qualimap) ----
depth_qualimap <- read_qualimap_depth_across_reference_batch(
  file.path(data_dir, proj, "qc", paste0(proj, "_qc_files"), "qualimap"),
  n_cpu
) %>% 
  arrange(sample)

if (bimodal) {
  depth_qualimap_summarized <- depth_qualimap %>% 
    group_by(sample) %>% 
    reframe(get_bimodal_means(depth))
} else{
  depth_qualimap_summarized <- depth_qualimap %>% 
    group_by(sample) %>% 
    reframe(get_unimodal_means(depth))
}

p5 <- depth_qualimap %>% 
  filter(depth <= quantile(depth, 0.99)) %>% 
  ggplot(aes(x = sample, y = depth)) +
    geom_violin() +
    geom_point(aes(x = sample, y = mean, color = mode),
               data = depth_qualimap_summarized) +
    coord_flip() +
    xlab("Sample") +
    ylab("Mapping Depth") +
    theme(legend.position = "none")

p7 <- depth_qualimap %>% 
  ggplot(aes(x = pos/10^6, y = depth, color = sample)) +
  geom_point(alpha = 0.5) +
  xlab("Position (Mbp)") +
  ylab("Mapping Depth") +
  theme_classic() +
  ylim(0, median(depth_qualimap_summarized$mean %>% max() * 1.5 ))

## ---- Coverage (bedtools) ----
coverage_bedtools <- read_bedtools_coverage_batch(
  file.path(data_dir, proj, "qc", paste0(proj, "_qc_files"), "coverage"),
  n_cpu
)

if (bimodal) {
  coverage_bedtools_summarized <- coverage_bedtools %>% 
    group_by(sample) %>% 
    reframe(get_bimodal_means(coverage))
} else{
  coverage_bedtools_summarized <- coverage_bedtools %>% 
    group_by(sample) %>% 
    reframe(get_unimodal_means(coverage))
}

p6 <- coverage_bedtools %>% 
  ggplot(aes(x = sample, y = coverage)) +
  geom_violin() +
  geom_point(aes(x = sample, y = mean, color = mode),
             data = coverage_bedtools_summarized) +
  coord_flip() +
  xlab("Sample") +
  ylab("Reference coverage") +
  theme(legend.position = "none")
  

p8 <- coverage_bedtools %>% 
  mutate(pos = start/2 + end/2) %>% 
  filter(str_detect(seqname, "chr")) %>% 
  ggplot(aes(x = pos/10^6, y = coverage, color = sample)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~seqname, scale = "free_x", nrow = 1) +
  xlab("Position (Mbp)") +
  ylab("Mapping Coverage") +
  theme_classic()

## ---- N_called sites ----
called_sites <- read_tsv(
  file.path(data_dir, proj, "qc", paste0(proj, "_qc_files"), paste0(proj, "_count_sites.tsv"))
)
  
called_sites_modif <- called_sites %>% 
  mutate(het_sites = ifelse(het_sites >= median(het_sites) * 1.5,
                            NA,
                            het_sites))

max_called_sites_M <- max(called_sites_modif$called_sites/1e6,  na.rm = TRUE)
max_het_sites_K    <- max(called_sites_modif$het_sites/1e3,  na.rm = TRUE)
scale_factor       <- max_called_sites_M / max_het_sites_K

called_sites_plot_df <- called_sites_modif %>%
  mutate(called_sites_M   = called_sites / 1e6,
         het_sites_K      = het_sites / 1e3,
         het_sites_scaled = het_sites_K * scale_factor)

p4 <- called_sites_plot_df %>% 
  ggplot(aes(x = sample)) +
  geom_col(aes(y = called_sites_M)) +
  geom_point(aes(y = het_sites_scaled), 
             size = 2, color = "grey") +
  geom_point(aes(y = max_het_sites_K * 1.05 * scale_factor), 
             size = 1, color = "red", 
             data = called_sites_plot_df %>% filter(is.na(het_sites))) +
  geom_point(aes(y = max_het_sites_K * 1.075 * scale_factor), 
             size = 1.5, color = "red", 
             data = called_sites_plot_df %>% filter(is.na(het_sites))) +
  geom_point(aes(y = max_het_sites_K * 1.1 * scale_factor), 
             size = 2, color = "red", 
             data = called_sites_plot_df %>% filter(is.na(het_sites))) +
  coord_flip() +
  xlab("Sample") +
  scale_y_continuous(
    name = "Total sites called (M)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Total heterozygous sites called as grey dots (K)")
  ) +
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme_classic()

## ---- Putting all together ----
p_1_to_6 <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3,legend = "none")

ggarrange(p_1_to_6, p7, p8, ncol = 1, nrow = 3, legend = "none", heights = c(3,1,1))

ggsave(file.path("mapping_qc", "out", paste0(proj, ".qc.tif")),
       height = 11, width = 8.5, units = "in")

filtered_counts %>% 
  left_join(insert_size %>% 
              group_by(sample) %>% 
              summarize(median_insert_size = median(insert_size)), 
            by = join_by(sample)) %>% 
  left_join(depth_qualimap_summarized %>% 
              dplyr::select(sample, mode, mean) %>% 
              mutate(mode = paste0(mode, "_mean_depth") %>% str_replace("-", "_")) %>% 
              pivot_wider(names_from = mode, values_from = mean),
            by = join_by(sample)) %>% 
  left_join(coverage_bedtools_summarized %>% 
              dplyr::select(sample, mode, mean) %>% 
              mutate(mode = paste0(mode, "_mean_coverage") %>% str_replace("-", "_")) %>% 
              pivot_wider(names_from = mode, values_from = mean),
            by = join_by(sample)) %>% 
  left_join(called_sites,
            by = join_by(sample)) %>% 
  write_tsv(file.path("mapping_qc", "out", paste0(proj, ".qc.tsv")))
  
