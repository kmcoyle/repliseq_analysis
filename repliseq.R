## Replication timing analysis

library(tidyverse)

repliseq_dir <- "/projects/rmorin_scratch/krysta_temp/repliseq/"
metadata <- read.table(paste0(repliseq_dir, "samples.tsv"), header=TRUE)

bed_names <- c("chr","start","end")
full_bed <- paste0(repliseq_dir, "B-merge_bam_coverage/50kb.hg38.bed")
bed_coverage <- read.table(full_bed, col.names = c(bed_names, "tags"))

sample_coverage <- apply(metadata, 1, function(row) {
  sample_id<- row[1]
  phase <- row[2]
  sample_table <- read.table(paste0(repliseq_dir, "A-150bp_coverage/",sample_id, ".50kb.hg38.bed"), col.names = c(bed_names, phase))
}) %>% purrr::reduce(left_join)

sample_coverage <- left_join(bed_coverage, sample_coverage) %>% mutate(
  G2.pctg = G2/tags,
  G1b.pctg = G1b/tags, 
  S1.pctg = S1/tags,
  S2.pctg = S2/tags,
  S3.pctg = S3/tags,
  S4.pctg = S4/tags,
) %>% select(-G2,-G1b, -tags, -S1, -S2, -S3, -S4)

coverage_long <- sample_coverage %>% gather(., key="phase", value="pctg", G2.pctg:S4.pctg)

ggplot(coverage_long %>% filter(chr == "chr1")) + 
  geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=pctg, fill=phase))+
  facet_wrap(facets=~phase, ncol = 1)+
  theme_bw()
