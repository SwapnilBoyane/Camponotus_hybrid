
#set wordking directory

setwd("/Volumes/T7/Chapter_1_2025/22_allele_freq")

#load library
library(tidyverse)

# Read allele frequency flies for each lineage and skip 1st line.
herc <- read_tsv("pure_Her.frq",
                 col_names = c("CHROM","POS","N_ALLELES","N_CHR_Herc","FREQ1_Herc","FREQ2_Herc"), skip = 1)
me   <- read_tsv("pure_ME.frq",
                 col_names = c("CHROM","POS","N_ALLELES","N_CHR_ME","FREQ1_ME","FREQ2_ME"), skip = 1)
mw   <- read_tsv("pure_MW.frq",
                 col_names = c("CHROM","POS","N_ALLELES","N_CHR_MW","FREQ1_MW","FREQ2_MW"), skip = 1)
nov  <- read_tsv("pure_Nov.frq",
                 col_names = c("CHROM","POS","N_ALLELES","N_CHR_Nov","FREQ1_Nov","FREQ2_Nov"), skip = 1)

# Merge all by CHROM and POS
af <- herc %>%
  select(CHROM, POS, N_CHR_Herc, FREQ2_Herc) %>%
  inner_join(select(me, CHROM, POS, N_CHR_ME, FREQ2_ME), by = c("CHROM","POS")) %>%
  inner_join(select(mw, CHROM, POS, N_CHR_MW, FREQ2_MW), by = c("CHROM","POS")) %>%
  inner_join(select(nov, CHROM, POS, N_CHR_Nov, FREQ2_Nov), by = c("CHROM","POS"))

# require ≥ 6 diploid individuals)
min_n <- 6      
fix_low  <- 0
fix_high <- 1

# Filter to sites with enough data in all lineages
af_filt <- af %>%
  filter(N_CHR_Herc >= min_n,
         N_CHR_ME   >= min_n,
         N_CHR_MW   >= min_n,
         N_CHR_Nov  >= min_n)


#  diagnostic SNPs per lineage

diag_herc <- af_filt %>% filter(
  (FREQ2_Herc >= fix_high & FREQ2_ME <= fix_low & FREQ2_MW <= fix_low & FREQ2_Nov <= fix_low) |
    (FREQ2_Herc <= fix_low & FREQ2_ME >= fix_high & FREQ2_MW >= fix_high & FREQ2_Nov >= fix_high)
)

diag_me <- af_filt %>% filter(
  (FREQ2_ME >= fix_high & FREQ2_Herc <= fix_low & FREQ2_MW <= fix_low & FREQ2_Nov <= fix_low) |
    (FREQ2_ME <= fix_low & FREQ2_Herc >= fix_high & FREQ2_MW >= fix_high & FREQ2_Nov >= fix_high)
)

diag_mw <- af_filt %>% filter(
  (FREQ2_MW >= fix_high & FREQ2_Herc <= fix_low & FREQ2_ME <= fix_low & FREQ2_Nov <= fix_low) |
    (FREQ2_MW <= fix_low & FREQ2_Herc >= fix_high & FREQ2_ME >= fix_high & FREQ2_Nov >= fix_high)
)

diag_nov <- af_filt %>% filter(
  (FREQ2_Nov >= fix_high & FREQ2_Herc <= fix_low & FREQ2_ME <= fix_low & FREQ2_MW <= fix_low) |
    (FREQ2_Nov <= fix_low & FREQ2_Herc >= fix_high & FREQ2_ME >= fix_high & FREQ2_MW >= fix_high)
)

# Make diagnostic sets disjoint 

diag_all <- bind_rows(
  diag_herc %>% mutate(Lineage = "Herc"),
  diag_me   %>% mutate(Lineage = "ME"),
  diag_mw   %>% mutate(Lineage = "MW"),
  diag_nov  %>% mutate(Lineage = "Nov")
) %>%
  mutate(ID = paste(CHROM, POS, sep=":")) %>%
  group_by(ID) %>%
  filter(n() == 1) %>%  # keep SNPs unique to one lineage
  ungroup()

table(diag_all$Lineage)

# Write BED per lineage
write.table(data.frame(CHROM=diag_herc$CHROM, START=diag_herc$POS-1, END=diag_herc$POS),
            "cutoff_1_diag_Herc.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(data.frame(CHROM=diag_me$CHROM, START=diag_me$POS-1, END=diag_me$POS),
            "cutoff_1_diag_ME.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(data.frame(CHROM=diag_mw$CHROM, START=diag_mw$POS-1, END=diag_mw$POS),
            "cutoff_1_diag_MW.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(data.frame(CHROM=diag_nov$CHROM, START=diag_nov$POS-1, END=diag_nov$POS),
            "cutoff_1_diag_Nov.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# write one table for all lineage
write.table(data.frame(CHROM=diag_all$CHROM, START=diag_all$POS-1, END=diag_all$POS, Lineage=diag_all$Lineage),
            "cutoff_1_diag_all_lineages.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
