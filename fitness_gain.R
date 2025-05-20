# fitness_gain.R
################################################################################
library(readxl)
library(ggsci)
library(scales)
library(dplyr)
library(tidyverse)
#library(gtx)
library(knitr)
library(kableExtra)
library(magrittr)
library(reporter)
library(ggstatsplot)
library(ggpubr)

################################################################################
# LOADING DATA AND FORMAT ORIGINAL DATA (MORLEY, 2015)
################################################################################
treatment.colors <- c("Control" = "#2E2A2BFF",
                      "Sudden" = "#00BFC4",
                      "Gradual" = "#F8766D",
                      "Ancestral" = "gray")

# Original data from Morley, 2015
fitness <- read.csv("Fitness_Assay_Data.csv",
                    sep = ",")

fitness$Treatment[is.na(fitness$Treatment)] <- 0

# Only interested treatments:
fitness <- fitness %>%
  filter(Treatment %in% c(0, 1, 5, 6) & !is.na(Passage)) %>%
  # Calculate mean values from block measurement
  rowwise() %>%
  mutate(BHK_titer.mean = mean(BHK_Titer_Block1, BHK_Titer_Block2, BHK_Titer_Block3),
         CHO_titer.mean = mean(CHO_Titer_Block1, CHO_Titer_Block2, CHO_Titer_Block3),
         Passage = factor(Passage),
         Treatment = factor(Treatment))

# Fitness of the ANCESTRAL VIRUS in each cell type
anc.BHK.fitness.mean <- mean(fitness %>%
                               filter(Treatment == 0) %>%
                               pull(BHK_titer.mean))
anc.CHO.fitness.mean <- mean(fitness %>%
                               filter(Treatment == 0) %>%
                               pull(CHO_titer.mean))

fitness <- fitness %>%
  rowwise() %>%
  # Add changes in fitness substracting values of ancestral fitness
  mutate(change.fitness.BHK = BHK_titer.mean - anc.BHK.fitness.mean,
         change.fitness.CHO = CHO_titer.mean - anc.CHO.fitness.mean)


################################################################################
# PLOT fitness gain per treatment
################################################################################
# By Treatment, paired by population
fitness %>%
  dplyr::select(c("Population", "Treatment", "BHK_titer.mean", "CHO_titer.mean", 
                  "change.fitness.BHK", "change.fitness.CHO")) %>%
  pivot_longer(cols = c("change.fitness.BHK", "change.fitness.CHO"),
               names_to = "cell.type",
               values_to = "changes") %>%
  filter(Treatment != 0) %>%
  mutate(Treatment = factor(Treatment, 
                            levels = c(0, 5, 1, 6), 
                            labels = c("Ancestral", "Gradual", "Sudden", "Control"))) %>%
  
  ggplot(aes(x = cell.type, 
             y = changes)) +
  geom_hline(yintercept = 0,
             linewidth = 0.2,
             lty = "dotted") +
  # Distribution by treatment
  geom_boxplot(aes(fill = Treatment), 
               color = "black",
               alpha = 0.3, 
               width = 0.3,
               linewidth = 0.2,
               outlier.size = 0) +
  # Lines connecting dots for each population (lineage)
  geom_line(aes(color = Treatment,
                group = Population),
            linewidth = 0.1) +
  # Points
  geom_jitter(aes(
    color = Treatment, 
    group = interaction(cell.type, Treatment),
    shape = Treatment), 
    position = position_jitterdodge(jitter.width = 0.1, 
                                    jitter.height = 0,
                                    dodge.width = 0.3), 
    alpha = 0.7, 
    size = 1, 
    show.legend = FALSE) +
  facet_wrap(~Treatment) +
  # color by Treatment
  scale_color_manual(values = treatment.colors) +
  scale_fill_manual(values = treatment.colors) +
  scale_x_discrete(labels = c("change.fitness.BHK" = "BHK", 
                              "change.fitness.CHO" = "CHO")) +  
  labs(#title = "Effect of the evolution in the fitness",
    #subtitle = "(fitness"[evolved] - "fitness"[ancestral] ~ ")",
    y = expression("Changes in fitness (pfu mL"^-1 ~ ")"),
    x = "Cell type",
    #caption = "Wilcoxon rank test (pasired by lineages)"
  ) + 
  theme(#panel.border = element_rect(color = "black", fill = "transparent"),
    axis.text = element_text(color = "black"),
  ) +
  stat_compare_means(aes(x = cell.type, 
                         y = changes), comparisons = list(c("change.fitness.BHK", "change.fitness.CHO")),
                     method = "wilcox.test", paired = TRUE, label = "p.signif", hide.ns = TRUE,
                   step.increase = 0, tip.length = 0, bracket.size = 0, vjust = 1) +
  common.theme + theme(strip.text = element_text(size = 12))

ggsave("3_fitness.pdf",
       width = 6.85,
       height = 2.5)

ggsave("3_fitness.tiff",
       dpi = 600,
       width = 6.85,
       height = 2.5)


################################################################################
# STATS fitness gain per treatment
################################################################################
## Wilcoxon signed rank 
wilcox.p.tb <- fitness %>%
  filter(Treatment != 0) %>%  # Exclude the ancestral group (Treatment == 0)
  mutate(Treatment = factor(Treatment, 
                            levels = c(6, 1, 5), 
                            labels = c("Control", "Sudden", "Gradual"))) %>%
  group_by(Treatment) %>%
  summarise(
    # Perform the paired Wilcoxon test and extract relevant components
    wilcox_result = list(wilcox.test(change.fitness.BHK, 
                                     change.fitness.CHO, 
                                     paired = TRUE)),
    
    # Extract the p-value from the Wilcoxon test result
    pvalue = wilcox_result[[1]]$p.value,
    
    median_inBHK = median(change.fitness.BHK),
    median_inCHO = median(change.fitness.CHO),
    
    # Median difference between BHK and CHO (Estimate)
    median_difference = median_inCHO - median_inBHK,
    
    # Calculate fold change (CHO/BHK)
    fold_change = mean(change.fitness.CHO) / mean(change.fitness.BHK),
    
    # Calculate log2 fold change
    log2_fold_change = log2(mean(change.fitness.CHO) / mean(change.fitness.BHK)),
    
    # Z-statistic from the Wilcoxon test result
    Z = wilcox_result[[1]]$statistic,  # Z-statistic
    
    # Sample size
    n = n(),
    
    # Effect size (r) from the Z statistic
    r = Z / sqrt(n),
    
    # Interquartile Range (IQR) of fitness change
    IQR_BHK = IQR(change.fitness.BHK, na.rm = TRUE),
    IQR_CHO = IQR(change.fitness.CHO, na.rm = TRUE),
    
    # Range of the fitness change
    range_BHK = max(change.fitness.BHK) - min(change.fitness.BHK),
    range_CHO = max(change.fitness.CHO) - min(change.fitness.CHO)
  ) %>% 
  dplyr::select(c("Treatment", "pvalue", "Z", "median_difference", "fold_change", "log2_fold_change", "r")) #%>%
  kable(col.names = c("Treatment", "pvalue", "Z", "median_difference", "fold_change", "log2_fold_change", "r"),
        caption = "Wilcoxon signed rank (repeated measures)", 
        escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "striped", "striped", "hover", "condensed"),
                full_width = FALSE)
