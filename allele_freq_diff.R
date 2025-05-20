source("source_me.R")

abs.diff.tb <- list(Gradual = gradual.mtx %>%
                      lapply(function(x) cbind("0" = 0, x)),
                    Sudden = sudden.mtx %>%
                      lapply(function(x) cbind("0" = 0, x))) %>% lapply(function(treatment) {
                        
                        df <- lapply(treatment, function(x) {
                          lapply(1:(ncol(x) - 1), function(y) {
                            x <- x[y:(y + 1)]
                            x %>% calc.afd()
                          }) %>% unlist()
                        }) %>% bind_rows() %>% t() %>% as.data.frame()
                        
                        colnames(df) <- c("0-4", "4-7", "7-10", "10-13", "13-16", "16-19", "19-22", "22-25")
                        df$rep <- row.names(df)
                        df %>% pivot_longer(-rep, names_to = "Passage", values_to = "Total AFD")
                        
                      }) %>% bind_rows(.id = "Treatment")

# number of alleles that are varying from one point to another

consec.nrow.tb <- list(Gradual = gradual.mtx %>%
                         lapply(function(x) cbind("0" = 0, x)),
                       Sudden = sudden.mtx %>%
                         lapply(function(x) cbind("0" = 0, x))) %>% lapply(function(treatment) {
                           
                           df <- lapply(treatment, function(x) {
                             lapply(1:(ncol(x) - 1), function(y) {
                               x <- x[y:(y + 1)]
                               x[rowSums(x) > 0,] %>% nrow() 
                             }) %>% unlist()
                           }) %>% bind_rows() %>% t() %>% as.data.frame()
                           
                           colnames(df) <- c("0-4", "4-7", "7-10", "10-13", "13-16", "16-19", "19-22", "22-25")
                           df$rep <- row.names(df)
                           df %>% pivot_longer(-rep, names_to = "Passage", values_to = "nrow")
                           
                         }) %>% bind_rows(.id = "Treatment")

abs.diff.tb$Passage <- factor(abs.diff.tb$Passage, levels = c("0-4", "4-7", "7-10", "10-13", "13-16", "16-19", "19-22", "22-25"))

afd.abs.gg <- abs.diff.tb %>% ggplot() +
  geom_boxplot(aes(Passage, `Total AFD`), outlier.size = 0) +
  geom_line(aes(Passage, `Total AFD`, group = rep, color = rep)) +
  facet_wrap(~Treatment) + theme(legend.position = "none")

# difference from ancestral

anc.diff.tb <- list(Gradual = gradual.mtx %>%
                      lapply(function(x) cbind("0" = 0, x)),
                    Sudden = sudden.mtx %>%
                      lapply(function(x) cbind("0" = 0, x))) %>% lapply(function(treatment) {
                        
                        df <- lapply(treatment, function(x) {
                          lapply(1:(ncol(x) - 1), function(y) {
                            x <- x[c(1, (y + 1))]
                            x %>% calc.afd()
                          }) %>% unlist()
                        }) %>% bind_rows() %>% t() %>% as.data.frame()
                        
                        colnames(df) <- c("0-4", "0-7", "0-10", "0-13", "0-16", "0-19", "0-22", "0-25")
                        df$rep <- row.names(df)
                        df %>% pivot_longer(-rep, names_to = "Passage", values_to = "Total AFD")
                        
                      }) %>% bind_rows(.id = "Treatment")

anc.diff.tb$Passage <- factor(anc.diff.tb$Passage, levels = c("0-4", "0-7", "0-10", "0-13", "0-16", "0-19", "0-22", "0-25"))

afd.anc.gg <- anc.diff.tb %>% ggplot() +
  geom_boxplot(aes(Passage, `Total AFD`), outlier.size = 0) +
  geom_line(aes(Passage, `Total AFD`, group = rep, color = rep)) +
  facet_wrap(~Treatment) + theme(legend.position = "none")

afd.gg <- ggarrange(afd.abs.gg + ggtitle("Concecutive passages") +
                      theme(axis.title = element_blank(), title = element_text(size = 6)),
                    afd.anc.gg + ggtitle("Difference to ancestral") +
                      theme(axis.title = element_blank(), title = element_text(size = 6)),
                    ncol = 1) %>% annotate_figure(left = "Total AFD", bottom = "Passage")

# consecutive vs difference to ancestral

afd.anc_consec.tb <- list(`Concecutive passages (total AFD)` = abs.diff.tb,
                          `Difference to ancestral (total AFD)` = anc.diff.tb) %>% bind_rows(.id = "data") %>%
  mutate(Passage := Passage %>% gsub(".*-", "", .) %>% as.numeric()) %>%
  pivot_wider(values_from = `Total AFD`, names_from = data)

afd.cor.gg <- afd.anc_consec.tb %>%
  ggplot(aes(`Difference to ancestral (total AFD)`, `Concecutive passages (total AFD)`, color = Treatment, group = rep)) +
  geom_point() +
  geom_path() +
  coord_equal() +
  geom_smooth(aes(`Difference to ancestral (total AFD)`, `Concecutive passages (total AFD)`, color = Treatment, group = Treatment),
              data = afd.anc_consec.tb,
              se = FALSE, method = "lm", linetype = "dashed", linewidth = .5)

# richness

richness.tb <- list(Gradual = gradual.mtx,
                    Sudden = sudden.mtx) %>% lapply(function(z) lapply(z, function(x) apply(x, 2, function(y) y[y != 0] %>% length()) %>%
                                                                         as.data.frame() %>% t() %>%
                                                                         as.data.frame() %>% pivot_longer(., cols = colnames(.))) %>%
                                                      bind_rows(.id = "Population")) %>% bind_rows(.id = "Treatment") %>%
  rename(Passage = name, Richness = value)

richness.boxplot.gg <- richness.tb %>%
  ggplot(aes(Passage %>% factor(levels = c("4", "7", "10", "13", "16", "19", "22", "25")), Richness, group = Passage)) +
  geom_boxplot(outlier.size = 0) +
  geom_line(aes(group = Population, color = Population)) +
  facet_wrap(~Treatment) +
  theme(legend.position = "none") +
  ylab("Richness") + xlab("Passage") +
  stat_compare_means(comparisons = list(c("4", "7"),
                                        c("7", "10"),
                                        c("10", "13"),
                                        c("13", "16"),
                                        c("16", "19"),
                                        c("19", "22"),
                                        c("22", "25")),
                     method = "wilcox.test", paired = TRUE, label = "p.signif", hide.ns = TRUE,
                     step.increase = 0, tip.length = 0, bracket.size = 0, vjust = 1) +
  common.theme + theme(strip.text = element_text(size = 12))

# normalize AFD by richness

afd.anc_consec.richness.tb <- bind_cols(afd.anc_consec.tb, richness.tb["Richness"]) %>%
  dplyr::mutate(`Concecutive passages (mean AFD)` = `Concecutive passages (total AFD)` / consec.nrow.tb$nrow,
                `Difference to ancestral (mean AFD)` = `Difference to ancestral (total AFD)` / Richness)

# make figures

mean.afd.consec.gg <- afd.anc_consec.richness.tb %>%
  dplyr::select(Treatment, rep, Passage, `Concecutive passages (mean AFD)`) %>%
  dplyr::mutate(Passage := factor(paste0(Passage - 3, "-", Passage) %>% gsub("1-4", "0-4", .),
                                  levels = c("0-4", "4-7", "7-10", "10-13", "13-16", "16-19", "19-22", "22-25"))) %>%
  rename(`Mean AFD` = `Concecutive passages (mean AFD)`) %>%
  ggplot(aes(Passage, `Mean AFD`)) +
  geom_boxplot(outlier.size = 0) +
  geom_line(aes(group = rep, color = rep)) +
  facet_wrap(~Treatment) + theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("0-4", "4-7"),
                                        c("4-7", "7-10"),
                                        c("7-10", "10-13"),
                                        c("10-13", "13-16"),
                                        c("13-16", "16-19"),
                                        c("16-19", "19-22"),
                                        c("19-22", "22-25")),
                     method = "wilcox.test", paired = TRUE, label = "p.signif", hide.ns = TRUE,
                     step.increase = 0, tip.length = 0, bracket.size = 0, vjust = 1) +
  common.theme + theme(axis.text.x = element_text(angle = -60, vjust = 1, hjust = 0), strip.text = element_text(size = 12))

mean.afd.anc.gg <- afd.anc_consec.richness.tb %>%
  dplyr::select(Treatment, rep, Passage, `Difference to ancestral (mean AFD)`) %>%
  dplyr::mutate(Passage := factor(paste0("0-", Passage),
                                  levels = c("0-4", "0-7", "0-10", "0-13", "0-16", "0-19", "0-22", "0-25"))) %>%
  rename(`Mean AFD` = `Difference to ancestral (mean AFD)`) %>%
  ggplot(aes(Passage, `Mean AFD`)) +
  geom_boxplot(outlier.size = 0) +
  geom_line(aes(group = rep, color = rep)) +
  facet_wrap(~Treatment) + theme(legend.position = "none")+
  stat_compare_means(comparisons = list(c("0-4", "0-7"),
                                        c("0-7", "0-10"),
                                        c("0-10", "0-13"),
                                        c("0-13", "0-16"),
                                        c("0-16", "0-19"),
                                        c("0-19", "0-22"),
                                        c("0-22", "0-25")),
                     method = "wilcox.test", paired = TRUE, label = "p.signif", hide.ns = TRUE,
                     step.increase = 0, tip.length = 0, bracket.size = 0, vjust = 1) +
  common.theme + theme(axis.text.x = element_text(angle = -60, vjust = 1, hjust = 0), strip.text = element_text(size = 12))
mean.afd.cor.gg <- afd.anc_consec.richness.tb %>%
  ggplot(aes(`Difference to ancestral (mean AFD)`, `Concecutive passages (mean AFD)`, color = Treatment, group = rep)) +
  geom_point(data = afd.anc_consec.richness.tb %>% filter(Passage == 4)) +
  geom_path(arrow = arrow(angle = 20, type = "closed", length = unit(.2, "cm"))) +
  coord_equal() +
  geom_smooth(aes(`Difference to ancestral (mean AFD)`, `Concecutive passages (mean AFD)`, color = Treatment, group = Treatment),
              data = afd.anc_consec.richness.tb,
              se = FALSE, method = "lm", linetype = "dashed", linewidth = .5)

total.afd.cor.gg <- afd.anc_consec.richness.tb %>%
  ggplot(aes(`Difference to ancestral (total AFD)`, `Concecutive passages (total AFD)`, color = Treatment, group = rep)) +
  geom_point(data = afd.anc_consec.richness.tb %>% filter(Passage == 4)) +
  geom_path(arrow = arrow(angle = 20, type = "closed", length = unit(.2, "cm"))) +
  coord_equal() +
  geom_smooth(aes(`Difference to ancestral (total AFD)`, `Concecutive passages (total AFD)`, color = Treatment, group = Treatment),
              data = afd.anc_consec.richness.tb,
              se = FALSE, method = "lm", linetype = "dashed", linewidth = .5)

mean.afd.gg <- ggarrange(mean.afd.consec.gg + ggtitle("Concecutive passages") +
                           theme(axis.title = element_blank(), title = element_text(size = 10)),
                         mean.afd.anc.gg + ggtitle("Difference to ancestral") +
                           theme(axis.title = element_blank(), title = element_text(size = 10)),
                         ncol = 1) %>% annotate_figure(left = "Mean AFD", bottom = "Passage")

mean.afd.richness.gg <- ggarrange(mean.afd.gg %>% annotate_figure(fig.lab = "a"),
                                  (richness.boxplot.gg + theme(axis.title = element_blank())) %>%
                                    annotate_figure(left = "Richness", bottom = "Passage") %>%
                                    annotate_figure(fig.lab = "b"), ncol = 1, heights = (c(1, .5)))

ggsave("2_afd.pdf", mean.afd.richness.gg, width = 6.85, height = 8)
ggsave("2_afd.tiff", mean.afd.richness.gg, width = 6.85, height = 8, dpi = 600)

# without last passage gradual

total.afd.cor.minus_gradual25.gg <- afd.anc_consec.richness.tb %>%
  filter(!(Passage == 25 & Treatment == "Gradual")) %>%
  ggplot(aes(`Difference to ancestral (total AFD)`, `Concecutive passages (total AFD)`, color = Treatment, group = rep)) +
  geom_point(data = afd.anc_consec.richness.tb %>% filter(Passage == 4)) +
  geom_path(arrow = arrow(angle = 20, type = "closed", length = unit(.2, "cm"))) +
  coord_equal() +
  geom_smooth(aes(`Difference to ancestral (total AFD)`, `Concecutive passages (total AFD)`, color = Treatment, group = Treatment),
              se = FALSE, method = "lm", linetype = "dashed", linewidth = .5)

mean.afd.cor.minus_gradual25.gg <- afd.anc_consec.richness.tb %>%
  filter(!(Passage == 25 & Treatment == "Gradual")) %>%
  ggplot(aes(`Difference to ancestral (mean AFD)`, `Concecutive passages (mean AFD)`, color = Treatment, group = rep)) +
  geom_point(data = afd.anc_consec.richness.tb %>% filter(Passage == 4)) +
  geom_path(arrow = arrow(angle = 20, type = "closed", length = unit(.2, "cm"))) +
  coord_equal() +
  geom_smooth(aes(`Difference to ancestral (mean AFD)`, `Concecutive passages (mean AFD)`, color = Treatment, group = Treatment),
              se = FALSE, method = "lm", linetype = "dashed", linewidth = .5)

# richness vs mean AFD

afd.anc_consec.richness.gg <- afd.anc_consec.richness.tb %>%
  ggplot(aes(`Concecutive passages (mean AFD)`, Richness, color = Passage)) +
  geom_point(data = afd.anc_consec.richness.tb %>% filter(Passage == 4)) +
  facet_wrap(~Treatment) +
  geom_path(aes(group = rep), arrow = arrow(angle = 20, type = "closed", length = unit(.2, "cm")))

afd.anc_consec.richness.density <- afd.anc_consec.richness.tb %>%
  ggplot(aes(`Concecutive passages (mean AFD)`, Richness, color = Treatment)) +
  geom_density_2d(aes(`Concecutive passages (mean AFD)`, Richness, color = Treatment)) +
  geom_point()

