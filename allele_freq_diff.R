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

abs.diff.tb$Passage <- factor(abs.diff.tb$Passage, levels = c("0-4", "4-7", "7-10", "10-13", "13-16", "16-19", "19-22", "22-25"))

abs.diff.aov <- aov(lm(`Total AFD` ~ Treatment + Passage + rep , abs.diff.tb))

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

afd.gg <- ggarrange(afd.abs.gg + ggtitle("Concecutive passages") + theme(axis.title = element_blank()),
                    afd.anc.gg + ggtitle("Difference to ancestral") + theme(axis.title = element_blank()),
                    ncol = 1) %>% annotate_figure(left = "Total AFD", bottom = "Passage")

#afd.gg <- lapply(list("Gradual", "Sudden"), function(x) {
#  y <- list(`Concecutive passages` = abs.diff.tb,
#            `Difference to ancestral` = anc.diff.tb) %>% bind_rows(.id = "data") %>%
#  dplyr::filter(Treatment == x) %>%
#  ggplot() +
#  geom_boxplot(aes(Passage, `Total AFD`), outlier.size = 0) +
#  geom_line(aes(Passage, `Total AFD`, group = rep, color = rep)) +
#  facet_wrap(~data, scales = "free") +
#  theme(legend.position = "none",
#        axis.title = element_blank()) +
#  ggtitle(x)
#}) %>% ggarrange(plotlist = ., ncol = 1) %>% annotate_figure(left = "Total AFD", bottom = "Passage")

# consecutive vs difference to ancestral

afd.anc_consec.tb <- list(`Concecutive passages (total AFD)` = abs.diff.tb,
                          `Difference to ancestral (total AFD)` = anc.diff.tb) %>% bind_rows(.id = "data") %>%
  mutate(Passage := Passage %>% gsub(".*-", "", .) %>% as.numeric()) %>%
  pivot_wider(values_from = `Total AFD`, names_from = data)

#afd.and.diff.tb %>%
#  ggplot(aes(`Difference to ancestral (AFD)`, `Concecutive passages (AFD)`, color = rep, group = rep)) +
#  geom_point() +
#  facet_wrap(~Treatment) +
#  theme(legend.position = "none") +
#  geom_line() +
#  coord_equal()

afd.cor.gg <- afd.anc_consec.tb %>%
  ggplot(aes(`Difference to ancestral (total AFD)`, `Concecutive passages (total AFD)`, color = Treatment, group = rep)) +
  geom_point() +
  #theme(legend.position = "none") +
  geom_line() +
  coord_equal() +
  geom_smooth(aes(`Difference to ancestral (total AFD)`, `Concecutive passages (total AFD)`, color = Treatment, group = Treatment),
              data = afd.anc_consec.tb,
              se = FALSE, method = "lm", linetype = "dashed", linewidth = .5)

ggsave("2_afd.pdf", afd.gg, width = 6.85, height = 6)
ggsave("2_afd.png", afd.gg, width = 6.85, height = 6, dpi = 600)

# richness

richness.tb <- list(Gradual = gradual.mtx,
                    Sudden = sudden.mtx) %>% lapply(function(z) lapply(z, function(x) apply(x, 2, function(y) y[y != 0] %>% length()) %>%
                                                                         as.data.frame() %>% t() %>%
                                                                         as.data.frame() %>% pivot_longer(., cols = colnames(.))) %>%
                                                      bind_rows(.id = "Population")) %>% bind_rows(.id = "Treatment") %>%
  rename(Passage = name, Richness = value)

richness.boxplot.gg <- richness.tb %>%
  ggplot(aes(Passage %>% as.numeric(), Richness, group = interaction(Treatment, Passage), color = Treatment)) +
  geom_boxplot() +
  geom_line(aes(group = Population)) +
  facet_wrap(~Treatment) +
  ylab("Passage")

richness.mean.sd.total.afd.gg <- list(`Concecutive passages (total AFD)` = abs.diff.tb,
     `Difference to ancestral (total AFD)` = anc.diff.tb) %>% bind_rows(.id = "data") %>%
  mutate(Passage := Passage %>% gsub(".*-", "", .) %>% as.numeric()) %>%
  pivot_wider(values_from = `Total AFD`, names_from = data) %>% group_by(Treatment, Passage) %>%
  summarise(mean_consecutive = mean(`Concecutive passages (total AFD)`),
            mean_ancestral = mean(`Difference to ancestral (total AFD)`),
            sd_consecutive = sd(`Concecutive passages (total AFD)`),
            sd_ancestral = sd(`Difference to ancestral (total AFD)`)) %>%
  ggplot(aes(mean_ancestral, mean_consecutive, color = Treatment)) +
  geom_point() + geom_path()

# normalize AFD by richness

afd.anc_consec.richness.tb <- bind_cols(afd.anc_consec.tb, richness.tb["Richness"]) %>%
  dplyr::mutate(`Concecutive passages (mean AFD)` = `Concecutive passages (total AFD)` / Richness,
                `Difference to ancestral (mean AFD)` = `Difference to ancestral (total AFD)` / Richness)

mean.afd.gg <- afd.anc_consec.richness.tb %>%
  dplyr::select(Treatment, rep, Passage, `Concecutive passages (mean AFD)`) %>%
  dplyr::mutate(Passage := factor(paste0(Passage - 3, "-", Passage) %>% gsub("1-4", "0-4", .),
                                  levels = c("0-4", "4-7", "7-10", "10-13", "13-16", "16-19", "19-22", "22-25"))) %>%
  rename(`Mean AFD` = `Concecutive passages (mean AFD)`) %>%
  ggplot() +
  geom_boxplot(aes(Passage, `Mean AFD`), outlier.size = 0) +
  geom_line(aes(Passage, `Mean AFD`, group = rep, color = rep)) +
  facet_wrap(~Treatment) + theme(legend.position = "none")

mean.afd.cor.gg <- afd.anc_consec.richness.tb %>%
  ggplot(aes(`Difference to ancestral (mean AFD)`, `Concecutive passages (mean AFD)`, color = Treatment, group = rep)) +
  geom_point() +
  geom_line() +
  coord_equal() +
  geom_smooth(aes(`Difference to ancestral (mean AFD)`, `Concecutive passages (mean AFD)`, color = Treatment, group = Treatment),
              data = afd.anc_consec.richness.tb,
              se = FALSE, method = "lm", linetype = "dashed", linewidth = .5)

# wilcox

richness.22.25.gradual.tb <- richness.tb %>% dplyr::filter(Treatment == "Gradual" & Passage %in% c(22, 25)) %>%
  dplyr::mutate(Richness = Richness %>% as.numeric()) %>%
  pivot_wider(names_from = Passage, values_from = Richness)

richness.tb.wilcox <- wilcox.test(richness.22.25.gradual.tb$`22`, richness.22.25.gradual.tb$`25`, alternative = "greater", paired = TRUE)

