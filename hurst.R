source("source_me.R")

# get positions all mutations

mut.pos <- list(gradual = gradual.mtx,
                sudden = sudden.mtx) %>% lapply(lapply, function(x) {
                  apply(x, 2, function(y) {
                    y <- y[!grepl("\\+|\\-", names(y))]
                    pos <- names(y)[ifelse(y > 0, TRUE, FALSE)] %>% gsub("\\D+", "", .) %>% as.numeric()
                    genome <- rep(0, 11703)
                    genome[pos] <- 1
                    genome
                  }) %>% t() %>% as.data.frame()
                })

h.exp <- list(gradual = fit.hurst(mut.pos$gradual, d = 2000) %>% bind_rows(.id = "time"),
              sudden = fit.hurst(mut.pos$sudden, d = 2000) %>% bind_rows(.id = "time")) %>% bind_rows(.id = "Treatment") %>%
  dplyr::mutate(Treatment := Treatment %>%
                  gsub("^g", "G", .) %>%
                  gsub("^s", "S", .))

# compute p-values

h.exp.wilcox.p <- h.exp %>% group_by(time) %>%
  summarise(wilcox.p = t.test(He ~ Treatment)$p.value) %>%
  dplyr::mutate(`p.signif` = case_when(wilcox.p <= 0.0001 ~ "****",
                                       wilcox.p <= 0.001 ~ "***",
                                       wilcox.p <= 0.01 ~ "**",
                                       wilcox.p <= 0.05 ~ "*",
                                       wilcox.p > 0.05 ~ NA),
                y = max(h.exp$He),
                `p.adj` = p.adjust(wilcox.p, method = "bonferroni"))

h.exp.gg <- ggplot(h.exp, aes(time %>% factor(levels = time %>% unique()), He, color = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_dodge(width = .75)) +
  xlab("Passage") + ylab( expression(italic("H"))) +
  geom_hline(aes(yintercept = .7), linetype = "dashed") +
  geom_text(aes(time, y = y, label = p.signif), h.exp.wilcox.p, inherit.aes = FALSE)

h.exp.facet_gg <- ggplot(h.exp, aes(time %>% factor(levels = time %>% unique()), He)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = data, color = data)) +
  xlab("Passage") + ylab( expression(italic("H"))) +
  geom_hline(aes(yintercept = .7), linetype = "dashed") +
  facet_wrap(~Treatment) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("4", "7"),
                                        c("7", "10"),
                                        c("10", "13"),
                                        c("13", "16"),
                                        c("16", "19"),
                                        c("19", "22"),
                                        c("22", "25")),
                     method = "wilcox.test", paired = TRUE, label = "p.signif", hide.ns = TRUE,
                     step.increase = 0, tip.length = 0, bracket.size = 0, vjust = 1)

# test effect of time + sample + data

h.exp.lm <- lm(He ~ Treatment + time + data, h.exp %>% mutate(time := time %>% as.numeric()))
h.exp.lm.random <- nlme::lme(He ~ Treatment + time, random = ~ 1|data, h.exp %>% mutate(time := time %>% as.numeric()))

# make a figure

ggsave("7_hurst1.pdf", h.exp.gg, width = 6.85, height = 4)
ggsave("7_hurst1.png", h.exp.gg, width = 6.85, height = 4, dpi = 600)

# hurst frequencies

h.freq.exp <- list(gradual = fit.hurst(gradual.mtx %>% lapply(function(x) as.data.frame(t(x))), d = 10, remove.zeros = TRUE) %>% bind_rows(.id = "time"),
                   sudden = fit.hurst(sudden.mtx %>% lapply(function(x) as.data.frame(t(x))), d = 10, remove.zeros = TRUE) %>% bind_rows(.id = "time")) %>% bind_rows(.id = "Treatment") %>%
  dplyr::mutate(Treatment := Treatment %>%
                  gsub("^g", "G", .) %>%
                  gsub("^s", "S", .))

h.freq.exp.wilcox.p <- h.freq.exp %>% group_by(time) %>%
  summarise(wilcox.p = wilcox.test(He ~ Treatment)$p.value) %>%
  dplyr::mutate(`p.signif` = case_when(wilcox.p <= 0.0001 ~ "****",
                                       wilcox.p <= 0.001 ~ "***",
                                       wilcox.p <= 0.01 ~ "**",
                                       wilcox.p <= 0.05 ~ "*",
                                       wilcox.p > 0.05 ~ NA),
                y = max(h.freq.exp$He))


# test effect of time + sample + data

h.freq.exp.lm <- lm(He ~ Treatment + time + data, h.freq.exp %>% mutate(time := time %>% as.numeric()))
h.freq.exp.lm.random <- nlme::lme(He ~ Treatment + time, random = ~ 1|data, h.freq.exp %>% mutate(time := time %>% as.numeric()))

h.freq.exp.gg <- ggplot(h.freq.exp, aes(time %>% factor(levels = time %>% unique()), He, color = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_dodge(width = .75)) +
  xlab("Passage") + ylab( expression(italic("H"))) +
  geom_hline(aes(yintercept = .7), linetype = "dashed") +
  geom_text(aes(time, y = y, label = p.signif), h.freq.exp.wilcox.p, inherit.aes = FALSE)
