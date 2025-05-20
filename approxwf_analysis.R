source("source_me.R")

# divide by cistron

cistrons <- list(UTR = c(1:59, 7602:7646, 11385:11703),
                 nsp1 = 60:1679, nsp2 = 1680:4100, nsp3 = 4101:5750, nsp4 = 5751:7601,
                 C = 7647:8438, E3 = 8439:8630, E2 = 8631:9899, `6K` = 9900:10064, E1 = 10065:11384)

approxwf.res$position <- approxwf.res$allele %>% gsub("[A-Z\\+\\*\\-]*", "", .) %>% as.numeric()
approxwf.res$cistron <- approxwf.res$position %>% lapply(function(x) names(cistrons)[lapply(cistrons, function(y) x %in% y) %>% unlist()]) %>% unlist()

# add p-values

wilcox.p.tb <- approxwf.res %>% 
  dplyr::select(treatment, mean_s, mean_Ne, cistron) %>%
  pivot_longer(c(mean_Ne, mean_s)) %>% group_by(cistron, name) %>%
  summarise(wilcox.p = wilcox.test(value ~ treatment)$p.value) %>%
  dplyr::mutate(`p.signif` = case_when(wilcox.p <= 0.0001 ~ "****",
                                       wilcox.p <= 0.001 ~ "***",
                                       wilcox.p <= 0.01 ~ "**",
                                       wilcox.p <= 0.05 ~ "*",
                                       wilcox.p > 0.05 ~ NA),
                y = case_when(name == "mean_Ne" ~ max(approxwf.res$mean_Ne),
                              name == "mean_s" ~ max(approxwf.res$mean_s)))

set.seed(255)

s.ne.cistron.gg <- approxwf.res %>% 
  dplyr::select(treatment, mean_s, mean_Ne, cistron) %>%
  pivot_longer(c(mean_Ne, mean_s)) %>%
  ggplot(data = .) +
  xlab("Cistron") +
  geom_boxplot(aes(cistron %>% factor(levels = names(cistrons)), value, color = treatment), outliers = FALSE) +
  geom_jitter(aes(cistron %>% factor(levels = names(cistrons)), value, color = treatment),
              position = position_jitterdodge(dodge.width = .75), alpha = .3, size = .5) +
  facet_wrap(~name %>% factor(levels = c("mean_s", "mean_Ne")), scales = "free_y",
             strip.position = "left",
             labeller = as_labeller(c(mean_Ne = "N\u2091", mean_s = "s")),
             ncol = 1) +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        #legend.position = "top",
        strip.text = element_text(face = "italic", size = 12)) +
  scale_color_manual(values = scales::hue_pal()(2), 
                     labels = c("Gradual", "Sudden"),
                     name = "Treatment") +
  geom_text(aes(cistron, y = y, label = p.signif), wilcox.p.tb) +
  common.theme

ggsave("1_s_ne.pdf", s.ne.cistron.gg, width = 6.85, height = 4.5)
ggsave("1_s_ne.tiff", s.ne.cistron.gg, width = 6.85, height = 4.5, dpi = 600)

# MANOVA

mean_s.Ne.manova <- manova(cbind(mean_s, mean_Ne) ~ treatment * cistron * as.factor(population), approxwf.res)

mean_s.Ne.manova.summary <- mean_s.Ne.manova %>% summary()

# s and Ne separately

mean_s.Ne.manova.summary.aov <- mean_s.Ne.manova %>% summary.aov()

# linear models for s and Ne

mean_s.lm <- lm(mean_s ~ treatment + cistron + as.factor(population), approxwf.res)
mean_Ne.lm <- lm(mean_Ne ~ treatment + cistron + as.factor(population), approxwf.res)

# random effects and interaction

mean_s.lm.random <- nlme::lme(mean_s ~ treatment * cistron, random = ~ 1|as.factor(population), approxwf.res)
mean_Ne.lm.random <- nlme::lme(mean_Ne ~ treatment * cistron, random = ~1|as.factor(population), approxwf.res)

# cistron effect by treatment

gradual.mean_s.anova <- lm(mean_s ~ cistron * population, approxwf.res %>% dplyr::filter(treatment == "gradual")) %>% anova()
sudden.mean_s.anova <- lm(mean_s ~ cistron * population, approxwf.res %>% dplyr::filter(treatment == "sudden")) %>% anova()

# effect sizes

manova.effect_size <- effectsize::effectsize(mean_s.Ne.manova)
mean_s.anova.effect_size <- lm(mean_s ~ treatment + cistron + as.factor(population), approxwf.res) %>% anova() %>%
  effectsize::effectsize()
mean_Ne.anova.effect_size <- lm(mean_Ne ~ treatment + cistron + as.factor(population), approxwf.res) %>% anova() %>%
  effectsize::effectsize()

gradual.mean_s.anova.effect_size <- effectsize::effectsize(gradual.mean_s.anova)
sudden.mean_s.anova.effect_size <- effectsize::effectsize(sudden.mean_s.anova)

# schematics of methodology

# proportion of CHO cells

cho.prop.gg <- data.frame(Passage = 0:25,
                          Sudden = c(0, rep(100, 25)),
                          Gradual = seq(0, 100, 4)) %>%
  pivot_longer(-Passage, names_to = "Treatment", values_to = "CHO cells (%)") %>%
  ggplot(aes(Passage, `CHO cells (%)`, color = Treatment)) +
  geom_line() +
  geom_vline(aes(xintercept = x), data = data.frame(x = c(4,7,10,13,16,19,22,25)), linetype = "dotted") +
  common.theme

# plot frequency table

method.theme <- ttheme_default(base_size = 7)

example.alleles <- c("2984C", "3494C", "3915A", "5383G")

allele.freq.table.gg <- grid.arrange(tableGrob(sudden.mtx$`2`[example.alleles,] %>% round(3),
                                               theme = method.theme))

set.seed(2)

approxwf.gg <- approxwf.res %>% filter(population == 2 & allele %in% example.alleles) %>%
  dplyr::select(treatment, mean_s, mean_Ne, cistron) %>%
  pivot_longer(c(mean_Ne, mean_s)) %>%
  ggplot(data = .) +
  geom_boxplot(aes(treatment, value), outliers = FALSE) +
  geom_jitter(aes(treatment, value),
              size = .5) +
  facet_wrap(~name %>% factor(levels = c("mean_s", "mean_Ne")), scales = "free_y",
             strip.position = "left",
             labeller = as_labeller(c(mean_Ne = "N\u2091", mean_s = "s")),
             ncol = 1) +
  ggtitle("approxwf") +
  scale_y_continuous(n.breaks = 4) +
  common.theme +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "italic"))

TL.gg <- fit.TL(list(s2 = sudden.mtx$`2`),
                zero.rate.threshold = NULL, normalize = FALSE,
                remove.zeros = FALSE)$log_log.mean_sd$s2 %>%
  mutate(., allele = rownames(.), is.example = ifelse(allele %in% example.alleles, "yes", "no")) %>%
  ggplot(aes(mean, sd)) +
  geom_point(shape = 20, stroke = 0, size = 1.5, color = "grey") +
  #scale_color_manual(values = c("grey", "black")) +
  geom_text_repel(data = fit.TL(list(s2 = sudden.mtx$`2`[example.alleles,]),
                                zero.rate.threshold = NULL, normalize = FALSE,
                                remove.zeros = FALSE)$log_log.mean_sd$s2 %>%
                    mutate(., allele = rownames(.),
                           is.example = ifelse(allele %in% example.alleles, "yes", "no")),
                  aes(mean, sd, label = allele), seed = 2) +
  geom_point(data = fit.TL(list(s2 = sudden.mtx$`2`[example.alleles,]),
                           zero.rate.threshold = NULL, normalize = FALSE,
                           remove.zeros = FALSE)$log_log.mean_sd$s2 %>%
               mutate(., allele = rownames(.),
                      is.example = ifelse(allele %in% example.alleles, "yes", "no")),
             aes(mean, sd),
             shape = 20, stroke = 0, size = 1.5, color = "black") +
  xlab("Allele frequency mean (log)") + ylab("Allele frequency variance (log)") +
  geom_smooth(aes(mean, sd), method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, inherit.aes = FALSE) +
  theme(legend.position = "none") +
  ggtitle("Taylor's law") +
  common.theme

afd.gg <- data.frame(`AFD consecutive` = lapply(1:8, function(y) {
  x <- cbind("0" = 0, sudden.mtx$`2`[example.alleles,])[y:(y + 1)]
  x %>% calc.afd()
}) %>% unlist(),
`AFD to ancestral` = lapply(1:8, function(y) {
  x <- cbind("0" = 0, sudden.mtx$`2`[example.alleles,])[c(1, (y + 1))]
  x %>% calc.afd()
}) %>% unlist(),
`Richness` = c(1, 2, 2, 1, 3, 2, 3, 3),
Passage = c(4, 7, 10, 13, 16, 19, 22, 25), check.names = FALSE) %>%
  mutate(`Mean AFD (consecutive passages)` = `AFD consecutive` / Richness,
         `Mean AFD (difference to ancestral)` = `AFD to ancestral` / Richness) %>%
  select(Passage, `Mean AFD (consecutive passages)`, `Mean AFD (difference to ancestral)`, Richness) %>%
  pivot_longer(-Passage) %>%
  ggplot(aes(Passage, value)) + geom_line() + facet_wrap(~name, scales = "free_y", ncol = 1) +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("Diversity measures") +
  common.theme

# hurst

hurst.0.gg <- sudden.mtx$`2` %>% mutate(., Position = rownames(.) %>% gsub("\\D+", "", .) %>% as.numeric()) %>%
  pivot_longer(-Position) %>% filter(value != 0) %>% mutate(value := 1) %>% unique() %>%
  ggplot(aes(Position, value)) +
  geom_col(color = "black") +
  facet_wrap(~name %>% factor(levels = c(4, 7, 10, 13, 16, 19, 22, 25)), ncol = 1, strip.position = "left") +
  theme(axis.title.y = element_blank(),
        axis.text = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle("Convert to binary sequences")

hurst.2.gg <- fit.hurst(list(s2 = sudden.mtx$`2` %>% t()), d = 2000) %>% bind_rows(.id = "Passage") %>%
  mutate(Passage := factor(Passage, levels =  c(4, 7, 10, 13, 16, 19, 22, 25))) %>%
  ggplot(aes(Passage, He)) +
  geom_point() + ylab(expression(italic("H"))) + ggtitle(expression(paste("Estimation of ", italic("H")))) +
  common.theme

# save methodology fig

method.gg <- ggarrange(cho.prop.gg %>% annotate_figure(fig.lab = "a"),
                       ggarrange(allele.freq.table.gg, NULL, approxwf.gg,
                                 nrow = 1, widths = c(1, .2, .8)) %>% annotate_figure(fig.lab = "b"),
                       ggarrange(NULL, afd.gg, TL.gg,
                                 widths = c(.1, 1, .8),
                                 nrow = 1),
                       ggarrange(hurst.0.gg, NULL, hurst.2.gg,
                                 widths = c(1, .1, .6),
                                 nrow = 1),
                       nrow = 4, heights = c(.45, .55, .8, .6))

ggsave("0_methodology.pdf", method.gg, width = 6.85, height = 9.21)
ggsave("0_methodology.tiff", method.gg, width = 6.85, height = 9.21, dpi = 600)

