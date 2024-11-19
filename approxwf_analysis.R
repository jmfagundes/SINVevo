source("source_me.R")

# divide by cistron

cistrons <- list(UTR = c(1:59, 7602:7646, 11385:11703),
                 nsp1 = 60:1679, nsp2 = 1680:4100, nsp3 = 4101:5750, nsp4 = 5751:7601,
                 C = 7647:8438, E3 = 8439:8630, E2 = 8631:9899, `6K` = 9900:10064, E1 = 10065:11384)

approxwf.res$position <- approxwf.res$allele %>% gsub("[A-Z\\+\\*\\-]*", "", .) %>% as.numeric()
approxwf.res$cistron <- approxwf.res$position %>% lapply(function(x) names(cistrons)[lapply(cistrons, function(y) x %in% y) %>% unlist()]) %>% unlist()

set.seed(255)

s.ne.cistron.gg <- approxwf.res %>% 
  dplyr::select(treatment, mean_s, mean_Ne, cistron) %>%
  pivot_longer(c(mean_Ne, mean_s)) %>%
  ggplot(aes(cistron %>% factor(levels = names(cistrons)), value, color = treatment)) +
  xlab("Cistron") +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(position = position_jitterdodge(dodge.width = .75), alpha = .3, size = .5) +
  facet_wrap(~name %>% factor(levels = c("mean_s", "mean_Ne")), scales = "free_y",
             strip.position = "left",
             labeller = as_labeller(c(mean_Ne = "N\u2091", mean_s = "s")),
             ncol = 1) +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        #legend.position = "top",
        strip.text = element_text(face = "italic")) +
  scale_color_manual(values = scales::hue_pal()(2), 
                     labels = c("Gradual", "Sudden"),
                     name = "Treatment")

ggsave("1_s_ne.pdf", s.ne.cistron.gg, width = 6.85, height = 4.5)
ggsave("1_s_ne.png", s.ne.cistron.gg, width = 6.85, height = 4.5, dpi = 600)

# MANOVA

mean_s.Ne.manova <- manova(cbind(mean_s, mean_Ne) ~ treatment * cistron * as.factor(population), approxwf.res)

mean_s.Ne.manova.summary <- mean_s.Ne.manova %>% summary()

# s and Ne separately

mean_s.Ne.manova.summary.aov <- mean_s.Ne.manova %>% summary.aov()

# linear models for s and Ne

mean_s.lm <- lm(mean_s ~ treatment + cistron + as.factor(population), approxwf.res)
mean_Ne.lm <- lm(mean_Ne ~ treatment + cistron + as.factor(population), approxwf.res)

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

# proportion of CHO cells

cho.prop.gg <- data.frame(Passage = 0:25,
                          Sudden = c(0, rep(100, 25)),
                          Gradual = seq(0, 100, 4)) %>%
  pivot_longer(-Passage, names_to = "Treatment", values_to = "CHO cells (%)") %>%
  ggplot(aes(Passage, `CHO cells (%)`, color = Treatment)) +
  geom_line() +
  geom_vline(aes(xintercept = x), data = data.frame(x = c(4,7,10,13,16,19,22,25)), linetype = "dotted")

ggsave("cho.png", cho.prop.gg, width = 6.85, height = 3, dpi = 600)


