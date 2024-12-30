source("source_me.R")

# pipe function

pipe <- function(gradual.mtx,
                 sudden.mtx,
                 remove.zeros = FALSE) {
  
  gradual.TL <- fit.TL(gradual.mtx,
                       zero.rate.threshold = NULL, normalize = FALSE,
                       remove.zeros = remove.zeros)
  sudden.TL <- fit.TL(sudden.mtx,
                      zero.rate.threshold = NULL, normalize = FALSE,
                      remove.zeros = remove.zeros)
  
  # plot parameters
  
  snps.TL.params <- rbind(gradual.TL$params[-1] %>% cbind(data = "gradual"),
                          sudden.TL$params[-1] %>% cbind(data = "sudden"))

  snps.TL.param.gg <- ggplot(snps.TL.params) +
    geom_point(aes(V, beta, color = data)) +
    geom_smooth(aes(V, beta, color = data), method = "lm", formula = y ~ x)
  
  # Taylor's plot
  
  snps.TP.gg <- mapply(function(x, y, z) {
    plot.power_law(x,
                   y %>% calc.0_rate(),
                   metadata.range = NULL, legend.title = "0 rate", l10 = TRUE) +
      xlab("mean (log)") + ylab("standard deviation (log)") +
      geom_smooth(method = "lm", se = FALSE, linewidth = .1) 
      
  },
  list(gradual = gradual.TL$log_log.mean_sd,
       sudden = sudden.TL$log_log.mean_sd),
  list(gradual = gradual.mtx,
       sudden = sudden.mtx),
  list("gradual", "sudden"),
  SIMPLIFY = FALSE)
  
  list(TL = list(gradual = gradual.TL,
                 sudden = sudden.TL),
       TL.gg = snps.TL.param.gg,
       snps.gg = snps.TP.gg)
}

# fit SNPs frequencies to Taylor's law

TL.classic <- pipe(gradual.mtx, sudden.mtx, FALSE)
TL.classic$snps.gg <- lapply(TL.classic$snps.gg, function(x) {
  x + xlim(-3, 0) + ylim(-3.2, -.45) +
    theme(legend.position = "none",
          strip.text = element_text(size = 6)) +
    geom_smooth(aes(mean, sd), method = "lm", se = FALSE, color = "black", linewidth = .1) +
    scale_color_discrete(direction = -1) 
})

# Taylor's law per 0 rate

pipe.per_0 <- function(n = 1:8,
                       remove.zeros = FALSE,
                       min.rows = 5) {
  
  lapply(setNames(n, nm = as.character(n)), function(x) {
    
    # fit SNPs frequencies to Taylor's law
    
    gradual.TL <- fit.TL(gradual.mtx %>% lapply(apply.filter, x, exact = TRUE),
                         zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = remove.zeros, min.rows = min.rows)
    sudden.TL <- fit.TL(sudden.mtx %>% lapply(apply.filter, x, exact = TRUE),
                        zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = remove.zeros, min.rows = min.rows)
    
    # plot parameters
    
    snps.TL.params <- rbind(gradual.TL$params[-1] %>% cbind(data = "gradual"),
                            sudden.TL$params[-1] %>% cbind(data = "sudden"))
    
    snps.TL.param.gg <- ggplot(snps.TL.params, aes(V, beta, color = data)) +
      geom_point() +
      geom_smooth(method = "lm", formula = y ~ x)
    
    # Taylor's plot
    
    snps.TP.gg <- mapply(function(x, y, z) {
      plot.power_law(x,
                     y %>% calc.0_rate(),
                     metadata.range = c(0, 1), legend.title = "0 rate", l10 = TRUE) +
        xlab("mean (log)") + ylab("standard deviation (log)") +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red") +
        ggtitle(z)
    },
    list(gradual = gradual.TL$log_log.mean_sd,
         sudden = sudden.TL$log_log.mean_sd),
    list(gradual = gradual.mtx,
         sudden = sudden.mtx),
    list("gradual", "sudden"), SIMPLIFY = FALSE)
    
    list(TL = list(gradual = gradual.TL,
                   sudden = sudden.TL),
         TL.gg = snps.TL.param.gg,
         snps.gg = snps.TP.gg)
  })
}

TL.per.zero <- pipe.per_0()

# some random walk matrices

random.walk.tb <- lapply(setNames(nm = 1:1000), function(x) {
  set.seed(x)
  y <- cumsum(rnorm(8, 0, .2))
  y[y < .01] <- 0
  y[y > .99] <- 1
  
  if (sum(y) == 0) return()
  y
}) %>% bind_rows() %>% t() %>% as.data.frame()

random.walk.tb.2 <- lapply(setNames(nm = 1:1000), function(x) {
  set.seed(x)
  y <- cumsum(rnorm(8, 0, .05))
  y[y < .01] <- 0
  y[y > .99] <- 1
  
  if (sum(y) == 0) return()
  y
}) %>% bind_rows() %>% t() %>% as.data.frame()

random.walk.tb.3 <- lapply(setNames(nm = 1:1000), function(x) {
  set.seed(x)
  y <- cumsum(rnorm(8, 0, .01))
  y[y < .01] <- 0
  y[y > .99] <- 1
  
  if (sum(y) == 0) return()
  y
}) %>% bind_rows() %>% t() %>% as.data.frame()

random.TL.classic <- fit.TL(list(`Random 0.2` = random.walk.tb,
                                 `Random 0.05` = random.walk.tb.2,
                                 `Random 0.01` = random.walk.tb.3),
                            zero.rate.threshold = NULL, normalize = FALSE,
                            remove.zeros = FALSE)

pipe.per_0.random <- function(n = 1:8,
                              remove.zeros = FALSE,
                              min.rows = 5) {
  
  lapply(setNames(n, nm = as.character(n)), function(x) {
    
    # fit SNPs frequencies to Taylor's law
    
    fit.TL(list(random.2 = random.walk.tb,
                random.05 = random.walk.tb.2,
                random.01 = random.walk.tb.3) %>% lapply(apply.filter, x, exact = TRUE),
           zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = remove.zeros, min.rows = min.rows)
  })
}

random.TL.per_0 <- pipe.per_0.random()

# join the plots

V.beta.per_0.tb <- lapply(TL.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% c(lapply(random.TL.per_0, function(x) {
  x$params %>% mutate(treatment = "random")
})) %>% bind_rows(.id = "n")

facet.V.beta.gg <- V.beta.per_0.tb %>% mutate(treatment := factor(treatment %>%
                                                                    gsub("^g", "G", .) %>%
                                                                    gsub("^s", "S", .) %>%
                                                                    gsub("^r", "R", .),
                                                                  levels = c("gradual", "sudden", "random") %>%
                                                                    gsub("^g", "G", .) %>%
                                                                    gsub("^s", "S", .) %>%
                                                                    gsub("^r", "R", .))) %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) + facet_wrap(~treatment) +
  ylab(expression(italic("\u03B2"[n]))) + xlab(expression(italic("V"[n])))

# where are the SNPs?

location.SNPs.gg <- lapply(list(gradual = gradual.mtx,
                                sudden = sudden.mtx), function(x) {
                                  lapply(x, function(y) {
                                   y <- as.matrix(y)
                                   y[y > 0] <- 1
                                   lapply(setNames(nm = 1:8), function(z) {
                                     y <- apply.filter(y %>% as.data.frame(), n.time = z, exact = TRUE)
                                     y %>% colSums()
                                   }) %>% bind_rows(.id = "n")
                                 }) %>% bind_rows(.id = "rep")
                               }) %>% bind_rows(.id = "treatment") %>%
  dplyr::mutate(treatment := treatment %>%
                  gsub("^g", "G", .) %>%
                  gsub("^s", "S", .)) %>%
  pivot_longer(-c(treatment, rep, n), names_to = "time", values_to = "Counts") %>%
  ggplot(aes(time %>% factor(levels = sort(time %>% unique() %>% as.numeric())), Counts, color = n)) +
  geom_boxplot(linewidth = .2, outlier.size = .5) + facet_wrap(~treatment) + xlab("Passage")

# Taylor's plots

snps.taylor.gg <- ggarrange(TL.classic$snps.gg$gradual + ggtitle("Gradual") +
                              theme(axis.title = element_blank()) +
                              geom_text(aes(x = -.5, y = -2.8, label = paste0(round(V, 2), ", ", round(beta, 2))),
                                        color = "black", size = 2,
                                        data = TL.classic$TL$gradual$params %>%
                                          dplyr::rename(title = data)),
                            TL.classic$snps.gg$sudden + ggtitle("Sudden") +
                              theme(axis.title = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.text.y = element_blank()) +
                              geom_text(aes(x = -.5, y = -2.8, label = paste0(round(V, 2), ", ", round(beta, 2))),
                                        color = "black", size = 2,
                                        data = TL.classic$TL$sudden$params %>%
                                          dplyr::rename(title = data)),
                            nrow = 1,
                            common.legend = TRUE,
                            legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                       theme(legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1)) +
                                                       scale_color_manual(values = scales::hue_pal()(8),
                                                                          name = expression(italic(n)))),
                            legend = "top", widths = c(1, .94)) %>% annotate_figure(left = "Allele frequency standard deviation (log)", bottom = "Allele frequency mean (log)")

ggsave("4_snps_TL.pdf", snps.taylor.gg, width = 6.85, height = 4)
ggsave("4_snps_TL.tiff", snps.taylor.gg, width = 6.85, height = 4, dpi = 600)

params.location.gg <- ggarrange(facet.V.beta.gg +
                                  labs(tag = "a") +
                                  theme(legend.position = "none") +
                                  theme(plot.margin = unit(c(5.5, 13.5, 5.5, 5.5), "pt")),
                                print(location.SNPs.gg +
                                        labs(tag = "b") +
                                        theme(legend.position = "none",
                                              axis.text.y.right = element_blank(),
                                              axis.ticks.y.right = element_blank()) +
                                        ylim(0, 85) +
                                        scale_y_break(breaks = c(39, 79),
                                                      ticklabels = c(80, 85))),
                                nrow = 2,
                                common.legend = TRUE,
                                legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                           theme(legend.direction = "horizontal") +
                                                           guides(colour = guide_legend(nrow = 1)) +
                                                           scale_color_manual(values = scales::hue_pal()(8),
                                                                              name = expression(italic(n)))),
                                legend = "bottom", heights = c(.8, 1))

ggsave("5_params_TL.pdf", params.location.gg, width = 6.85, height = 4.75)
ggsave("5_params_TL.tiff", params.location.gg, width = 6.85, height = 4.75, dpi = 600)

# random fluctuations

random.TL.gg <- plot.power_law(random.TL.classic$log_log.mean_sd,
                               list(`Random 0.2` = random.walk.tb,
                                    `Random 0.05` = random.walk.tb.2,
                                    `Random 0.01` = random.walk.tb.3) %>% calc.0_rate(),
                               metadata.range = NULL, legend.title = "0 rate", l10 = TRUE) +
  xlab("Mean (log)") + ylab("Standard deviation (log)") +
  geom_smooth(method = "lm", se = FALSE, linewidth = .1) +
  geom_smooth(aes(mean, sd), method = "lm", se = FALSE, color = "black", linewidth = .1) +
  scale_color_discrete(direction = -1) +
  geom_text(aes(x = -.5, y = -2.8, label = paste0(round(V, 2), ", ", round(beta, 2))),
            color = "black", size = 2,
            data = random.TL.classic$params %>%
              dplyr::rename(title = data))

walks.gg <- list(`Random 0.2` = random.walk.tb,
                 `Random 0.05` = random.walk.tb.2,
                 `Random 0.01` = random.walk.tb.3) %>% lapply(function(x) {
                   
                   t(x) %>% as.data.frame() %>% dplyr::mutate(t = 1:8) %>%
                     pivot_longer(-t)
                   
                 }) %>% bind_rows(.id = "data") %>% ggplot(aes(t, value, group = name)) +
  geom_line(linewidth = .05) +
  facet_wrap(~data) +
  xlab("Passage") + ylab("Frequency")

# save partial supplementary figure

set.seed(101)

example.n8.distributions <- data.frame(a = rnorm(8, .01, .001),
                                       b.1 = c(rep(.1, 7), 1),
                                       b.4 = c(rep(.1, 4), rep(1, 4)),
                                       b.7 = c(.1, rep(1, 7)),
                                       
                                       
                                       c.1 = c(rep(.5, 7), 1),
                                       c.4 = c(rep(.5, 4), rep(1, 4)),
                                       c.7 = c(.5, rep(1, 7)),
                                       
                                       d.1 = c(rep(.9, 7), 1),
                                       d.4 = c(rep(.9, 4), rep(1, 4)),
                                       d.7 = c(.9, rep(1, 7))) %>% t() %>% as.data.frame()

example.n8.distributions.TL <- fit.TL(list(Distributions = example.n8.distributions), zero.rate.threshold = NULL, normalize = FALSE)

example.n8.distributions.gg <- ggarrange(example.n8.distributions.TL$log_log.mean_sd$Distributions %>%
                                           mutate(., Allele = rownames(.),
                                                  Group = rownames(.) %>% gsub("\\..*", "", .),
                                                  n1 = rownames(.) %>% gsub(".*\\.", "", .)) %>%
                                           ggplot(aes(mean, sd, color = Group)) +
                                           geom_point() +
                                           ylab("Standard deviation (log)") + xlab("Mean (log)"),
                                         example.n8.distributions %>% t() %>% as.data.frame() %>% dplyr::mutate(t = 1:8) %>%
                                           pivot_longer(-t) %>%
                                           mutate(., Allele = name,
                                                  Group = name %>% gsub("\\..*", "", .),
                                                  n1 = name %>% gsub(".*\\.", "", .)) %>%
                                           ggplot(aes(t, value, group = name, color = Group)) +
                                           geom_line(linewidth = .2) +
                                           xlab("Passage") + ylab("Frequency"),
                                         common.legend = TRUE, legend = "top")

random.gg <- ggarrange(example.n8.distributions.gg %>% annotate_figure(fig.lab = "a"),
                       walks.gg %>% annotate_figure(fig.lab = "b"),
                       ggarrange(random.TL.gg, common.legend = TRUE, legend = "bottom",
                                 legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                            theme(legend.direction = "horizontal") +
                                                            guides(colour = guide_legend(nrow = 1)) +
                                                            scale_color_manual(values = scales::hue_pal()(8),
                                                                               name = expression(italic(n))))) %>%
                         annotate_figure(fig.lab = "c"),
                       ncol = 1, heights = c(.8, .6, 1))

ggsave("s1_random.pdf", random.gg, width = 6.85, height = 7)
ggsave("s1_random.tiff", random.gg, width = 6.85, height = 7, dpi = 600)

# MANOVA

V.beta.manova <- manova(cbind(V, beta) ~ treatment * n * data, V.beta.per_0.tb %>%
                          mutate(n := as.numeric(n)) %>% dplyr::filter(treatment %in% c("gradual", "sudden")))

V.beta.manova.summary <- V.beta.manova %>% summary()

# V and beta separately

V.beta.manova.summary.aov <- V.beta.manova %>% summary.aov()

# effect size

V.beta.manova.effect_size <- effectsize::effectsize(V.beta.manova)

# lm for V and beta and contrasts

V.beta.lm.lst <- list(V = lm(V ~ treatment * n * data + beta,
                             V.beta.per_0.tb %>%
                               mutate(n := as.numeric(n)) %>%
                               dplyr::filter(treatment %in% c("gradual", "sudden"))),
                      beta = lm(beta ~ treatment * n * data + V,
                                V.beta.per_0.tb %>%
                                  mutate(n := as.numeric(n)) %>%
                                  dplyr::filter(treatment %in% c("gradual", "sudden"))))

V.beta.lm.contrasts.n <- list(V = lapply(setNames(nm = 2:8), function(x) {
  comparisons(V.beta.lm.lst$V,
              variables = list(treatment = c("sudden", "gradual")),
              newdata = datagrid(n = x))
}),
beta = lapply(setNames(nm = 2:8), function(x) {
  comparisons(V.beta.lm.lst$beta,
              variables = list(treatment = c("sudden", "gradual")),
              newdata = datagrid(n = x))
}))

# random effects

V.beta.lm.random.lst <- list(V = nlme::lme(V ~ treatment * n + beta,
                                           random = ~ 1|data,
                                           V.beta.per_0.tb %>%
                                             mutate(n := as.numeric(n)) %>%
                                             dplyr::filter(treatment %in% c("gradual", "sudden"))),
                             beta = nlme::lme(beta ~ treatment * n + V,
                                              random = ~ 1|data,
                                              V.beta.per_0.tb %>%
                                                mutate(n := as.numeric(n)) %>%
                                                dplyr::filter(treatment %in% c("gradual", "sudden"))))

# wilcoxon test

n.wilcox.V <- lapply(setNames(nm = 2:8), function(x) {
  wilcox.test(V.beta.per_0.tb[V.beta.per_0.tb$treatment == "gradual" & V.beta.per_0.tb$n == x,]$V,
              V.beta.per_0.tb[V.beta.per_0.tb$treatment == "sudden" & V.beta.per_0.tb$n == x,]$V)$p.value
}) %>% unlist() %>% p.adjust("bonferroni")

n.wilcox.beta <- lapply(setNames(nm = 2:8), function(x) {
  wilcox.test(V.beta.per_0.tb[V.beta.per_0.tb$treatment == "gradual" & V.beta.per_0.tb$n == x,]$beta,
              V.beta.per_0.tb[V.beta.per_0.tb$treatment == "sudden" & V.beta.per_0.tb$n == x,]$beta)$p.value
}) %>% unlist() %>% p.adjust("bonferroni")

# pull all alleles together

cat.mtx <- lapply(list(gradual = gradual.mtx,
                       sudden = sudden.mtx), function(x) {
                         mapply(function(y, z) {
                           rownames(y) <- paste0(z, "_", rownames(y))
                           y
                         }, x, names(x), SIMPLIFY = FALSE) %>% bind_rows()
                       })

cat.mtx.TL <- fit.TL(cat.mtx, zero.rate.threshold = NULL, normalize = FALSE)
  
# Taylor's law with mutations pulled together but comparing alleles with s < 0 and s > 0

get.mut.s <- function(x, treat, s_ht_0,
                      filter.range = 0,
                      within.range = FALSE) {
  
  if (within.range) {
    
    y <- approxwf.res[approxwf.res$mean_s < filter.range &
                        approxwf.res$mean_s > -filter.range,]
    
    sel.muts <- y[c("population", "allele")]
    sel.muts <- paste(sel.muts$population, sel.muts$allele, sep = "_")
    
    return(x[rownames(x) %in% sel.muts,])
  }
  
  y <- approxwf.res[approxwf.res$mean_s >= filter.range |
                      approxwf.res$mean_s <= -filter.range,]
  
  sel.muts <- y[y$treatment == treat &
                  (y$mean_s > 0) == s_ht_0, c("population", "allele")]
  sel.muts <- paste(sel.muts$population, sel.muts$allele, sep = "_")
  x[rownames(x) %in% sel.muts,]
}

# pull all SNPs together by s

# discern neutral alleles

# distribution of selection coefficients

approxwf.res.s_distribution.gg <- approxwf.res %>%
  mutate(Treatment = treatment %>% gsub("^g", "G", .) %>% gsub("^s", "S", .)) %>%
  ggplot(aes(mean_s, color = Treatment, fill = Treatment)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -0.16, linetype = "dashed") +
  geom_vline(xintercept = 0.16, linetype = "dashed") +
  ylab("Distribution") + xlab(expression("Mean"~italic(s)))

cat.mtx.s_16 <- list(`Gradual positive` = cat.mtx$gradual %>% get.mut.s("gradual", TRUE, filter.range = .16),
                     `Gradual negative` = cat.mtx$gradual %>% get.mut.s("gradual", FALSE, filter.range = .16),
                     `Gradual neutral` = cat.mtx$gradual %>% get.mut.s("gradual", filter.range = .16, within.range = TRUE),
                     `Sudden positive` = cat.mtx$sudden %>% get.mut.s("sudden", TRUE, filter.range = .16),
                     `Sudden negative` = cat.mtx$sudden %>% get.mut.s("sudden", FALSE, filter.range = .16),
                     `Sudden neutral` = cat.mtx$sudden %>% get.mut.s("sudden", filter.range = .16, within.range = TRUE))

cat.s_16.TL <- fit.TL(cat.mtx.s_16,
                      zero.rate.threshold = NULL, normalize = FALSE)

cat.s_16.TL.per_0 <- lapply(setNames(1:8, nm = as.character(1:8)), function(x) {
  
  # fit SNPs frequencies to Taylor's law
  
  TL <- fit.TL(cat.mtx.s_16 %>% lapply(apply.filter, x, exact = TRUE),
               zero.rate.threshold = NULL, normalize = FALSE, min.rows = 5)
  
  list(TL = TL)
})

cat.s_16.TL.gg <- plot.power_law(cat.s_16.TL$log_log.mean_sd,
                                 cat.mtx.s_16 %>%
                                   calc.0_rate() %>% lapply(function(x) (1 - x) * 8),
                                 metadata.range = NULL, legend.title = "n", l10 = TRUE,
                                 shadow = list(`Gradual positive` = cat.mtx.TL$log_log.mean_sd$gradual,
                                               `Gradual negative` = cat.mtx.TL$log_log.mean_sd$gradual,
                                               `Gradual neutral` = cat.mtx.TL$log_log.mean_sd$gradual,
                                               `Sudden positive` = cat.mtx.TL$log_log.mean_sd$sudden,
                                               `Sudden negative` = cat.mtx.TL$log_log.mean_sd$sudden,
                                               `Sudden neutral` = cat.mtx.TL$log_log.mean_sd$sudden)) +
  xlab("Allele frequency mean (log)") + ylab("Allele frequency standard deviation (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2) +
  scale_color_manual(values = scales::hue_pal()(8),
                     name = expression(italic(n))) +
  facet_wrap(~title) +
  geom_smooth(aes(mean, sd), method = "lm", formula = y ~ x, se = FALSE, linewidth = .1, color = "black") +
  geom_text(aes(x = -.5, y = -2.8, label = paste0(round(V, 2), ", ", round(beta, 2))),
            color = "black", size = 2,
            data = cat.s_16.TL$params %>%
              dplyr::rename(title = data))

cat.s_16.TL.per_0.tb <- cat.s_16.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
  bind_rows(.id = "n") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data), s = gsub(".* ", "", data))

cat.s_16.TL.params.gg <- cat.s_16.TL.per_0.tb %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(aes(shape = data)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .2)

cat.s_16.TL.per_0.V.beta.gg <- cat.s_16.TL.per_0.tb %>%
  pivot_longer(c(V, beta)) %>%
  mutate(name := name %>% factor(levels = c("V", "beta")),
         s := s %>% gsub("^n", "N", .) %>% gsub("^p", "P", .)) %>%
  ggplot(aes(s, value)) +
  facet_grid(name ~ treatment, scales = "free_y", switch = "y",
             labeller = as_labeller(c(V = "italic(V[n])", beta = "italic(\u03B2[n])", Gradual = "Gradual", Sudden = "Sudden"),  label_parsed)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(aes(color = n)) +
  theme(axis.title.y = element_blank(),
        strip.background.y = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), legend.position = "none") +
  geom_pwc(aes(group = s), method = "wilcox.test", label = "p.signif", hide.ns = "p.signif", vjust = .68, tip.length = 0)

# lm test

cat.s_16.TL.per_0.lm <- list(V = lm(V ~ s * treatment * sparsity + beta, cat.s_16.TL.per_0.tb),
                             beta = lm(beta ~ s * treatment * sparsity + V, cat.s_16.TL.per_0.tb))

# MANOVA

cat.s_16.manova <- manova(cbind(V, beta) ~ s * treatment * sparsity, cat.s_16.TL.per_0.tb)

cat.s_16.manova.summary <- cat.s_16.manova %>% summary()

# effect size

cat.s_16.manova.effect_size <- effectsize::effectsize(cat.s_16.manova)

# V and beta separately

cat.s_16.manova.summary.aov <- cat.s_16.manova %>% summary.aov()

# save fig

TL.s.gg <- ggarrange((approxwf.res.s_distribution.gg + theme(legend.position = "top")) %>%
                       annotate_figure(fig.lab = "a"),
                     ggarrange(cat.s_16.TL.gg, common.legend = TRUE,
                               legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                          theme(legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1)) +
                                                          scale_color_manual(values = scales::hue_pal()(8),
                                                                             name = expression(italic(n)))),
                               legend = "top") %>%
                       annotate_figure(fig.lab = "b"),
                     cat.s_16.TL.per_0.V.beta.gg %>% annotate_figure(fig.lab = "c"),
                     ncol = 1, heights = c(.5, 1, .75))

ggsave("6_s_fluctuation.pdf", TL.s.gg, width = 6.85, height = 9)
ggsave("6_s_fluctuation.tiff", TL.s.gg, width = 6.85, height = 9, dpi = 600)
