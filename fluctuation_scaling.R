source("source_me.R")

# pipe function

pipe <- function(gradual.mtx,
                 sudden.mtx,
                 remove.zeros = FALSE,
                 sd.from.max = FALSE) {
  
  gradual.TL <- fit.TL(gradual.mtx,
                       zero.rate.threshold = NULL, normalize = FALSE,
                       remove.zeros = remove.zeros, sd.from.max = sd.from.max)
  sudden.TL <- fit.TL(sudden.mtx,
                      zero.rate.threshold = NULL, normalize = FALSE,
                      remove.zeros = remove.zeros, sd.from.max = sd.from.max)
  
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

TL.classic <- pipe(gradual.mtx, sudden.mtx, FALSE, FALSE)
TL.classic$snps.gg <- lapply(TL.classic$snps.gg, function(x) {
  x + xlim(-3, 0) + ylim(-2.5, -.45) +
    theme(legend.position = "none",
          strip.text = element_text(size = 6)) +
    geom_smooth(aes(mean, sd), method = "lm", se = FALSE, color = "black", linewidth = .1) +
    scale_color_discrete(direction = -1) 
})

# Taylor's law per 0 rate

pipe.per_0 <- function(n = 1:8,
                       remove.zeros = FALSE,
                       min.rows = 5,
                       sd.from.max = FALSE) {
  
  lapply(setNames(n, nm = as.character(n)), function(x) {
    
    # fit SNPs frequencies to Taylor's law
    
    gradual.TL <- fit.TL(gradual.mtx %>% lapply(apply.filter, x, exact = TRUE),
                         zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = remove.zeros, min.rows = min.rows, sd.from.max = sd.from.max)
    sudden.TL <- fit.TL(sudden.mtx %>% lapply(apply.filter, x, exact = TRUE),
                        zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = remove.zeros, min.rows = min.rows, sd.from.max = sd.from.max)
    
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

random.TL.classic <- fit.TL(list(random.2 = random.walk.tb,
                                 random.05 = random.walk.tb.2,
                                 random.01 = random.walk.tb.3),
                            zero.rate.threshold = NULL, normalize = FALSE,
                            remove.zeros = FALSE)

pipe.per_0.random <- function(n = 1:8,
                              remove.zeros = FALSE,
                              min.rows = 5,
                              sd.from.max = FALSE) {
  
  lapply(setNames(n, nm = as.character(n)), function(x) {
    
    # fit SNPs frequencies to Taylor's law
    
    fit.TL(list(random.2 = random.walk.tb,
                random.05 = random.walk.tb.2,
                random.01 = random.walk.tb.3) %>% lapply(apply.filter, x, exact = TRUE),
           zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = remove.zeros, min.rows = min.rows, sd.from.max = sd.from.max)
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
  ylab(expression(italic("\u03B2"))) + xlab(expression(italic("V")))

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

# make figure

snps.taylor.gg <- ggarrange(TL.classic$snps.gg$gradual + ggtitle("Gradual") +
                              theme(axis.title = element_blank()) ,
                            TL.classic$snps.gg$sudden + ggtitle("Sudden") +
                              theme(axis.title = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.text.y = element_blank()),
                            nrow = 1,
                            common.legend = TRUE,
                            legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                       theme(legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1)) +
                                                       scale_color_manual(values = scales::hue_pal()(8),
                                                                          name = expression(italic(n)))),
                            legend = "top", widths = c(1, .94)) %>% annotate_figure(left = "Standard deviation (log)", bottom = "Mean (log)")

ggsave("3_snps_TL.pdf", snps.taylor.gg, width = 6.85, height = 4)
ggsave("3_snps_TL.png", snps.taylor.gg, width = 6.85, height = 4, dpi = 600)

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

ggsave("4_params_TL.pdf", params.location.gg, width = 6.85, height = 4.75)
ggsave("4_params_TL.png", params.location.gg, width = 6.85, height = 4.75, dpi = 600)

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

# Taylor's law with mutations pulled together but comparing alleles with s < 0 and s > 0

get.mut.s <- function(x, treat, s_ht_0, filter.range = 0) {
  y <- approxwf.res[approxwf.res$mean_s >= filter.range |
                      approxwf.res$mean_s <= -filter.range,]
  sel.muts <- y[y$treatment == treat &
                  (y$mean_s > 0) == s_ht_0, c("population", "allele")]
  sel.muts <- paste(sel.muts$population, sel.muts$allele, sep = "_")
  x[rownames(x) %in% sel.muts,]
}

# pull all SNPs together by s

cat.mtx.s <- list(`gradual s > 0` = cat.mtx$gradual %>% get.mut.s("gradual", TRUE),
                  `gradual s < 0` = cat.mtx$gradual %>% get.mut.s("gradual", FALSE),
                  `sudden s > 0` = cat.mtx$sudden %>% get.mut.s("sudden", TRUE),
                  `sudden s < 0` = cat.mtx$sudden %>% get.mut.s("sudden", FALSE))

cat.s.TL <- fit.TL(cat.mtx.s,
                   zero.rate.threshold = NULL, normalize = FALSE)

cat.s.TL.per_0 <- lapply(setNames(1:8, nm = as.character(1:8)), function(x) {
  
  # fit SNPs frequencies to Taylor's law
  
  TL <- fit.TL(cat.mtx.s %>% lapply(apply.filter, x, exact = TRUE),
               zero.rate.threshold = NULL, normalize = FALSE, min.rows = 5)
  
  list(TL = TL)
})

cat.s.TL.gg <- plot.power_law(cat.s.TL$log_log.mean_sd %>% setNames(nm = c("Gradual~italic(s)~'>'~0",
                                                                           "Gradual~italic(s)~'<'~0",
                                                                           "Sudden~italic(s)~'>'~0",
                                                                           "Sudden~italic(s)~'<'~0")),
                              cat.mtx.s %>% setNames(nm = c("Gradual~italic(s)~'>'~0",
                                                            "Gradual~italic(s)~'<'~0",
                                                            "Sudden~italic(s)~'>'~0",
                                                            "Sudden~italic(s)~'<'~0")) %>%
                                calc.0_rate() %>% lapply(function(x) (1 - x) * 8),
                              metadata.range = NULL, legend.title = "n", l10 = TRUE) +
  xlab("Mean (log)") + ylab("Standard deviation (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2) +
  scale_color_manual(values = scales::hue_pal()(8),
                     name = expression(italic(n))) +
  facet_wrap(~title, labeller = label_parsed) +
  geom_smooth(aes(mean, sd), method = "lm", formula = y ~ x, se = FALSE, linewidth = .1, color = "black")

cat.s.TL.params.gg <- cat.s.TL.per_0 %>% lapply(function(x) x$TL$params) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(aes(shape = data)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .2)

cat.s.TL.per_0.tb <- cat.s.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
  bind_rows(.id = "n") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data), s = gsub(".*s ", "s ", data))

# MANOVA

cat.s.manova <- manova(cbind(V, beta) ~ s * treatment * sparsity, cat.s.TL.per_0.tb)

cat.s.manova.summary <- cat.s.manova %>% summary()

# effect size

cat.s.manova.effect_size <- effectsize::effectsize(cat.s.manova)

# V and beta separately

cat.s.manova.summary.aov <- cat.s.manova %>% summary.aov()

# save fig

ggsave("s1_cats.pdf", cat.s.TL.gg, width = 6.85, height = 5)
ggsave("s1_cats.png", cat.s.TL.gg, width = 6.85, height = 5, dpi = 600)

ggsave("params_s.pdf", cat.s.TL.params.gg, width = 6.85, height = 5)

# investigate systemic behavior in fluctuations

# by dataset/treatment/s

s.mtx <- c(mapply(function(x, y) {
  rownames(x) <- paste(y, rownames(x), sep = "_")
  x %>% get.mut.s("gradual", TRUE)
}, gradual.mtx %>% setNames(paste0("gradual s > 0 ", names(gradual.mtx))), names(gradual.mtx), SIMPLIFY = FALSE),

mapply(function(x, y) {
  rownames(x) <- paste(y, rownames(x), sep = "_")
  x %>% get.mut.s("gradual", FALSE)
}, gradual.mtx %>% setNames(paste0("gradual s < 0 ", names(gradual.mtx))), names(gradual.mtx), SIMPLIFY = FALSE),

mapply(function(x, y) {
  rownames(x) <- paste(y, rownames(x), sep = "_")
  x %>% get.mut.s("sudden", TRUE)
}, sudden.mtx %>% setNames(paste0("sudden s > 0 ", names(sudden.mtx))), names(sudden.mtx), SIMPLIFY = FALSE),

mapply(function(x, y) {
  rownames(x) <- paste(y, rownames(x), sep = "_")
  x %>% get.mut.s("sudden", FALSE)
}, sudden.mtx %>% setNames(paste0("sudden s < 0 ", names(sudden.mtx))), names(sudden.mtx), SIMPLIFY = FALSE))

# all time points except the ones that lead to loss/fixation

s.TL.time.point.not.fix.entropy <- fit.time.sd(s.mtx,
                                               return.mean = FALSE, not.fix = TRUE, 
                                               use.entropy = TRUE)

s.TL.time.point.not.fix.entropy.tb <- s.TL.time.point.not.fix.entropy$log_log.mean_sd %>% bind_rows(.id = "data") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .))

s.TL.time.point.not.fix.entropy.params <- s.TL.time.point.not.fix.entropy$params %>%
  dplyr::mutate(treatment = gsub(" .*0 [0-9]*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .),
                population = gsub(" s . 0", "", data) %>% gsub(" 1:7", "", .))

# t-test slopes

s.TL.time.point.not.fix.entropy.t_test.gamma <- s.TL.time.point.not.fix.entropy$lm.log_log.mean_sd %>% lapply(function(x) {
  summary(x)$coefficient[8]
}) %>% bind_rows() %>% t() %>% as.data.frame()

s.TL.time.point.not.fix.entropy.t_test.gamma$treatment <- rownames(s.TL.time.point.not.fix.entropy.t_test.gamma) %>% gsub(" .*", "", .)
s.TL.time.point.not.fix.entropy.t_test.gamma$s <- rownames(s.TL.time.point.not.fix.entropy.t_test.gamma) %>%
  gsub(" [0-9]*$", "", .) %>% gsub(".* s", "s", .)
s.TL.time.point.not.fix.entropy.t_test.gamma$p.adj <- p.adjust(s.TL.time.point.not.fix.entropy.t_test.gamma$V1, "bonferroni")

s.TL.time.point.not.fix.entropy.params$slope.significant <- s.TL.time.point.not.fix.entropy.t_test.gamma$p.adj < .05

s.TL.time.point.not.fix.entropy.params.slope.sig <- s.TL.time.point.not.fix.entropy.params[s.TL.time.point.not.fix.entropy.params$slope.significant,]

# MANOVA
# remove outlier?

#s.TL.time.point.no.fix.entropy.I.gamma.manova <- manova(cbind(V, beta) ~ treatment * s + population,
#                                                        s.TL.time.point.not.fix.entropy.params %>% dplyr::filter(data != "gradual s > 0 116"))

s.TL.time.point.no.fix.entropy.I.gamma.manova <- manova(cbind(V, beta) ~ treatment * s + population,
                                                        s.TL.time.point.not.fix.entropy.params)

s.TL.time.point.no.fix.entropy.I.gamma.manova.summary <- s.TL.time.point.no.fix.entropy.I.gamma.manova %>% summary()

# effect size

s.TL.time.point.no.fix.entropy.I.gamma.manova.effect_size <- effectsize::effectsize(s.TL.time.point.no.fix.entropy.I.gamma.manova)

# gamma and I separately 

s.TL.time.point.no.fix.entropy.I.gamma.manova.summary.aov <- s.TL.time.point.no.fix.entropy.I.gamma.manova %>% summary.aov()

# post-hoc test

s.TL.time.point.no.fix.entropy.I.gamma.lm.lst <- list(gamma = lm(beta ~ s * treatment + population + V, s.TL.time.point.not.fix.entropy.params),
                                                      I = lm(V ~ s * treatment + population + beta, s.TL.time.point.not.fix.entropy.params))

# mixed model with population as random effect

#s.TL.time.point.no.fix.entropy.I.gamma.lm.lst <- list(gamma = lme4::lmer(beta ~ s * treatment + V + (1|population), s.TL.time.point.not.fix.entropy.params),
#                                                      I = lme4::lmer(V ~ s * treatment + beta + (1|population), s.TL.time.point.not.fix.entropy.params))

# contrasts

s.TL.time.point.no.fix.entropy.I.gamma.lm.contrasts <- list(gamma = list(`s < 0, gradual - sudden` = comparisons(s.TL.time.point.no.fix.entropy.I.gamma.lm.lst$gamma,
                                                                                                                                  variables = list(treatment = c("sudden", "gradual")),
                                                                                                                                  newdata = datagrid(s = "s < 0")),
                                                                         `s > 0, gradual - sudden` = comparisons(s.TL.time.point.no.fix.entropy.I.gamma.lm.lst$gamma,
                                                                                                                                  variables = list(treatment = c("sudden", "gradual")),
                                                                                                                                  newdata = datagrid(s = "s > 0")),
                                                                         `gradual, s < 0 - s > 0` = comparisons(s.TL.time.point.no.fix.entropy.I.gamma.lm.lst$gamma,
                                                                                                                                 variables = list(s = c("s > 0", "s < 0")),
                                                                                                                                 newdata = datagrid(treatment = "gradual")),
                                                                         `sudden, s < 0 - s > 0` = comparisons(s.TL.time.point.no.fix.entropy.I.gamma.lm.lst$gamma,
                                                                                                                                variables = list(s = c("s > 0", "s < 0")),
                                                                                                                                newdata = datagrid(treatment = "sudden"))),
                                                            I = list(`s < 0, gradual - sudden` = comparisons(s.TL.time.point.no.fix.entropy.I.gamma.lm.lst$I,
                                                                                                                              variables = list(treatment = c("sudden", "gradual")),
                                                                                                                              newdata = datagrid(s = "s < 0")),
                                                                     `s > 0, gradual - sudden` = comparisons(s.TL.time.point.no.fix.entropy.I.gamma.lm.lst$I,
                                                                                                                              variables = list(treatment = c("sudden", "gradual")),
                                                                                                                              newdata = datagrid(s = "s > 0")),
                                                                     `gradual, s < 0 - s > 0` = comparisons(s.TL.time.point.no.fix.entropy.I.gamma.lm.lst$I,
                                                                                                                             variables = list(s = c("s > 0", "s < 0")),
                                                                                                                             newdata = datagrid(treatment = "gradual")),
                                                                     `sudden, s < 0 - s > 0` = comparisons(s.TL.time.point.no.fix.entropy.I.gamma.lm.lst$I,
                                                                                                                           variables = list(s = c("s > 0", "s < 0")),
                                                                                                                           newdata = datagrid(treatment = "sudden"))))

# plots

s.TL.time.point.not.fix.entropy.regression_plot <- s.TL.time.point.not.fix.entropy.tb %>%
  dplyr::mutate(treatment := treatment %>%
                  gsub("^g", "G", .) %>%
                  gsub("^s", "S", .)) %>%
  ggplot(aes(mean, sd, color = s, group = data)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .2) +
  xlab(expression(paste(italic(S), " (log)"))) + ylab(expression(paste(italic(d), " (log)"))) +
  facet_wrap(~treatment) +
  scale_color_manual(values = scales::hue_pal()(2), 
                     labels = expression(paste(italic(s), " < 0"),
                                         paste(italic(s), " > 0")),
                     name = expression(italic(s)))

s.TL.time.point.not.fix.entropy.treatment.gg <- s.TL.time.point.not.fix.entropy.tb %>%
  dplyr::mutate(treatment := treatment %>%
                  gsub("^g", "G", .) %>%
                  gsub("^s", "S", .)) %>%
  ggplot(aes(mean, sd)) + geom_point(shape = 20, stroke = 0, size = .75) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + geom_density_2d(linewidth = .2) +
  xlab(expression(paste(italic(S), " (log)"))) + ylab(expression(paste(italic(d), " (log)"))) +
  facet_wrap(~treatment)

s.TL.time.point.not.fix.entropy.treatment.s.gg <- s.TL.time.point.not.fix.entropy.tb %>%
  ggplot(aes(mean, sd)) + geom_point(shape = 20, stroke = 0, size = .75) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + geom_density_2d(linewidth = .2) +
  xlab("entropy (log)") + ylab("absolute difference from next value (log)") + facet_wrap(~paste0(treatment, " ", s))

s.TL.time.point.no.fix.entropy.I.gamma.log.gg <- s.TL.time.point.not.fix.entropy.params %>% 
  ggplot(aes(log10(V), log10(beta), color = paste0(treatment, s))) +
  ylab("gamma (log)") + xlab("I (log)") +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + theme(legend.title = element_blank())

# no outlier

s.TL.time.point.no.fix.entropy.I.gamma.gg <- s.TL.time.point.not.fix.entropy.params %>%
  dplyr::filter(data != "gradual s > 0 116") %>%
  ggplot(aes(V, beta, color = paste0(treatment, s))) +
  ylab("gamma") + xlab("I") +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + theme(legend.title = element_blank())

s.TL.time.point.no.fix.entropy.I.gamma.boxplot.gg <- s.TL.time.point.not.fix.entropy.params %>% 
  dplyr::filter(data != "gradual s > 0 116") %>%
  rename(I = V, gamma = beta) %>%
  dplyr::select(-c(R.squared, n.cells, n.alleles, model, slope.significant, sparsity, population)) %>%
  pivot_longer(c(I, gamma)) %>%
  ggplot(aes(s, value)) + geom_boxplot() + facet_grid(rows = vars(treatment), cols = vars(name), scales = "free")

# compare with shuffled

#sudden.shuffle.mtx <- lapply(setNames(1:1000, nm = paste0("sudden shuffled rep ", 1:1000)), function(x) {
#  
#  set.seed(x)
#  y <- s.mtx[grepl("sudden", names(s.mtx))] %>% bind_rows()
#  
#  y[sample(8)]
#})
#
#gradual.shuffle.mtx <- lapply(setNames(1:1000, nm = paste0("gradual shuffled rep ", 1:1000)), function(x) {
#  
#  set.seed(x)
#  y <- s.mtx[grepl("gradual", names(s.mtx))] %>% bind_rows()
#  
#  y[sample(8)]
#})
#
#sudden.shuffle.mtx <- c(list(`sudden unshuffled` = s.mtx[grepl("sudden", names(s.mtx))] %>% bind_rows(),
#                             `gradual unshuffled` = s.mtx[grepl("gradual", names(s.mtx))] %>% bind_rows()),
#                        sudden.shuffle.mtx,
#                        gradual.shuffle.mtx)
#
#sudden.shuffle.time.point.not.fix.entropy <- fit.time.sd(sudden.shuffle.mtx, return.mean = FALSE, not.fix = TRUE, use.entropy = TRUE)

# plot parameters

sudden.shuffle.time.point.not.fix.entropy$params$shuf <- sudden.shuffle.time.point.not.fix.entropy$params$data %>%
  gsub(" rep.*", "", .) %>% gsub(".* ", "", .)

sudden.shuffle.time.point.not.fix.entropy$params$treatment <- sudden.shuffle.time.point.not.fix.entropy$params$data %>%
  gsub(" .*", "", .)

sudden.shuffle.time.point.not.fix.entropy.gg <- sudden.shuffle.time.point.not.fix.entropy$params %>%
  ggplot(aes(V, beta, color = shuf)) + geom_point() + theme(legend.title = element_blank()) +
  xlab("I") + ylab("gamma") + facet_wrap(~treatment)

# save figs

hid.gg <- ggarrange(s.TL.time.point.not.fix.entropy.treatment.gg +
                      labs(tag = "a"),
                      #theme(axis.ticks.x = element_blank(),
                      #      axis.title.x = element_blank(),
                      #      axis.text.x = element_blank()),
                    s.TL.time.point.not.fix.entropy.regression_plot +
                      theme(legend.position = "bottom") +
                      labs(tag = "b"),
                    ncol = 1,
                    heights = c(.9, 1))

ggsave("5_Hid.pdf", hid.gg,  width = 6.85, height = 6)
ggsave("5_Hid.png", hid.gg,  width = 6.85, height = 6, dpi = 600)

ggsave("Hidgamma_treatment.pdf", s.TL.time.point.not.fix.entropy.treatment.gg, width = 6.85, height = 5)
ggsave("Hidgamma_regressions.pdf", s.TL.time.point.not.fix.entropy.regression_plot, width = 6.85, height = 4.5)
