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
  #random.TL <- fit.TL(random.mtx,
  #                    zero.rate.threshold = NULL, normalize = FALSE,
  #                    remove.zeros = remove.zeros, sd.from.max = sd.from.max)
  
  # plot parameters
  
  snps.TL.params <- rbind(gradual.TL$params[-1] %>% cbind(data = "gradual"),
                          sudden.TL$params[-1] %>% cbind(data = "sudden"))
                          #random.TL$params[-1] %>% cbind(data = "random"))
  
  snps.TL.param.gg <- ggplot(snps.TL.params, aes(V, beta, color = data)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x)
  
  # Taylor's plot
  
  snps.TP.gg <- mapply(function(x, y, z) {
    plot.power_law(x,
                   y %>% calc.0_rate(),
                   metadata.range = NULL, legend.title = "0 rate", l10 = TRUE) +
      xlab("mean (log)") + ylab("standard deviation (log)") +
      geom_smooth(method = "lm", se = FALSE, linewidth = .2) +
      ggtitle(z)
      
  },
  list(gradual = gradual.TL$log_log.mean_sd,
       sudden = sudden.TL$log_log.mean_sd),
       #random = random.TL$log_log.mean_sd),
  list(gradual = gradual.mtx,
       sudden = sudden.mtx),
       #random = random.mtx),
  list("gradual", "sudden"),
       #"random"),
  SIMPLIFY = FALSE)
  
  list(TL = list(gradual = gradual.TL,
                 sudden = sudden.TL),
                 #random = random.TL),
       TL.gg = snps.TL.param.gg,
       snps.gg = snps.TP.gg)
}

# create random matrices

random.mtx <- lapply(setNames(nm = 1:9), function(x) {
  set.seed(x)
  random <- matrix(runif(8 * 500, min = 0, max = 1),
                   nrow = 500) %>% as.data.frame()
  random <- cbind("V0" = 0, random)
  random[random <= 0.01] <- 0
  rownames(random) <- paste0("SNP", 1:500)
  return(random)
})

# fit SNPs frequencies to Taylor's law

TL.classic <- pipe(gradual.mtx, sudden.mtx, FALSE, FALSE)
TL.classic$snps.gg <- lapply(TL.classic$snps.gg, function(x) {
  x + xlim(-3, 0) + ylim(-2.5, -.45) +
    theme(legend.position = "none",
          strip.text = element_text(size = 6)) +
    geom_smooth(aes(mean, sd), method = "lm", linetype = "dashed", se = FALSE, color = "black", linewidth = .2) +
    scale_color_discrete(direction = -1) 
})

# how beta scales with V

all.snps.V.beta.scale <- lapply(list(gradual = TL.classic$TL$gradual$params,
                                     sudden = TL.classic$TL$sudden$params),
                                function(x) {
                                  lm(beta ~ V, x) %>% summary()
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

# some random matrices

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

facet.V.beta.gg <- V.beta.per_0.tb %>% mutate(treatment := factor(treatment, levels = c("gradual", "sudden", "random"))) %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) + facet_wrap(~treatment)

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
  pivot_longer(-c(treatment, rep, n), names_to = "time", values_to = "counts") %>%
  ggplot(aes(time %>% factor(levels = sort(time %>% unique() %>% as.numeric())), counts, color = n)) +
  geom_boxplot(linewidth = .2, outlier.size = .5) + facet_wrap(~treatment) + xlab("time")

# make figure

snps.taylor.gg <- ggarrange(TL.classic$snps.gg$gradual + theme(axis.title = element_blank()),
                            TL.classic$snps.gg$sudden + theme(axis.title = element_blank(),
                                                              axis.ticks.y = element_blank(),
                                                              axis.text.y = element_blank()),
                            nrow = 1,
                            common.legend = TRUE,
                            legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                       theme(legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1))),
                            legend = "top", widths = c(1, .94)) %>% annotate_figure(left = "standard deviation (log)", bottom = "mean (log)")

ggsave("snps_TL.pdf", snps.taylor.gg, width = 6.85, height = 4)

params.location.gg <- ggarrange(facet.V.beta.gg + labs(tag = "a") + theme(legend.position = "none"),
                         location.SNPs.gg + labs(tag = "b") + theme(legend.position = "none"),
                         nrow = 2,
                         common.legend = TRUE,
                         legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                    theme(legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1))),
                         legend = "bottom")

ggsave("params_TL.pdf", params.location.gg, width = 6.85, height = 6)

# V and beta power relationship

V.beta_log.per_0.tb <- V.beta.per_0.tb %>% dplyr::mutate(`V (log)` = log10(V), `beta (log)` = log10(beta))

snps.taylor_log.gg <- V.beta_log.per_0.tb %>%
  ggplot(aes(`V (log)`, `beta (log)`, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) + facet_wrap(~treatment)

# test for the effect of treatment * zero-rate * param + data

V.lm <- lm(V ~ treatment + n:beta + data , V.beta.per_0.tb %>% mutate(n := as.numeric(n)))
beta.lm <- lm(beta ~ treatment + n:V + data, V.beta.per_0.tb %>% mutate(n := as.numeric(n)))

V.beta.manova <- manova(cbind(V, beta) ~ treatment * n * data, V.beta.per_0.tb %>% mutate(n := as.numeric(n)))

# test for n

n.wilcox.V <- lapply(setNames(nm = 2:8), function(x) {
  wilcox.test(V.beta.per_0.tb[V.beta.per_0.tb$treatment == "gradual" & V.beta.per_0.tb$n == x,]$V,
              V.beta.per_0.tb[V.beta.per_0.tb$treatment == "sudden" & V.beta.per_0.tb$n == x,]$V)$p.value
}) %>% unlist() %>% p.adjust("bonferroni")

n.wilcox.beta <- lapply(setNames(nm = 2:8), function(x) {
  wilcox.test(V.beta.per_0.tb[V.beta.per_0.tb$treatment == "gradual" & V.beta.per_0.tb$n == x,]$beta,
              V.beta.per_0.tb[V.beta.per_0.tb$treatment == "sudden" & V.beta.per_0.tb$n == x,]$beta)$p.value
}) %>% unlist() %>% p.adjust("bonferroni")

# test if the V ~ beta slope is higher depending on host replacement regime

V.beta.slope.per_0.tb <- lapply(setNames(nm = c("gradual", "sudden")), function(x) {
  lapply(setNames(nm = 2:8), function(y) {
    data <- V.beta.per_0.tb[V.beta.per_0.tb$n == y & V.beta.per_0.tb$treatment == x,]
    slope <- summary(lm(beta ~ V, data))$coefficients[2]
    data.frame(treatment = x, n = y, slope = slope)
  }) %>% bind_rows()
}) %>% bind_rows()

V.beta.slope.lm <- lm(slope ~ n + treatment, V.beta.slope.per_0.tb)

V.beta.slope.wilcox <- wilcox.test(slope ~ treatment,
                                   V.beta.slope.per_0.tb,
                                   paired = TRUE)

# test if the V ~ beta POWER relationship slope is higher depending on host replacement regime

V.beta_log.slope.per_0.tb <- lapply(setNames(nm = c("gradual", "sudden")), function(x) {
  lapply(setNames(nm = 2:8), function(y) {
    data <- V.beta_log.per_0.tb[V.beta_log.per_0.tb$n == y & V.beta_log.per_0.tb$treatment == x,]
    slope <- summary(lm(`beta (log)` ~ `V (log)`, data))$coefficients[2]
    data.frame(treatment = x, n = y, slope = slope)
  }) %>% bind_rows()
}) %>% bind_rows()

V.beta_log.slope.lm <- lm(slope ~ n + treatment, V.beta_log.slope.per_0.tb)

V.beta_log.slope.wilcox <- wilcox.test(slope ~ treatment,
                                       V.beta_log.slope.per_0.tb,
                                       paired = TRUE)

# pull all alleles together

cat.mtx <- lapply(list(gradual = gradual.mtx,
                       sudden = sudden.mtx), function(x) {
                         mapply(function(y, z) {
                           rownames(y) <- paste0(z, "_", rownames(y))
                           y
                         }, x, names(x), SIMPLIFY = FALSE) %>% bind_rows()
                       })

cat.TL <- fit.TL(cat.mtx,
                 zero.rate.threshold = NULL, normalize = FALSE)

cat.TL.per_0 <- lapply(setNames(1:8, nm = as.character(1:8)), function(x) {
  
  # fit SNPs frequencies to Taylor's law
  
  TL <- fit.TL(cat.mtx %>% lapply(apply.filter, x, exact = TRUE),
               zero.rate.threshold = NULL, normalize = FALSE,)
  
  list(TL = TL)
})

cat.TL.gg <- plot.power_law(cat.TL$log_log.mean_sd,
               cat.mtx %>% calc.0_rate(),
               metadata.range = NULL, legend.title = "0 rate", l10 = TRUE) +
  xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE)

cat.TL.params.gg <- cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n, shape = data)) + geom_point()

cat.manova <- manova(cbind(V, beta) ~ data * n, cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
                       bind_rows(.id = "n") %>% mutate(n := as.numeric(n)))

# Taylor's law with mutations pulled together but comparing alleles with s < 0 and s > 0

approxwf.res <- list(gradual = lapply(setNames(nm = c(98, 101, 102,
                                                      106, 107, 109,
                                                      112, 115, 116) %>% as.character()), function(x) {
                                                        y <- read.table(paste0("SE_results/Results/Gradual.", x, "_res.csv"),
                                                                        sep = ",", header = TRUE)
                                                        y$allele <- rownames(gradual.mtx[[as.character(x)]])
                                                        y
                                                      }) %>% bind_rows(.id = "population"),
                     sudden = lapply(setNames(nm = c(2, 5, 6,
                                                     10, 11, 13,
                                                     16, 19, 20) %>% as.character()), function(x) {
                                                       y <- read.table(paste0("SE_results/Results/Sudden.", x, "_res.csv"),
                                                                       sep = ",", header = TRUE)
                                                       y$allele <- rownames(sudden.mtx[[as.character(x)]])
                                                       y
                                                     }) %>% bind_rows(.id = "population")) %>%
  bind_rows(.id = "treatment")

get.mut.s <- function(x, treat, s_ht_0) {
  sel.muts <- approxwf.res[approxwf.res$treatment == treat &
                           (approxwf.res$mean_s > 0) == s_ht_0, c("population", "allele")]
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

cat.s.TL.gg <- plot.power_law(cat.s.TL$log_log.mean_sd,
                              cat.mtx.s %>% calc.0_rate() %>% lapply(function(x) (1 - x) * 8),
                              metadata.range = NULL, legend.title = "n", l10 = TRUE) +
  xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "black", linetype = "dashed")

cat.s.TL.params.gg <- cat.s.TL.per_0 %>% lapply(function(x) x$TL$params) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(aes(shape = data)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .2)

cat.s.TL.per_0.tb <- cat.s.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
  bind_rows(.id = "n") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data), s = gsub(".*s ", "s ", data))

# manova

cat.s.manova <- manova(cbind(V, beta) ~ s * treatment * sparsity, cat.s.TL.per_0.tb)

#ggsave("params_s_TL.pdf", cat.s.TL.gg, width = 6.85, height = 5)
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

# test with entropy (higher fluctuation on higher entropy)

# compare H = I*d**gamma vs other transformations

#list(H.n = fit.time.sd(s.mtx,return.mean = FALSE, not.fix = TRUE,use.entropy = TRUE)$params$R.squared,
#     qH.p = fit.time.sd(s.mtx,return.mean = FALSE, #next.point = FALSE,
#                        not.fix = TRUE, use.quasi.entropy = TRUE, logt = FALSE)$params$R.squared) %>%
#  bind_rows() %>% ggplot(aes(H.n, qH.p)) + geom_point() + xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) 
#
#list(H.n = fit.time.sd(s.mtx,return.mean = FALSE, not.fix = TRUE,use.entropy = TRUE)$lm.log_log.mean_sd,
#     qH.p = fit.time.sd(s.mtx,return.mean = FALSE, #next.point = FALSE,
#                        not.fix = TRUE, use.quasi.entropy = TRUE, logt = FALSE)$lm.log_log.mean_sd) %>% 
#  lapply(function(x) {
#    lapply(x, function(y) summary(y)$coefficient[8])
#  }) %>% bind_rows() %>% t() %>% as.data.frame() %>%
#  ggplot(aes(V1, V2)) + geom_point() + xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25)
#
#wilcox.test(fit.time.sd(s.mtx,return.mean = FALSE, not.fix = TRUE,use.entropy = TRUE)$params$R.squared,
#            fit.time.sd(s.mtx,return.mean = FALSE, #next.point = FALSE,
#                        not.fix = TRUE,use.quasi.entropy = TRUE, logt = FALSE)$params$R.squared,
#            paired = TRUE)

# only alleles that are fixed
# too few points

s.TL.time.point.not.fix.fix1.entropy <- fit.time.sd(cat.mtx.s %>% lapply(function(x) x[1:7]),
                                                    return.mean = FALSE, not.fix = TRUE, fix = TRUE,
                                                    use.entropy = TRUE)

s.TL.time.point.not.fix.fix1.entropy.tb <- s.TL.time.point.not.fix.fix1.entropy$log_log.mean_sd %>% bind_rows(.id = "data") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .))

s.TL.time.point.not.fix.fix1.entropy.gg <- s.TL.time.point.not.fix.fix1.entropy.tb %>%
  ggplot(aes(mean, sd)) + geom_point(shape = 20, stroke = 0, size = .75) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + geom_density_2d(linewidth = .2) +
  xlab("entropy (log)") + ylab("absolute difference from next value (log)") + facet_wrap(~treatment)

# all time points except the one that lead to loss/fixation
# add gradual vs gradual 1:7

s.TL.time.point.not.fix.entropy <- fit.time.sd(c(s.mtx,
                                                 s.mtx[grep("gradual", names(s.mtx))] %>%
                                                   setNames(nm = paste0(names(s.mtx[grep("gradual", names(s.mtx))]), " 1:7")) %>%
                                                   lapply(function(x) x[1:7])),
                                               return.mean = FALSE, not.fix = TRUE, 
                                               use.entropy = TRUE)

s.TL.time.point.not.fix.entropy.tb <- s.TL.time.point.not.fix.entropy$log_log.mean_sd %>% bind_rows(.id = "data") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .) %>% gsub(" [0-9]* 1:7", "", .))

s.TL.time.point.not.fix.entropy.params <- s.TL.time.point.not.fix.entropy$params %>%
  dplyr::mutate(treatment = gsub(" .*0 [0-9]*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .) %>% gsub(" [0-9]* 1:7", "", .),
                population = gsub(" s . 0", "", data) %>% gsub(" 1:7", "", .))

# MANOVA

gradual_vs_gradual1.7_s.TL.time.point.no.fix.entropy.I.gamma.manova <- manova(cbind(V, beta) ~ treatment * s + population,
                                                                              s.TL.time.point.not.fix.entropy.params %>%
                                                                                dplyr::filter(grepl("gradual", data)))

# t-test slopes

s.TL.time.point.not.fix.entropy.t_test.gamma <- s.TL.time.point.not.fix.entropy$lm.log_log.mean_sd %>% lapply(function(x) {
  summary(x)$coefficient[8]
}) %>% bind_rows() %>% t() %>% as.data.frame()

s.TL.time.point.not.fix.entropy.t_test.gamma$treatment <- rownames(s.TL.time.point.not.fix.entropy.t_test.gamma) %>% gsub(" .*", "", .)
s.TL.time.point.not.fix.entropy.t_test.gamma$s <- rownames(s.TL.time.point.not.fix.entropy.t_test.gamma) %>%
  gsub(" [0-9]*$", "", .) %>% gsub(".* s", "s", .) %>% gsub(" [0-9]* 1:7", "", .)
s.TL.time.point.not.fix.entropy.t_test.gamma$p.adj <- p.adjust(s.TL.time.point.not.fix.entropy.t_test.gamma$V1, "bonferroni")

s.TL.time.point.not.fix.entropy.params$slope.significant <- s.TL.time.point.not.fix.entropy.t_test.gamma$p.adj < .05

s.TL.time.point.not.fix.entropy.params.slope.sig <- s.TL.time.point.not.fix.entropy.params[s.TL.time.point.not.fix.entropy.params$slope.significant,]

# filter gradual 1:7 out 

s.TL.time.point.not.fix.entropy.tb <- s.TL.time.point.not.fix.entropy.tb[!grepl("1:7", s.TL.time.point.not.fix.entropy.tb$data),]
s.TL.time.point.not.fix.entropy.params <- s.TL.time.point.not.fix.entropy.params[!grepl("1:7", s.TL.time.point.not.fix.entropy.params$data),]

# MANOVA

s.TL.time.point.no.fix.entropy.I.gamma.manova <- manova(cbind(V, beta) ~ treatment * s + population,
                                                        s.TL.time.point.not.fix.entropy.params)

# post-hoc test

s.TL.time.point.no.fix.entropy.I.gamma.manova.lst <- list(`s > 0` = manova(cbind(V, beta) ~ treatment, s.TL.time.point.not.fix.entropy.params %>%
                                                                             filter(s == "s > 0")),
                                                          `s < 0` = manova(cbind(V, beta) ~ treatment, s.TL.time.point.not.fix.entropy.params %>%
                                                                             filter(s == "s < 0")))

s.TL.time.point.no.fix.entropy.I.gamma.lm.lst <- list(gamma = list(`s > 0` = lm(beta ~ treatment, s.TL.time.point.not.fix.entropy.params %>%
                                                                                  filter(s == "s > 0")),
                                                                   `s < 0` = lm(beta ~ treatment, s.TL.time.point.not.fix.entropy.params %>%
                                                                                  filter(s == "s < 0")),
                                                                   s = lm(beta ~ s * treatment + population, s.TL.time.point.not.fix.entropy.params)),
                                                      I = list(`s > 0` = lm(V ~ treatment + beta, s.TL.time.point.not.fix.entropy.params %>%
                                                                              filter(s == "s > 0")),
                                                               `s < 0` = lm(V ~ treatment + beta, s.TL.time.point.not.fix.entropy.params %>%
                                                                              filter(s == "s < 0")),
                                                               s = lm(V ~ s * treatment + population + beta, s.TL.time.point.not.fix.entropy.params)))

# cat

#cat.s.TL.time.point.not.fix.entropy <- fit.time.sd(c(cat.mtx.s,
#                                                     cat.mtx.s[grep("gradual", names(cat.mtx.s))] %>%
#                                                       setNames(nm = paste0(names(cat.mtx.s[grep("gradual", names(cat.mtx.s))]), " 1:7")) %>%
#                                                       lapply(function(x) x[1:7])),
#                                                   return.mean = FALSE, not.fix = TRUE, 
#                                                   use.entropy = TRUE)
#
#cat.s.TL.time.point.not.fix.entropy$params %>% ggplot(aes(beta, V, label = data)) + geom_point() + geom_label()

# plots

s.TL.time.point.not.fix.entropy.regression_plot <- s.TL.time.point.not.fix.entropy.tb %>%
  ggplot(aes(mean, sd, color = treatment, group = data)) + geom_smooth(method = "lm", se = FALSE) +
  xlab("entropy (log)") + ylab("absolute difference from next value (log)") + facet_wrap(~s)

s.TL.time.point.not.fix.entropy.regression_plot_2 <- s.TL.time.point.not.fix.entropy.tb %>%
  ggplot(aes(mean, sd, color = treatment)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .2, linetype = "dashed", color = "black") + geom_density_2d(linewidth = .2) +
  geom_point(shape = 20, stroke = 0, size = .75) +
  xlab("entropy (log)") + ylab("absolute difference from next value (log)") + facet_wrap(~data)

s.TL.time.point.not.fix.entropy.treatment.gg <- s.TL.time.point.not.fix.entropy.tb %>%
  ggplot(aes(mean, sd)) + geom_point(shape = 20, stroke = 0, size = .75) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + geom_density_2d(linewidth = .2) +
  xlab("entropy (log)") + ylab("absolute difference from next value (log)") + facet_wrap(~treatment)

s.TL.time.point.not.fix.entropy.treatment.s.gg <- s.TL.time.point.not.fix.entropy.tb %>%
  ggplot(aes(mean, sd)) + geom_point(shape = 20, stroke = 0, size = .75) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + geom_density_2d(linewidth = .2) +
  xlab("entropy (log)") + ylab("absolute difference from next value (log)") + facet_wrap(~paste0(treatment, " ", s))

s.TL.time.point.no.fix.entropy.V.beta.gg <- s.TL.time.point.not.fix.entropy.params %>% 
  ggplot(aes(log10(V), log10(beta), color = treatment)) +
  ylab("gamma (log)") + xlab("I (log)") +
  geom_point() + facet_wrap(~s)

# compare with shuffled

sudden.shuffle.mtx <- lapply(setNames(1:1000, nm = paste0("sudden shuffled rep ", 1:1000)), function(x) {
  
  set.seed(x)
  y <- s.mtx[grepl("sudden", names(s.mtx))] %>% bind_rows()
  
  y[sample(8)]
})

gradual.shuffle.mtx <- lapply(setNames(1:1000, nm = paste0("gradual shuffled rep ", 1:1000)), function(x) {
  
  set.seed(x)
  y <- s.mtx[grepl("gradual", names(s.mtx))] %>% bind_rows()
  
  y[sample(8)]
})

sudden.shuffle.mtx <- c(list(`sudden unshuffled` = s.mtx[grepl("sudden", names(s.mtx))] %>% bind_rows(),
                             `gradual unshuffled` = s.mtx[grepl("gradual", names(s.mtx))] %>% bind_rows()),
                        sudden.shuffle.mtx,
                        gradual.shuffle.mtx)

sudden.shuffle.time.point.not.fix.entropy <- fit.time.sd(sudden.shuffle.mtx, return.mean = FALSE, not.fix = TRUE, use.entropy = TRUE)

# plot parameters

sudden.shuffle.time.point.not.fix.entropy$params$shuf <- sudden.shuffle.time.point.not.fix.entropy$params$data %>%
  gsub(" rep.*", "", .) %>% gsub(".* ", "", .)

sudden.shuffle.time.point.not.fix.entropy$params$treatment <- sudden.shuffle.time.point.not.fix.entropy$params$data %>%
  gsub(" .*", "", .)
  
sudden.shuffle.time.point.not.fix.entropy.gg <- sudden.shuffle.time.point.not.fix.entropy$params %>%
  ggplot(aes(V, beta, color = shuf)) + geom_point() + theme(legend.title = element_blank()) +
  xlab("I") + ylab("gamma") + facet_wrap(~treatment)

# make figure

ggsave("Hidgamma_treatment.pdf", s.TL.time.point.not.fix.entropy.treatment.gg, width = 6.85, height = 5)
ggsave("Hidgamma_regressions.pdf", s.TL.time.point.not.fix.entropy.regression_plot, width = 6.85, height = 4.5)

# permutation entropy

perm.ent <- function(matrices,
                     n = 3) {
  
  ent.lst <- list()
  
  for (m in names(matrices)) {
    
    # extract every n point
    
    x <- matrices[[m]]
    
    tri.lst <- apply(x, 1, function(y) {
      lapply(1:(length(y) - (n - 1)), function(z) {
        tri <- y[z:(z + (n - 1))]
        
        # only analyze fluctuations
        
        if (length(tri[tri %in% c(0, 1)]) > 0) return()
        if (length(diff(tri)[diff(tri) == 0]) > 0) return()
        
        ifelse(diff(tri) > 0, "up", "down") %>% paste0(collapse = "_")
      })
    }) %>% unlist() %>% shannon.entropy.lst()
    ent.lst[[m]] <- tri.lst
  }
  return(ent.lst)
}

# permutation entropy

s.cat.perm.entropy <- perm.ent(s.mtx, 3) %>% bind_cols() %>% t() %>% as.data.frame()
colnames(s.cat.perm.entropy) <- "H"
s.cat.perm.entropy$data <- rownames(s.cat.perm.entropy) %>% gsub(".* ", "", .)
s.cat.perm.entropy$s <- rownames(s.cat.perm.entropy) %>% gsub(".*s ", "s ", .) %>% gsub("[0-9]*$", "", .)
s.cat.perm.entropy$treatment <- rownames(s.cat.perm.entropy) %>% gsub(" .*", "", .)

#s.cat.perm.entropy %>%
#  ggplot(aes(treatment, H, color = s)) + geom_boxplot()

cat.perm.entropy <- perm.ent(cat.mtx.s, 3)
s.perm.entropy <- perm.ent(s.mtx, 3) %>%
  bind_rows() %>% t() %>% as.data.frame() %>% rename(H = V1) %>% tibble::rownames_to_column("data") %>%
  mutate(treatment = data %>% gsub(" .*", "", .),
         population = data %>% gsub(" .* ", " ", .),
         s = data %>% gsub(".*s ", "s ", .) %>% gsub(" [0-9]*$", "", .))

s.perm.entropy.anova <- anova(lm(H ~ treatment * s + population, s.perm.entropy))
