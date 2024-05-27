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

# join the plots

V.beta.per_0.tb <- lapply(TL.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n")

facet.V.beta.gg <- V.beta.per_0.tb %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) + facet_wrap(~treatment)

# where are the SNPs?

location.SNPs.gg <- lapply(list(gradual = gradual.mtx,
                                sudden = sudden.mtx), function(x) {
                                  lapply(x, function(y) {
                                   y <- as.matrix(y)
                                   y[y > 0] <- 1
                                   lapply(setNames(nm = 1:7), function(z) {
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

combined.gg <- ggarrange(facet.V.beta.gg + labs(tag = "a") + theme(legend.position = "none"),
                         location.SNPs.gg + labs(tag = "b") + theme(legend.position = "none"),
                         nrow = 2,
                         common.legend = TRUE,
                         legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                    theme(legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1))),
                         legend = "bottom")

ggsave("params_TL.pdf", combined.gg, width = 6.85, height = 6)

# V and beta power relationship

V.beta_log.per_0.tb <- V.beta.per_0.tb %>% dplyr::mutate(`V (log)` = log10(V), `beta (log)` = log10(beta))

snps.taylor_log.gg <- V.beta_log.per_0.tb %>%
  ggplot(aes(`V (log)`, `beta (log)`, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) + facet_wrap(~treatment)

# test for the effect of treatment + zero-rate:param + data

V.lm <- lm(V ~ treatment + n:beta + data , V.beta.per_0.tb %>% mutate(n := as.numeric(n)))
beta.lm <- lm(beta ~ treatment + n:V + data, V.beta.per_0.tb %>% mutate(n := as.numeric(n)))

V.beta.manova <- manova(cbind(V, beta) ~ treatment * n * data, V.beta.per_0.tb %>% mutate(n := as.numeric(n)))

# test n = 8

n.8.wilcox.V <- wilcox.test(V.beta.per_0.tb[V.beta.per_0.tb$treatment == "gradual" & V.beta.per_0.tb$n == 8,]$V,
                            V.beta.per_0.tb[V.beta.per_0.tb$treatment == "sudden" & V.beta.per_0.tb$n == 8,]$V)

n.8.wilcox.beta <- wilcox.test(V.beta.per_0.tb[V.beta.per_0.tb$treatment == "gradual" & V.beta.per_0.tb$n == 8,]$beta,
                               V.beta.per_0.tb[V.beta.per_0.tb$treatment == "sudden" & V.beta.per_0.tb$n == 8,]$beta)

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

# with no zero

TL.no_zero <- pipe(gradual.mtx, sudden.mtx, remove.zeros = TRUE)
TL.no_zero$snps.gg <- lapply(TL.no_zero$snps.gg, function(x) x + xlab("mean (log)") + ylab("standard deviation (log)"))

TL.no_zero.per.zero <- pipe.per_0(n = 2:8, remove.zeros = TRUE)

V.beta.no_zero.per_0.tb <- lapply(TL.no_zero.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n")

facet.V.beta.no_zero.gg <- V.beta.no_zero.per_0.tb %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) + facet_wrap(~treatment)

snps.taylor.no_zero.gg <- ggarrange(TL.no_zero$snps.gg$gradual + theme(axis.title = element_blank()),
                                    TL.no_zero$snps.gg$sudden + theme(axis.title = element_blank(),
                                                                      axis.ticks.y = element_blank(),
                                                                      axis.text.y = element_blank()),
                                    nrow = 1,
                                    common.legend = TRUE,
                                    legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                               theme(legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1))),
                                    legend = "top", widths = c(1, .94)) %>% annotate_figure(left = "standard deviation (log)", bottom = "mean (log)")

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

cat.wilcox.V <- wilcox.test(cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
                              bind_rows(.id = "n") %>% filter(data == "gradual") %>% dplyr::select(V) %>% unlist(),
                            cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
                              bind_rows(.id = "n") %>% filter(data == "sudden") %>% dplyr::select(V) %>% unlist(),
                            paired = TRUE)

cat.wilcox.beta <- wilcox.test(cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
                                 bind_rows(.id = "n") %>% filter(data == "gradual") %>% dplyr::select(beta) %>% unlist(),
                               cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
                                 bind_rows(.id = "n") %>% filter(data == "sudden") %>% dplyr::select(beta) %>% unlist(),
                               paired = TRUE)

# pull mutations that reach fixation together by treatment

fix.mtx <- lapply(list(gradual = gradual.mtx,
                       sudden = sudden.mtx), function(x) {
                         mapply(function(y, z) {
                           y <- y %>% apply.filter(freq.thresh = 1)
                           if (nrow(y) == 0) return()
                           rownames(y) <- paste0(z, "_", rownames(y))
                           y
                         }, x, names(x), SIMPLIFY = FALSE) %>% bind_rows()
                       })

fix.TL <- fit.TL(fix.mtx,
                 zero.rate.threshold = NULL, normalize = FALSE)

fix.TL.gg <- fix.TL$log_log.mean_sd %>% bind_rows(.id = "data") %>% dplyr::mutate(mut = rownames(.)) %>%
  ggplot(aes(x = mean, y = sd, color = data, label = mut)) + geom_point() + geom_label()

fix.TL.per_0 <- lapply(setNames(1:8, nm = as.character(1:8)), function(x) {
  
  # fit SNPs frequencies to Taylor's law
  
  TL <- fit.TL(fix.mtx %>% lapply(apply.filter, x, exact = TRUE),
               zero.rate.threshold = NULL, normalize = FALSE, analyze.fluctuation.only = TRUE, min.rows = 4, sd.from.max = sd.from.max)
  
  list(TL = TL)
})

fix.TL.per_0.gg <- plot.power_law(fix.TL$log_log.mean_sd,
               fix.mtx %>% calc.0_rate(),
               metadata.range = NULL, legend.title = "0 rate", l10 = TRUE) +
  xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE)

# get positions all mutations

mut.pos <- list(gradual = gradual.mtx,
                sudden = sudden.mtx) %>% lapply(lapply, function(x) {
                  apply(x, 2, function(y) {
                    y <- y[!grepl("\\+|\\-", names(y))]
                    pos <- names(y)[ifelse(y > 0, TRUE, FALSE)] %>% gsub("\\D+", "", .) %>% as.numeric()
                    genome <- rep(0, 11685)
                    genome[pos] <- 1
                    genome
                  }) %>% t() %>% as.data.frame()
                })

h.exp <- list(gradual = fit.hurst(mut.pos$gradual, d = 1950) %>% bind_rows(.id = "time"),
              sudden = fit.hurst(mut.pos$sudden, d = 1950) %>% bind_rows(.id = "time")) %>% bind_rows(.id = "treatment")

h.exp.gg <- ggplot(h.exp, aes(time %>% factor(levels = time %>% unique()), He, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_dodge(width = .75)) +
  xlab("time point") + ylab("He") +
  geom_hline(aes(yintercept = .7), linetype = "dashed")

# test effect of time + sample + data

h.exp.lm <- lm(He ~ treatment + time + data, h.exp %>% mutate(time := time %>% as.numeric()))

# only below a threshold

mut.pos.filter <- list(gradual = gradual.mtx %>% lapply(apply.filter, freq.thresh = .3, rev = TRUE),
                       sudden = sudden.mtx %>% lapply(apply.filter, freq.thresh = .3, rev = TRUE)) %>% lapply(lapply, function(x) {
                         apply(x, 2, function(y) {
                           y <- y[!grepl("\\+|\\-", names(y))]
                           pos <- names(y)[ifelse(y > 0, TRUE, FALSE)] %>% gsub("\\D+", "", .) %>% as.numeric()
                           genome <- rep(0, 11685)
                           genome[pos] <- 1
                           genome
                         }) %>% t() %>% as.data.frame()
                       })

h.exp.filter <- list(gradual = fit.hurst(mut.pos.filter$gradual, d = 2000) %>% bind_rows(.id = "time"),
                     sudden = fit.hurst(mut.pos.filter$sudden, d = 2000) %>% bind_rows(.id = "time")) %>% bind_rows(.id = "treatment")

# NA where He could not be calculated

h.exp.filter[h.exp.filter == 0] <- NA

h.exp.filter.gg <- ggplot(h.exp.filter, aes(time %>% factor(levels = time %>% unique()), He, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_dodge(width = .75)) +
  xlab("time point") + ylab("He") +
  geom_hline(aes(yintercept = .7), linetype = "dashed")

h.exp.filter.lm <- lm(He ~ treatment + time + data, h.exp.filter %>% mutate(time := time %>% as.numeric()))

# pull mutations together by time

mut.pos.time <- mut.pos %>% lapply(function(x) {
  pos.time <- lapply(setNames(c(4, 7, 10, 13, 16, 19, 22, 25), nm = c("4", "7", "10", "13", "16", "19", "22", "25")), function(y) {
    pos.vector <- lapply(x, function(z) z[rownames(z) == y,]) %>% bind_rows() %>% colSums()
    pos.vector[pos.vector > 1] <- 1
    pos.vector
  }) %>% bind_rows(.id = "time") %>% as.data.frame()
  rownames(pos.time) <- pos.time$time
  pos.time[-1]
})

h.exp.time <- fit.hurst(mut.pos.time, d = 500) %>% bind_rows(.id = "time")

h.exp.time.gg <- ggplot(h.exp.time, aes(time %>% factor(levels = time %>% unique()), He, color = data, group = data)) +
  geom_point() +
  geom_line() +
  xlab("time point") + ylab("He") +
  geom_hline(aes(yintercept = .7), linetype = "dashed")

# pull mutations together by population

mut.pos.pop <- mut.pos %>% lapply(function(x) {
  pos.pop <- lapply(x, function(y) {
    pos.vector <- colSums(y)
    pos.vector[pos.vector > 1] <- 1
    pos.vector
    }) %>% bind_rows(.id = "rep") %>% as.data.frame()
  rownames(pos.pop) <- pos.pop$rep %>% as.character()
  pos.pop[-1]
})

h.exp.pop <- fit.hurst(mut.pos.pop, d = 550) %>% bind_rows(.id = "population")

h.exp.pop.gg <- ggplot(h.exp.pop, aes(data, He, color = data, group = data)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  xlab("") + ylab("He") +
  geom_hline(aes(yintercept = .7), linetype = "dashed")

h.exp.pop.wilcox <- wilcox.test(He ~ data,
                                h.exp.pop,
                                paired = FALSE)

# pull all mutations by treatment and time

mut.cat.pos <- list(gradual = gradual.mtx,
                           sudden = sudden.mtx) %>% lapply(lapply, function(x) {
                             apply(x, 2, function(y) {
                               y <- y[!grepl("\\+|\\-", names(y))]
                               pos <- names(y)[ifelse(y > 0, TRUE, FALSE)] %>% gsub("\\D+", "", .) %>% as.numeric()
                               genome <- rep(0, 11685)
                               genome[pos] <- 1
                               genome
                             }) %>% t() %>% as.data.frame() %>% colSums()
                           }) %>% lapply(function(x) {
                             pos.vector <- bind_rows(x) %>% colSums()
                             pos.vector[pos.vector > 1] <- 1
                             pos.vector
                           }) %>% bind_rows(.id = "treatment") %>% as.data.frame()

rownames(mut.cat.pos) <- mut.cat.pos$treatment
mut.cat.pos <- mut.cat.pos[-1]

h.exp.cat <- fit.hurst(list(data = mut.cat.pos), d = 100) %>% bind_rows(.id = "treatment")

# pull only mutations that reach a threshold by treatment and time

mut.pos.thresh <- list(gradual = gradual.mtx %>% lapply(apply.filter, freq.thresh = .1),
                       sudden = sudden.mtx %>% lapply(apply.filter, freq.thresh = .1)) %>% lapply(lapply, function(x) {
                         apply(x, 2, function(y) {
                           y <- y[!grepl("\\+|\\-", names(y))]
                           pos <- names(y)[ifelse(y > 0, TRUE, FALSE)] %>% gsub("\\D+", "", .) %>% as.numeric()
                           genome <- rep(0, 11685)
                           genome[pos] <- 1
                           genome
                         }) %>% t() %>% as.data.frame() %>% colSums()
                       }) %>% lapply(function(x) {
                         pos.vector <- bind_rows(x) %>% colSums()
                         pos.vector[pos.vector > 1] <- 1
                         pos.vector
                       }) %>% bind_rows(.id = "treatment") %>% as.data.frame()

rownames(mut.pos.thresh) <- mut.pos.thresh$treatment
mut.pos.thresh <- mut.pos.thresh[-1]

h.exp.thresh <- fit.hurst(list(data = mut.pos.thresh), d = 1000) %>% bind_rows(.id = "treatment")

# make a figure

hurst.gg <- ggarrange(h.exp.gg + labs(tag = "a"),
                      h.exp.pop.gg + theme(axis.ticks.x = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.title.y = element_blank()) + labs(tag = "b"),
                      common.legend = TRUE, widths = c(1, .5))

ggsave("hurst.pdf", hurst.gg, width = 6.85, height = 4)

# calculate SI (same as RSI but with frequencies; max_hop = ncol - 1) 

cat.SI <- calc.RSI(cat.mtx %>% lapply(apply.filter, n.time = 8), shuffle.rep = 1000)

SI.gg <- ggplot(cat.SI %>% bind_rows(.id = "treatment"), aes(RSI, RSI.boot)) +
  xlab("RSI") + ylab("PSI") +
  geom_point(aes(color = ifelse(p.adjust < .05, "a", "b")), size = 1) +
  facet_wrap(~treatment) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("tomato", "grey60"))

# Taylor's law with mutations pulled together but comparing alleles with s < 0 and s > 0

WF_ABC.res <- list(gradual = lapply(setNames(nm = c(98, 101, 102,
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
  sel.muts <- WF_ABC.res[WF_ABC.res$treatment == treat &
                           (WF_ABC.res$mean_s > 0) == s_ht_0, c("population", "allele")]
  sel.muts <- paste(sel.muts$population, sel.muts$allele, sep = "_")
  x[rownames(x) %in% sel.muts,]
}

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
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2)

cat.s.TL.params.gg <- cat.s.TL.per_0 %>% lapply(function(x) x$TL$params) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(aes(shape = data)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .2)

cat.s.TL.per_0.tb <- cat.s.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
  bind_rows(.id = "n") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data), s = gsub(".*s ", "s ", data))

cat.s.V.lm <- lm(V ~ s + treatment + sparsity,
                 cat.s.TL.per_0.tb)

cat.s.beta.lm <- lm(beta ~ s + treatment + sparsity,
                    cat.s.TL.per_0.tb)

# take into account how one parameter influences another

cat.s.V.acc.beta.lm <- lm(V ~ s + treatment + sparsity + beta,
                          cat.s.TL.per_0.tb)

cat.s.beta.acc.beta.lm <- lm(beta ~ s + treatment + sparsity + V,
                             cat.s.TL.per_0.tb)

# manova

cat.s.manova <- manova(cbind(V, beta) ~ s * treatment * sparsity, cat.s.TL.per_0.tb)

ggsave("params_s_TL.pdf", cat.s.TL.gg, width = 6.85, height = 5)
ggsave("params_s.pdf", cat.s.TL.params.gg, width = 6.85, height = 5)

# investigate systemic behavior in fluctuations

# by treatment

cat.TL.time.point <- fit.time.sd(cat.mtx, return.mean = FALSE)

cat.TL.time.point.gg <- plot.power_law(cat.TL.time.point$log_log.mean_sd,
                                       metadata.range = NULL, l10 = FALSE) +
  xlab("value (log)") + ylab("absolute difference from next value (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2)# + geom_density_2d()

# do not include differences that lead to elimination/fixation

cat.TL.time.point.not.fix <- fit.time.sd(cat.mtx, return.mean = FALSE, not.fix = TRUE)

cat.TL.time.point.not.fix.gg <- plot.power_law(cat.TL.time.point.not.fix$log_log.mean_sd,
                                       metadata.range = NULL, l10 = TRUE) +
  xlab("value (log)") + ylab("absolute difference from next value (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + geom_density_2d(linewidth = .2)

# by allele

cat.TL.time.point.log_log.allele <- cat.TL.time.point$log_log.mean_sd %>%
  bind_rows(.id = "data") %>% dplyr::rename(allele = name)

cat.TL.time.point.by.gene.gg <- cat.TL.time.point.log_log.allele %>%
  ggplot(aes(mean, sd, group = allele)) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) +
  facet_wrap(~data) +
  xlab("value (log)") + ylab("absolute difference from next value (log)")

# lm by allele

cat.TL.time.point.allele.params <- cat.TL.time.point.log_log.allele$allele %>% unique() %>%
  setNames(nm = .) %>% lapply(function(x) {
    y <- cat.TL.time.point.log_log.allele[cat.TL.time.point.log_log.allele$allele == x,]
    
    if(nrow(y) < 5) return()
    
    fit.summary <- lm(sd ~ mean, y) %>% summary()
    list(allele = x,
         V = fit.summary$coefficients[1] %>% exp(),
         beta = fit.summary$coefficients[2],
         Rsquared = fit.summary$r.squared,
         `mean (log)` = exp(y$mean) %>% mean() %>% log10(),
         treatment = y$data %>% unique())
  }) %>% bind_rows()

cat.TL.time.point.allele.beta.gg <- cat.TL.time.point.allele.params %>% ggplot(aes(treatment, beta)) +
  geom_boxplot() 

cat.TL.time.point.allele.beta.lm <- lm(beta ~ treatment,
                                       cat.TL.time.point.allele.params %>% separate(allele, c("data", "allele"), "_"))

cat.TL.time.point.allele.params.manova <- manova(cbind(V, beta) ~ data + allele,
                                                 cat.TL.time.point.allele.params %>% separate(allele, c("data", "allele"), "_"))

# by s/treatment

cat.s.TL.time.point <- fit.time.sd(cat.mtx.s, return.mean = FALSE)

cat.s.TL.time.point.gg <- plot.power_law(cat.s.TL.time.point$log_log.mean_sd,
                                         metadata.range = NULL, l10 = FALSE) +
  xlab("value (log)") + ylab("absolute difference from next value (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2)

# by dataset/treatment

TL.time.point <- lapply(list(gradual = gradual.mtx,
                             sudden = sudden.mtx), fit.time.sd, return.mean = FALSE)

TL.time.point.V.beta <- lapply(TL.time.point, function(x) x$params) %>%
  bind_rows(.id = "treatment")

TL.time.point.V.beta.gg <- ggplot(TL.time.point.V.beta, aes(V, beta, color = treatment)) + geom_point()

TL.time.point.V.beta.manova <- manova(cbind(V, beta) ~ treatment, TL.time.point.V.beta)

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

s.TL.time.point <- fit.time.sd(s.mtx, return.mean = FALSE)

s.TL.time.point.regression_plot <- s.TL.time.point$log_log.mean_sd %>% bind_rows(.id = "data") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .)) %>%
  ggplot(aes(mean, sd, color = treatment, linetype = s, group = data)) + geom_smooth(method = "lm", se = FALSE) +
  xlab("value (log)") + ylab("absolute difference from next value (log)")

s.TL.time.point.params <- s.TL.time.point$params %>%
  dplyr::mutate(treatment = gsub(" .*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .))

s.TL.time.point.V.beta.gg <- ggplot(s.TL.time.point.params, aes(V, beta, color = treatment, shape = s)) + geom_point()

s.TL.time.point.V.beta.manova <- manova(cbind(V, beta) ~ treatment * s, s.TL.time.point.params)

# do not include differences that lead to elimination/fixation

s.TL.time.point.not.fix <- fit.time.sd(s.mtx, return.mean = FALSE, not.fix = TRUE)

s.TL.time.point.not.fix.regression_plot <- s.TL.time.point.not.fix$log_log.mean_sd %>% bind_rows(.id = "data") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .)) %>%
  ggplot(aes(mean, sd, color = treatment, linetype = s, group = data)) + geom_smooth(method = "lm", se = FALSE) +
  xlab("value (log)") + ylab("absolute difference from next value (log)")

s.TL.time.point.not.fix.params <- s.TL.time.point.not.fix$params %>%
  dplyr::mutate(treatment = gsub(" .*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .))

s.TL.time.point.no.fix.V.beta.gg <- ggplot(s.TL.time.point.not.fix.params, aes(V, beta, color = treatment, shape = s)) + geom_point()

s.TL.time.point.no.fix.V.beta.manova <- manova(cbind(V, beta) ~ treatment * s, s.TL.time.point.not.fix.params)

# post-hoc test



# test with entropy (higher fluctuation on higher entropy)

s.TL.time.point.not.fix.entropy <- fit.time.sd(s.mtx, return.mean = FALSE, not.fix = TRUE, use.entropy = TRUE)

s.TL.time.point.not.fix.entropy.tb <- s.TL.time.point.not.fix.entropy$log_log.mean_sd %>% bind_rows(.id = "data") %>%
  dplyr::mutate(treatment = gsub(" .*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .))

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

s.TL.time.point.not.fix.entropy.params <- s.TL.time.point.not.fix.entropy$params %>%
  dplyr::mutate(treatment = gsub(" .*", "", data),
                s = gsub(".*s ", "s ", data) %>% gsub(" [0-9]*$", "", .),
                population = gsub(" s . 0", "", data))
 
s.TL.time.point.no.fix.entropy.V.beta.gg <- ggplot(s.TL.time.point.not.fix.entropy.params, aes(log(V), log(beta), color = treatment, shape = s)) + geom_point()

s.TL.time.point.no.fix.entropy.V.beta.manova <- manova(cbind(V, beta) ~ treatment * s + population, s.TL.time.point.not.fix.entropy.params)

# post-hoc test

s.TL.time.point.no.fix.entropy.V.beta.manova.lst <- list(`s > 0` = manova(cbind(V, beta) ~ treatment, s.TL.time.point.not.fix.entropy.params %>%
                                                                            filter(s == "s > 0")),
                                                         `s < 0` = manova(cbind(V, beta) ~ treatment, s.TL.time.point.not.fix.entropy.params %>%
                                                                            filter(s == "s < 0")))

# permutation entropy

perm.ent <- function(matrices) {
  
  for (m in names(matrices)) {
    
    x <- matrices[[m]]
    
    d.x <- apply(x, 1, diff) %>% t() %>% as.data.frame()
    
    perm <- d.x %>% mutate()
    
    print(perm)
    
    # ignore permutation where one element doesn't change
  }
  
}



