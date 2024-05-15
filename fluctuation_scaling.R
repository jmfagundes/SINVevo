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
                       min.rows = 4,
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

# to check samples

lapply(TL.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n)) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) +
  geom_point(size = .5) + facet_wrap(~treatment) + geom_label(aes(label = data))

joint.gg <- lapply(TL.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n, linetype = treatment)) +
  geom_smooth(method = "lm", se = FALSE)

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

# test for the effect of treatment + zero-rate + data

beta.lm <- lm(beta ~ treatment + n + data, V.beta.per_0.tb %>% mutate(n := as.numeric(n)))
V.lm <- lm(V ~ treatment + n + data, V.beta.per_0.tb %>% mutate(n := as.numeric(n)))

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

# Taylor's law from the max

TL.max <- pipe(gradual.mtx, sudden.mtx, FALSE, TRUE)
TL.max$snps.gg <- lapply(TL.max$snps.gg, function(x) {
  x + xlab("max (log)") + ylab("standard deviation of the max (log)") +
    xlim(-2.2, 0) + ylim(-2.2, 0) +
    theme(legend.position = "none") +
    geom_smooth(aes(mean, sd), method = "lm", linetype = "dashed", se = FALSE, color = "black", linewidth = .2) +
    scale_color_discrete(direction = -1)
})

# figures

snps.taylor.max.gg <- ggarrange(TL.max$snps.gg$gradual + theme(axis.title = element_blank()),
                                TL.max$snps.gg$sudden + theme(axis.title = element_blank(),
                                                              axis.ticks.y = element_blank(),
                                                              axis.text.y = element_blank()),
                                nrow = 1,
                                common.legend = TRUE,
                                legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                           theme(legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1))),
                                legend = "top", widths = c(1, .94)) %>% annotate_figure(left = "standard deviation of the max (log)", bottom = "max (log)")

TL.max.per.zero <- pipe.per_0(sd.from.max = TRUE)

max.facet.V.beta.gg <- lapply(TL.max.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) + facet_wrap(~treatment)

# to check samples

lapply(TL.max.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n)) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) +
  geom_point(size = .5) + facet_wrap(~treatment) + geom_label(aes(label = data))

max.joint.gg <- lapply(TL.max.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n, linetype = treatment)) +
  geom_smooth(method = "lm", se = FALSE)

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

TL.max.no_zero <- pipe(gradual.mtx, sudden.mtx, TRUE, TRUE)
TL.max.no_zero$snps.gg <- lapply(TL.max.no_zero$snps.gg, function(x) x + xlab("max (log)") + ylab("standard deviation of the max (log)"))

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

plot.power_law(fix.TL$log_log.mean_sd,
               fix.mtx %>% calc.0_rate(),
               metadata.range = NULL, legend.title = "0 rate", l10 = TRUE) +
  xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE)

# rank stability, all SNPs

#all.snps.rank <- lapply(list(gradual = gradual.mtx,
#                             sudden = sudden.mtx),
#                        function(x) process.rank(matrices = x, method = "random"))
#
#all.snps.RSI <- all.snps.rank %>% lapply(calc.RSI, shuffle.rep = 1000)
#
#all.snps.RSI.gg <- mapply(function(x, y) {
#  
#  x <- x %>% bind_rows(.id = "replicate")
#  
#  ggplot(x, aes(RSI, RSI.boot)) +
#    xlab("RSI") + ylab("PSI") +
#    geom_point(aes(color = ifelse(p.adjust.survival < .05, "a", "b")), size = 1) +
#    facet_wrap(~replicate) +
#    theme(legend.position = "none") +
#    scale_color_manual(values = c("tomato", "grey60")) +
#    ggtitle(y)
#},
#all.snps.RSI,
#list("gradual", "sudden"), SIMPLIFY = FALSE)

# Hurst SNP location. Are they more random on the genomes depending on treatment?

# unknown contig length, set to max (11685)

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

h.exp <- list(gradual = fit.hurst(mut.pos$gradual, d = 2000) %>% bind_rows(.id = "time"),
              sudden = fit.hurst(mut.pos$sudden, d = 2000) %>% bind_rows(.id = "time")) %>% bind_rows(.id = "treatment")

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

# pull only mutations that reach a threshold by treatment and time

mut.pos.thresh <- list(gradual = gradual.mtx %>% lapply(apply.filter, freq.thresh = .3),
                       sudden = sudden.mtx %>% lapply(apply.filter, freq.thresh = .3)) %>% lapply(lapply, function(x) {
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

h.exp.thresh <- fit.hurst(list(data = mut.pos.thresh), d = 2400) %>% bind_rows(.id = "treatment")

# make a figure

hurst.gg <- ggarrange(h.exp.gg + labs(tag = "a"),
                      h.exp.pop.gg + theme(axis.ticks.x = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.title.y = element_blank()) + labs(tag = "b"),
                      common.legend = TRUE, widths = c(1, .5))

ggsave("hurst.pdf", hurst.gg, width = 6.85, height = 4)

# play with binary sequences from a Markov chain

bin.mc <- function(pi.12, pi.21) inverse.rle(list(values = rep(0:1, 10000),
                                                  lengths = 1 + rgeom(2*10000, rep(c(pi.12, pi.21), 10000))))

H.test <- lapply(setNames(seq(.01, .99, .01), nm = 1:99), function(x) {
  set.seed(x*10)
  lst <- list(bin.mc(x, 1 - x)[1:11685],
              bin.mc(x, x)[1:11685])
  names(lst) <- c(paste0(1, ":", x, "_", 1 - x), paste0(2, ":", x, "_", x))
  lst %>% bind_rows() %>% t() %>% as.data.frame()
}) %>% bind_rows()

hurst.test <- fit.hurst(list(test = H.test), d = 2000) %>% bind_rows(.id = "rep")

