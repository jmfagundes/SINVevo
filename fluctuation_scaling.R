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
TL.classic$snps.gg <- lapply(TL.classic$snps.gg, function(x) {x + xlim(-3, 0) + ylim(-2.5, -.45) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    geom_smooth(aes(mean, sd), method = "lm", linetype = "dashed", se = FALSE, color = "black", linewidth = .2) +
    scale_color_discrete(direction = -1) 
})

# Taylor's law with SNPs above threshold

TL.above.thresh <- pipe(gradual.mtx %>% lapply(apply.filter, n.time, .1),
                        sudden.mtx %>% lapply(apply.filter, n.time, .1),
                        FALSE, FALSE)

# Taylor's law with SNPs below threshold

TL.below.thresh <- pipe(gradual.mtx %>% lapply(apply.filter, n.time, .1, rev = TRUE),
                        sudden.mtx %>% lapply(apply.filter, n.time, .1, rev = TRUE),
                        FALSE, FALSE)

# how beta scales with V

all.snps.V.beta.scale <- lapply(list(gradual = TL.classic$TL$gradual$params,
                                     sudden = TL.classic$TL$sudden$params),
                                function(x) {
                                  lm(beta ~ V, x) %>% summary()
                                })

snps.filter.V.beta.scale <- lapply(list(gradual = TL.above.thresh$TL$gradual$params,
                                        sudden = TL.above.thresh$TL$sudden$params),
                                   function(x) {
                                     lm(beta ~ V, x) %>% summary()
                                   })

snps.rev.V.beta.scale <- lapply(list(gradual = TL.below.thresh$TL$gradual$params,
                                     sudden = TL.below.thresh$TL$sudden$params),
                                function(x) {
                                  lm(beta ~ V, x) %>% summary()
                                })

# Taylor's law per 0 rate

pipe.per_0 <- function(sd.from.max = FALSE) {
  
  lapply(setNames(1:8, nm = as.character(1:8)), function(x) {
    
    # fit SNPs frequencies to Taylor's law
    
    gradual.TL <- fit.TL(gradual.mtx %>% lapply(apply.filter, x, exact = TRUE),
                         zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = FALSE, min.rows = 3, sd.from.max = sd.from.max)
    sudden.TL <- fit.TL(sudden.mtx %>% lapply(apply.filter, x, exact = TRUE),
                        zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = FALSE, min.rows = 3, sd.from.max = sd.from.max)
    
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

facet.V.beta.gg <- lapply(TL.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) + facet_wrap(~treatment)

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
                               }) %>% bind_rows(.id = "treatment") %>% dplyr::select(!`0`) %>%
  pivot_longer(-c(treatment, rep, n), names_to = "time", values_to = "counts") %>%
  ggplot(aes(time %>% factor(levels = sort(time %>% unique() %>% as.numeric())), counts, color = n)) +
  geom_boxplot(linewidth = .2, outlier.size = .5) + facet_wrap(~treatment) + xlab("time")

# Taylor's law removing zeros

TL.no_zero <- pipe(gradual.mtx, sudden.mtx, TRUE, FALSE)

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

ggsave("snps_TL.pdf", snps.taylor.gg, width = 6.85, height = 3.5)

combined.gg <- ggarrange(facet.V.beta.gg + labs(tag = "a") + theme(legend.position = "none"),
                         location.SNPs.gg + labs(tag = "b") + theme(legend.position = "none"),
                         nrow = 2,
                         common.legend = TRUE,
                         legend.grob = get_legend(facet.V.beta.gg + geom_point(shape = 15, size = 6) +
                                                    theme(legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1))),
                         legend = "bottom")

ggsave("params_TL.pdf", combined.gg, width = 6.85, height = 6)

# Taylor's law from the max

TL.max <- pipe(gradual.mtx, sudden.mtx, FALSE, TRUE)
TL.max$snps.gg <- lapply(TL.max$snps.gg, function(x) {
  x + xlab("max (log)") + ylab("standard deviation of the max (log)") +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
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
                                legend = "top", widths = c(1, .94)) %>% annotate_figure(left = "standard deviation (log)", bottom = "mean (log)")

TL.max.per.zero <- pipe.per_0(sd.from.max = TRUE)

max.facet.V.beta.gg <- lapply(TL.max.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2) + facet_wrap(~treatment)

max.joint.gg <- lapply(TL.max.per.zero, function(x) {
  bind_rows(lapply(x$TL, function(y) y$params), .id = "treatment")
}) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n, linetype = treatment)) +
  geom_smooth(method = "lm", se = FALSE)

# with no zero

TL.max.no_zero <- pipe(gradual.mtx, sudden.mtx, TRUE, TRUE)
TL.max.no_zero$snps.gg <- lapply(TL.max.no_zero$snps.gg, function(x) x + xlab("max (log)") + ylab("standard deviation of the max (log)"))

# sd from the max ~ stability

gradual.JL <- fit.JL(gradual.mtx,
                     zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = FALSE)
sudden.JL <- fit.JL(sudden.mtx,
                    zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = FALSE)

# plot parameters

snps.JL.params <- rbind(gradual.JL$params[-1] %>% cbind(data = "gradual"),
                        sudden.JL$params[-1] %>% cbind(data = "sudden"))

all.snps.JL.param.gg <- ggplot(snps.JL.params, aes(V, beta, color = data)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)

# Taylor's plot

all.snps.JP.gg <- mapply(function(x, y) {
  ggplot(x %>% bind_rows(.id = "rep"),
         aes(max.time %>% factor(levels = max.time %>% unique()),
             sd.from.max %>% log10())) + geom_boxplot() +
    xlab("time point") + ylab("standard deviation from the max (log)") +
    facet_wrap(~rep) +
    ggtitle(y)
},
list(gradual = gradual.JL$max_sd,
     sudden = sudden.JL$max_sd),
list("gradual", "sudden"), SIMPLIFY = FALSE)

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

