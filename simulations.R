source("source_me.R")

# test H = I*d**gamma in Gompertz and random walk

gompz <- function(t, alpha, beta, gamma) alpha * (exp(1) ** (-beta * exp(1) ** (-gamma * t)))

gompz.H.tb <- data.frame(t = 1:100,
                         freq = gompz(1:100, 1, 10, .1),
                         d = gompz(1:101, 1, 10, .1) %>% diff() %>% abs())#,
                         #dprevious = gompz(0:100, 1, 10, .1) %>% diff() %>% abs())

gompz.H.tb$H <- lapply(gompz.H.tb$freq, function(x) shannon.entropy(c(x, 1 - x))) %>% unlist()

#gompz.H.tb$d <- gompz.H.tb$d/max(gompz.H.tb$d)
#gompz.H.tb$dprevious <- gompz.H.tb$dprevious/max(gompz.H.tb$dprevious)

gompz.H.tb$quasiH <- lapply(1:nrow(gompz.H.tb), function(x) {
  freq <- gompz.H.tb$freq[x]
  -freq*log(freq, 2)
}) %>% unlist()

#gompz.H.tb$quasiH <- gompz.H.tb$quasiH / max(gompz.H.tb$quasiH)

gompz.gg <- gompz.H.tb %>% pivot_longer(-t) %>% ggplot(aes(t, value, color = name)) + geom_line() +
  geom_vline(aes(xintercept = gompz.H.tb[gompz.H.tb$H %>% which.max(), "t"]))

ggsave("gompz.pdf", gompz.gg, width = 6.85, height = 5)

#gompz.H.tb %>% 
#  dplyr::mutate(d := log(d), H := log(H), mean.d := log(mean.d)) %>% dplyr::select(-freq) %>%
#  pivot_longer(-t) %>%
#  ggplot(aes(t, value, color = name)) + geom_line() +
#  geom_vline(aes(xintercept = gompz.H.tb[gompz.H.tb$H %>% which.max(), "t"]))

#gompz.H.tb %>% 
#  dplyr::mutate(d := d, H := H) %>% dplyr::select(-freq) %>%
#  pivot_longer(-t) %>%
#  ggplot(aes(t, value, color = name)) + geom_line() +
#  geom_vline(aes(xintercept = gompz.H.tb[gompz.H.tb$H %>% which.max(), "t"]))

# test gompertz vs random

#gompz.tb <- data.frame(gompz.beta.5 = gompz(0:8, 1, 5, 1),
#                       gompz.beta.20 = gompz(0:8, 1, 20, 1),
#                       gompz.beta.50.gamma.2 = gompz(0:8, 1, 50, 2),
#                       gompz.beta.10.gamma.2 = gompz(0:8, 1, 10, 1),
#                       gompz.beta.10.gamma..6 = gompz(0:8, 1, 10, 1)) %>% t()

# lots of early fixations

#gompz.tb <- lapply(seq(.5, 5, .05), function(y) {
#  lapply(seq(1, 10, .5) %>% setNames(nm = .), function(x) gompz(1:8, 1, x * 10, y)) %>% bind_rows()
#}) %>% bind_cols()

gompz.tb <- lapply(seq(1.1, 1.2, .05), function(y) {
  lapply(seq(1, 10, .5) %>% setNames(nm = .), function(x) gompz(1:8, 1, x * 10, y)) %>% bind_rows()
}) %>% bind_cols()

gompz.tb[gompz.tb < .01] <- 0
gompz.tb[gompz.tb > .99] <- 1

gompz.tb <- as.data.frame(gompz.tb) %>% t() %>% as.data.frame()

gompz.tb.2 <- lapply(c(1:10, 2:10 * 10) %>% setNames(nm = .), function(x) gompz(1:100, 1, x * 10, .1)) %>% bind_rows()

gompz.tb.2[gompz.tb.2 < .01] <- 0
gompz.tb.2[gompz.tb.2 > .99] <- 1

gompz.tb.2 <- as.data.frame(gompz.tb.2) %>% t() %>% as.data.frame()

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

# plot walks

walks.gg <- list(gompertz.8 = gompz.tb,
     random.2 = random.walk.tb,
     random.05 = random.walk.tb.2,
     random.01 = random.walk.tb.3) %>% lapply(function(x) {
       
       t(x) %>% as.data.frame() %>% dplyr::mutate(t = 1:8) %>%
         pivot_longer(-t)
       
     }) %>% bind_rows(.id = "data") %>% ggplot(aes(t, value, group = name)) +
  geom_line(linewidth = .05) + facet_wrap(~data)

fit.entropy <- fit.time.sd(list(gompertz.8 = gompz.tb,
                                #gompertz.100 = gompz.tb.2,
                                random.2 = random.walk.tb,
                                random.05 = random.walk.tb.2,
                                random.01 = random.walk.tb.3),
                           return.mean = FALSE, not.fix = TRUE, use.entropy = TRUE)

fit.test.gg <- fit.entropy$log_log.mean_sd %>% bind_rows(.id = "treatment") %>%
  ggplot(aes(mean, sd)) + geom_point(shape = 20, stroke = 0, size = .75) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + 
  #xlab("entropy (log)") + ylab("absolute difference from next value (log)") +
  xlab("H (log)") + ylab("d (log)") +
  facet_wrap(~treatment)

ggsave("test_Hd.pdf", ggarrange(walks.gg, fit.test.gg, nrow = 2), width = 6.85, height = 8)

#fit.entropy$log_log.mean_sd %>% bind_rows(.id = "treatment") %>%
#  ggplot(aes(exp(mean), exp(sd))) + geom_point(shape = 20, stroke = 0, size = .75) + 
#  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + 
#  xlab("entropy (log)") + ylab("absolute difference from next value (log)") + facet_wrap(~treatment)

#fit.entropy$log_log.mean_sd$gompertz.100 %>%
#  ggplot(aes(mean, sd)) + geom_point(shape = 20, stroke = 0, size = .75) + 
#  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + 
#  xlab("entropy (log)") + ylab("absolute difference from next value (log)")

#fit.entropy$log_log.mean_sd$gompertz.100 %>%
#  ggplot(aes(exp(mean), exp(sd))) + geom_point(shape = 20, stroke = 0, size = .75) + 
#  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2, color = "red") + 
#  xlab("entropy (log)") + ylab("absolute difference from next value (log)")

#fit.entropy$log_log.mean_sd %>% bind_rows(.id = "treatment") %>%
#  ggplot(aes(mean, sd, color =treatment)) +
#  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .2) + 
#  xlab("entropy (log)") + ylab("absolute difference from next value (log)")

# test TL

gompz.TL.classic <- fit.TL(list(gompertz.8 = gompz.tb,
                                gompertz.100 = gompz.tb.2,
                                random.2 = random.walk.tb,
                                random.05 = random.walk.tb.2,
                                random.01 = random.walk.tb.3,
                                mix.2 = bind_rows(gompz.tb, random.walk.tb),
                                mix.05 = bind_rows(gompz.tb, random.walk.tb.2),
                                mix.01 = bind_rows(gompz.tb, random.walk.tb.3)),
                           zero.rate.threshold = NULL, normalize = FALSE,
                           remove.zeros = FALSE)

gompz.random.TP.gg <- plot.power_law(gompz.TL.classic$log_log.mean_sd,
               list(gompertz.8 = gompz.tb,
                    gompertz.100 = gompz.tb.2,
                    random.2 = random.walk.tb,
                    random.05 = random.walk.tb.2,
                    random.01 = random.walk.tb.3,
                    mix.2 = bind_rows(gompz.tb, random.walk.tb),
                    mix.05 = bind_rows(gompz.tb, random.walk.tb.2),
                    mix.01 = bind_rows(gompz.tb, random.walk.tb.3)) %>% calc.0_rate(),
               metadata.range = NULL, legend.title = "0 rate", l10 = TRUE) +
  xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_smooth(method = "lm", se = FALSE, linewidth = .2)

# per 0 random walks

pipe.per_0.random <- function(n = 1:8,
                              remove.zeros = FALSE,
                              min.rows = 5,
                              sd.from.max = FALSE) {
  
  lapply(setNames(n, nm = as.character(n)), function(x) {
    
    # fit SNPs frequencies to Taylor's law
    
    fit.TL(list(mix.1 = bind_rows(gompz.tb, random.walk.tb),
                mix.5 = bind_rows(gompz.tb, random.walk.tb.2),
                random.2 = random.walk.tb,
                random.05 = random.walk.tb.2,
                random.01 = random.walk.tb.3) %>% lapply(apply.filter, x, exact = TRUE),
           zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = remove.zeros, min.rows = min.rows, sd.from.max = sd.from.max)
  })
}

random.TL.per_0 <- pipe.per_0.random()

gompz.random.pre_0.params <- lapply(random.TL.per_0, function(x) {
  x$params
}) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se = FALSE, linewidth = .2)
