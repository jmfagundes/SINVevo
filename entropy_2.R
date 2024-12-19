source("source_me.R")

entropy.mtx <- lapply(list(Gradual = gradual.mtx,
                           Sudden = sudden.mtx), function(x) {
                             lapply(x, function(y) {
                               df <- apply(y, 2, calc.entropy, return.sum = FALSE) %>%
                                 bind_rows() %>% t() %>% as.data.frame()
                               colnames(df) <- c("4", "7", "10", "13", "16", "19", "22", "25")
                               df %>% t() %>% as.data.frame()
                             })
                           })

entropy.TL <- entropy.mtx %>% lapply(fit.TL, normalize = FALSE)

entropy.TL.gg <- entropy.TL %>% lapply(function(x) x$log_log.mean_sd %>% bind_rows(.id = "Population")) %>%
  bind_rows(.id = "Treatment") %>%
  mutate(mean = log10(exp(mean)), sd = log10(exp(sd))) %>%
  ggplot(aes(mean, sd, color = Treatment)) +
  xlab("Mean (log)") + ylab("Standard deviation (log)") +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y ~ x,
              show.legend = FALSE,
              se = FALSE) +
  geom_smooth(aes(group = Population),
              method = "lm",
              formula = y ~ x,
              show.legend = FALSE,
              se = FALSE,
              linewidth = .2) +
  theme(legend.title = element_blank())

# test beta

entropy.TL.params <- entropy.TL %>% lapply(function(x) x$params) %>% bind_rows(.id = "Treatment")

entropy.TL.params %>% ggplot(aes(Treatment, beta)) + geom_boxplot()
entropy.TL.params %>% ggplot(aes(Treatment, V)) + geom_boxplot()

wilcox.test(entropy.TL$Gradual$params$beta, entropy.TL$Sudden$params$beta)
t.test(entropy.TL$Gradual$params$beta, entropy.TL$Sudden$params$beta)

lm(beta ~ Treatment + V, entropy.TL.params)
lm(V ~ Treatment + beta, entropy.TL.params)

manova(cbind(V, beta) ~ Treatment, entropy.TL.params)
