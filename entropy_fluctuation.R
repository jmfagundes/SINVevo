source("source_me.R")

entropy.mtx <- lapply(list(gradual = gradual.mtx,
                           sudden = sudden.mtx), function(x) {
                             lapply(x, function(y) {
                               apply(y, 2, calc.entropy)
                             }) %>% bind_rows() %>% as.data.frame()
                           })

entropy.TL <- fit.TL(entropy.mtx, normalize = FALSE)

entropy.boxplot.gg <- entropy.mtx %>% lapply(function(x) {
  x$`0` <- 0
  pivot_longer(x, -`0`, names_to = "time", values_to = "H")[-1]
}) %>% bind_rows(.id = "data") %>%
  ggplot(aes(time %>% factor(levels = c("4", "7", "10", "13", "16", "19", "22", "25")), H, color = data)) +
  geom_boxplot() +
  xlab("time") +
  theme(legend.title = element_blank()) +
  labs(tag = "a")

entropy.TL.gg <- entropy.TL$log_log.mean_sd %>%
  bind_rows(.id = "title") %>%
  mutate(mean = log10(exp(mean)), sd = log10(exp(sd))) %>%
  ggplot(aes(mean, sd, color = title)) +
  xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y ~ x,
              show.legend = FALSE) +
  theme(legend.title = element_blank())  +
  labs(tag = "b")

entropy.full.genome.gg <- ggarrange(entropy.boxplot.gg, entropy.TL.gg, ncol = 1,
                                    legend.grob = get_legend(entropy.TL.gg), legend = "right")

# only SNPs above filter parameters

entropy.filter.mtx <- lapply(list(gradual = gradual.mtx %>% lapply(apply.filter, n.time, freq.thresh),
                                  sudden = sudden.mtx  %>% lapply(apply.filter, n.time, freq.thresh)), function(x) {
                                    lapply(x, function(y) {
                                      apply(y, 2, calc.entropy)
                                    }) %>% bind_rows() %>% as.data.frame()
                                  })

entropy.filter.TL <- fit.TL(entropy.filter.mtx, normalize = FALSE)

entropy.filter.boxplot.gg <- entropy.filter.mtx %>% lapply(function(x) {
  x$`0` <- 0
  pivot_longer(x, -`0`, names_to = "time", values_to = "H")[-1]
}) %>% bind_rows(.id = "data") %>%
  ggplot(aes(time %>% factor(levels = c("4", "7", "10", "13", "16", "19", "22", "25")), H, color = data)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_dodge(width = .75)) +
  xlab("time") +
  theme(legend.title = element_blank()) +
  labs(tag = "a")

entropy.filter.TL.gg <- entropy.filter.TL$log_log.mean_sd %>%
  bind_rows(.id = "title") %>%
  mutate(mean = log10(exp(mean)), sd = log10(exp(sd))) %>%
  ggplot(aes(mean, sd, color = title)) +
  xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y ~ x,
              show.legend = FALSE) +
  theme(legend.title = element_blank())  +
  labs(tag = "b")

entropy.filter.full.genome.gg <- ggarrange(entropy.filter.boxplot.gg, entropy.filter.TL.gg, ncol = 1,
                                           legend.grob = get_legend(entropy.TL.gg), legend = "right")

ggsave("entropy_TL.pdf",entropy.filter.full.genome.gg, width = 6.85, height = 6)

# calculate entropy per ORF

