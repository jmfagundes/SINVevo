source("source_me.R")

# calculate entropy

shannon.entropy <- function(freq.lst) -sum(lapply(freq.lst, function(x) x * log(x, 2)) %>% unlist())

calc.entropy <- function(alt.freq.lst) {
  
  if (sum(alt.freq.lst) == 0) return(0)
  
  # remove indels
  
  alt.freq.lst <- alt.freq.lst[!grepl("\\+|\\-", names(alt.freq.lst))]
  
  # remove 0 or 1 variants
  
  alt.freq.lst <- alt.freq.lst[!alt.freq.lst %in% c(0, 1)]
  
  # get snps at same site
  
  sites <- names(alt.freq.lst) %>%
    gsub("[A-Z+-]*", "", .) %>% gsub("\\*", "", .) %>% unique()
  
  # calcute entropy
  
  entropies <- lapply(setNames(nm = sites), function(x) {
    
    snps <- alt.freq.lst[grepl(paste0("^", x, "[ACTG]|^", x, "\\*"), names(alt.freq.lst))]
    ref.freq <- 1 - sum(snps)
    shannon.entropy(c(snps, ref.freq))
  })
  entropy.sum <- sum(entropies %>% unlist())
  
  return(entropy.sum)
}

entropy.mtx <- lapply(list(gradual = gradual.mtx,
                           sudden = sudden.mtx), function(x) {
                             lapply(x, function(y) {
                               apply(y, 2, calc.entropy)
                             }) %>% bind_rows() %>% as.data.frame()
                           })

entropy.TL <- fit.TL(entropy.mtx, normalize = FALSE)

entropy.boxplot.gg <- entropy.mtx %>% lapply(function(x) {
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
  pivot_longer(x, -`0`, names_to = "time", values_to = "H")[-1]
}) %>% bind_rows(.id = "data") %>%
  ggplot(aes(time %>% factor(levels = c("4", "7", "10", "13", "16", "19", "22", "25")), H, color = data)) +
  geom_boxplot() +
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

