source("source_me.R")

# load matrices

gradual.mtx <- lapply(setNames(1:9, nm = c(98, 101, 102,
                                           106, 107, 109,
                                           112, 115, 116) %>% as.character()), function(x) {
                                             y <- read_xlsx("data/Gradual.xlsx", sheet = x) %>% as.data.frame()
                                             rownames(y) <- y$'...1'
                                             y[-1]
                                           })

sudden.mtx <- lapply(setNames(1:9, nm = c(2, 5, 6,
                                          10, 11, 13,
                                          16, 19, 20) %>% as.character()), function(x) {
  y <- read_xlsx("data/Sudden.xlsx", sheet = x) %>% as.data.frame()
  rownames(y) <- y$'...1'
  y[-1]
})

# fit SNPs frequencies to Taylor's law

gradual.TL <- fit.TL(gradual.mtx, zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = TRUE)
sudden.TL <- fit.TL(sudden.mtx, zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = TRUE)

# plot parameters

snps.TL.params <- rbind(gradual.TL$params[-1] %>% cbind(data = "gradual"),
                        sudden.TL$params[-1] %>% cbind(data = "sudden"))

ggplot(snps.TL.params, aes(V, beta, color = data)) + geom_point()

# Taylor's plot

plot.power_law(gradual.TL$log_log.mean_sd,
               gradual.TL$log_log.mean_sd %>%
                 lapply(function(x) setNames(rownames(x) %>%
                                               gsub("[A-Z+-]*", "", .) %>% gsub("\\*", "", .) %>% as.numeric(),
                                             nm = rownames(x))),
               metadata.range = c(0, 11703), legend.title = "position", l10 = TRUE) +
  ggtitle("gradual") + xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_density_2d(color = "tomato")

plot.power_law(sudden.TL$log_log.mean_sd,
               sudden.TL$log_log.mean_sd %>%
                 lapply(function(x) setNames(rownames(x) %>%
                                               gsub("[A-Z+-]*", "", .) %>% gsub("\\*", "", .) %>% as.numeric(),
                                             nm = rownames(x))),
               metadata.range = c(0, 11703), legend.title = "position", l10 = TRUE) +
  ggtitle("sudden") + xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_density_2d(color = "tomato", )

# test if distances of points to the centroid are greater in sudden

# calculate euclidean distances to centroid

distances.to.centroid <- lapply(list(gradual = gradual.TL$log_log.mean_sd,
                                     sudden = sudden.TL$log_log.mean_sd), function(x) {
                                       lapply(x, function(y) {
                                         
                                         y <- y[!is.na(y$sd),]
                                         
                                         centroid <- list(mean = mean(y$mean),
                                                          sd = mean(y$sd))
                                         
                                         distances <- data.frame(mean.dist = y$mean - centroid$mean,
                                                                sd.dist = y$sd - centroid$sd)
                                         sqrt((distances$mean.dist + distances$sd.dist)^2)
                                       })
                                     })

# plot distributions

lapply(distances.to.centroid, function(x) {
  lapply(x, function(y) {
    data.frame(dist = y)
  }) %>% bind_rows(.id = "rep")
}) %>% bind_rows(.id = "data") %>% ggplot(aes(dist, color = data, group = paste0(data, rep))) +
  geom_density()

# calculate entropy

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
    -sum(c(lapply(snps, function(y) y * log(y, 2)), ref.freq * log(ref.freq, 2)) %>% unlist())
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

entropy.mtx %>% lapply(function(x) {
  pivot_longer(x, -`0`, names_to = "time", values_to = "H")[-1]
}) %>% bind_rows(.id = "data") %>%
  ggplot(aes(time %>% factor(levels = c("4", "7", "10", "13", "16", "19", "22", "25")), H, color = data)) +
  geom_boxplot() +
  xlab("time") +
  theme(legend.title = element_blank())

entropy.TL$log_log.mean_sd %>%
  bind_rows(.id = "title") %>%
  mutate(mean = log10(exp(mean)), sd = log10(exp(sd))) %>% ggplot(aes(mean, sd, color = title)) +
  ggtitle("entropy") + xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y ~ x) +
  theme(legend.title = element_blank())

# calculate entropy per ORF


