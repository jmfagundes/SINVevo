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

# calculate entropy per ORF


# fit SNPs frequencies to Taylor's law

gradual.TL <- fit.TL(gradual.mtx, zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = FALSE)
sudden.TL <- fit.TL(sudden.mtx, zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = FALSE)

# plot parameters

snps.TL.params <- rbind(gradual.TL$params[-1] %>% cbind(data = "gradual"),
                        sudden.TL$params[-1] %>% cbind(data = "sudden"))

all.snps.TL.param.gg <- ggplot(snps.TL.params, aes(V, beta, color = data)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)

# Taylor's plot

all.snps.TP.gradual.gg <- plot.power_law(gradual.TL$log_log.mean_sd,
               gradual.mtx %>% calc.0_rate(),
               metadata.range = c(0, 1), legend.title = "0 rate", l10 = TRUE) +
  ggtitle("gradual") + xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red")

all.snps.TP.sudden.gg <- plot.power_law(sudden.TL$log_log.mean_sd,
               sudden.mtx %>% calc.0_rate(),
               metadata.range = c(0, 1), legend.title = "0 rate", l10 = TRUE) +
  ggtitle("sudden") + xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red")

# reconstruct haplotypes from euclidean distance matrices calculated from a Pearson's coefficients matrix
# only include SNPs present in at least n.time points

n.time <- 2

snps.cor <- lapply(list(gradual = gradual.mtx,
                        sudden = sudden.mtx),
                   lapply, function(x) {
                     x <- x[apply(x, 1, function(y) length(y[y > 0]) >= n.time),]
                     sqrt(2 * (1 - cor(t(x))))
                   })

# estimate number of clusters

n.clusters <- lapply(snps.cor, lapply, function(x) fviz_nbclust(x, pam, k.max = nrow(x) - 1)$plot_env$data$y %>% which.max())

# generate clusters

haplotypes <- mapply(function(x, y) {
  mapply(function(i, j) {
    pam(i, j)
  }, x, y, SIMPLIFY = FALSE)
}, snps.cor, n.clusters,
SIMPLIFY = FALSE)

# plot SNPs colored by haplotype

haplo.SNP.color.gg <- mapply(function(x, y) {
  mapply(function(i, j) {
    
    i <- i %>% filter(rownames(i) %in% names(j$clustering))
    
    i %>% mutate(SNP = rownames(i), haplotype = paste0("haplotype ", j$clustering)) %>%
      pivot_longer(-c(SNP, haplotype), names_to = "timepoint", values_to = "frequency")
    
  }, x, y, SIMPLIFY = FALSE) %>% bind_rows(.id = "rep")
}, list(gradual = gradual.mtx,
        sudden = sudden.mtx),
haplotypes, SIMPLIFY = FALSE) %>% 
  lapply(function(x) {
    ggplot(x, aes(timepoint, frequency, color = haplotype)) + geom_point() +
      facet_wrap(~rep) +
      theme(legend.position = "none")
    })

# calculate frequency of each haplotype as the mean of its SNPs frequencies at each timepoint

haplo.freqs <- mapply(function(x, y) {
  mapply(function(i, j) {
    
    i <- i[apply(i, 1, function(z) length(z[z > 0]) >= n.time),] # only SNPs present at >= n.time timepoints
    
    freqs <- apply(i, 2, function(z) {
      aggregate(z, by = list(j$clustering), mean)$x
    }) %>% as.data.frame()
    
    rownames(freqs) <- paste0("haplotype ", 1:nrow(freqs))
    
    freqs
    
  }, x, y, SIMPLIFY = FALSE)
}, list(gradual = gradual.mtx,
        sudden = sudden.mtx),
haplotypes,
SIMPLIFY = FALSE)

# Taylor's law with haplotypes
# don't normalize because of background mutations

haplo.freqs.TL <- lapply(haplo.freqs, function(x) {
  fit.TL(x, zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = FALSE, min.rows = 10)
})

# plot

haplo.TL.params <- rbind(haplo.freqs.TL$gradual$params[-1] %>% cbind(data = "gradual"),
                         haplo.freqs.TL$sudden$params[-1] %>% cbind(data = "sudden"))

ggplot(haplo.TL.params, aes(V, beta, color = data)) + geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)

mapply(function(x, y, z) {
  
  plot.power_law(x$log_log.mean_sd,
                 z %>% calc.0_rate(),
                 metadata.range = c(0, 1), legend.title = "0 rate",
                 l10 = TRUE) +
    ggtitle(y) + xlab("mean (log)") + ylab("standard deviation (log)") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red")
}, haplo.freqs.TL, c("gradual", "sudden"), haplo.freqs,
SIMPLIFY = FALSE)

# entropy based on haplotype frequencies
# not accurate because background mutations are counted as different haplotypes

haplo.entropy.mtx <- lapply(haplo.freqs, function(x) lapply(x, function(y) {
  apply(y, 2, function(z) {
    
    # remove haplotypes with frequency = 0 and normalize
    
    z <- z[z > 0]
    z <- z / sum(z)
    
    shannon.entropy(z)
    
  })
}) %>% bind_rows() %>% as.data.frame() )

haplo.entropy.TL <- fit.TL(haplo.entropy.mtx, normalize = FALSE)

haplo.entropy.mtx %>% lapply(function(x) {
  pivot_longer(x, -`0`, names_to = "time", values_to = "H")[-1]
}) %>% bind_rows(.id = "data") %>%
  ggplot(aes(time %>% factor(levels = c("4", "7", "10", "13", "16", "19", "22", "25")), H, color = data)) +
  geom_boxplot() +
  xlab("time") +
  theme(legend.title = element_blank())

haplo.entropy.TL$log_log.mean_sd %>%
  bind_rows(.id = "title") %>%
  mutate(mean = log10(exp(mean)), sd = log10(exp(sd))) %>%
  ggplot(aes(mean, sd, color = title)) +
  ggtitle("entropy") + xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y ~ x) +
  theme(legend.title = element_blank())
