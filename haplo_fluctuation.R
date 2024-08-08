source("source_me.R")

# reconstruct haplotypes from euclidean distance matrices calculated from a Pearson's coefficients matrix

snps.cor <- lapply(list(gradual = gradual.mtx,
                        sudden = sudden.mtx),
                   lapply, function(x) {
                     sqrt(2 * (1 - cor(t(x %>% apply.filter(n.time, freq.thresh)))))
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
    
    i <- i %>% apply.filter(n.time, freq.thresh)
    
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
  fit.TL(x, zero.rate.threshold = NULL, normalize = FALSE, remove.zeros = FALSE, min.rows = 5)
})

# plot

haplo.TL.params <- rbind(haplo.freqs.TL$gradual$params[-1] %>% cbind(data = "gradual"),
                         haplo.freqs.TL$sudden$params[-1] %>% cbind(data = "sudden"))

haplo.TL.param.gg <- ggplot(haplo.TL.params, aes(V, beta, color = data)) + geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)

haplo.TP.gg <- mapply(function(x, y, z) {
  
  plot.power_law(x$log_log.mean_sd,
                 z %>% calc.0_rate(),
                 metadata.range = c(0, 1), legend.title = "0 rate",
                 l10 = TRUE) +
    ggtitle(y) + xlab("mean (log)") + ylab("standard deviation (log)") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red")
}, haplo.freqs.TL, c("gradual", "sudden"), haplo.freqs,
SIMPLIFY = FALSE)

# reconstruct haplotypes based on PCA

snps.pca <- lapply(list(gradual = gradual.mtx,
                        sudden = sudden.mtx),
                   lapply, function(x) {
                     
                     prcomp(x %>% apply.filter(n.time, freq.thresh), scale = FALSE)
                   })

snps.pca.gg <- snps.pca %>% lapply(lapply, function(x) {
  fviz_pca_ind(x,
               repel = TRUE)
})

# RSI haplotypes

haplo.rank <- haplo.freqs %>% lapply(function(x) process.rank(matrices = x, method = "random"))

haplo.RSI <- haplo.rank %>% lapply(calc.RSI, shuffle.rep = 1000)

# read CALDER output (clone proportion matrix)

gradual.calder <- lapply(setNames(nm = c("98", "101", "102",
                                         "106", "107", "109",
                                         "112", "115", "116") %>% as.character()), function(x) {
                                           y <- read.table(paste0("SE_results/calder/Gradual.", x, "/Gradual_soln1.csv"),
                                                           sep = ",", row.names = 1, skip = 15, nrows = 8, fill = TRUE) %>% t()
                                           y[is.na(y)] <- 0
                                           y %>% as.data.frame()
                                         })

colnames(gradual.calder$`112`) <- c("4", "7", "10", "13", "16", "19", "22", "25")

sudden.calder <- lapply(setNames(nm = c("2", "5", "6",
                                        "10", "11", "13",
                                        "16", "19", "20") %>% as.character()), function(x) {
                                          y <- read.table(paste0("SE_results/calder/Sudden.", x, "/Sudden_soln1.csv"),
                                                          sep = ",", row.names = 1, skip = 15, nrows = 8, fill = TRUE) %>% t()
                                          y[is.na(y)] <- 0
                                          y %>% as.data.frame()
                                        })

calder.TL <- list(gradual = gradual.calder,
                  sudden = sudden.calder) %>% lapply(function(x) fit.TL(lapply(x, function(y) {
                    y[rowSums(y) > 0,]
                  }), normalize = FALSE, min.rows = 3))

calder.TL %>% lapply(function(x) x$params) %>% bind_rows(.id = "treatment") %>%
  ggplot(aes(V, beta, color = treatment)) + geom_point() + geom_smooth(method = "lm")

calder.TL %>% lapply(function(x) x$log_log.mean_sd) %>% lapply(bind_rows, .id = "population") %>% bind_rows(.id = "treatment") %>%
  ggplot(aes(mean, sd, color = population)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + facet_wrap(~treatment)

# pull variants together

calder.cat <- list(gradual = bind_rows(gradual.calder),
                   sudden = bind_rows(sudden.calder))

calder.cat.TL <- fit.TL(calder.cat %>% lapply(function(x) x[rowSums(x) > 0,]),
                        normalize = FALSE, min.rows = 3)

calder.cat.TL$log_log.mean_sd %>% bind_rows(.id = "treatment") %>%
  ggplot(aes(mean, sd, color = treatment)) + geom_point() + geom_smooth(method = "lm", se = FALSE)

# per 0 rate

calder.cat.TL.per_0 <- lapply(setNames(1:8, nm = as.character(1:8)), function(x) {
  
  # fit SNPs frequencies to Taylor's law
  
  TL <- fit.TL(calder.cat %>% lapply(apply.filter, x, exact = TRUE),
               zero.rate.threshold = NULL, normalize = FALSE, min.rows = 5)
  
  list(TL = TL)
})

calder.cat.TL.per_0.gg <- plot.power_law(calder.cat.TL$log_log.mean_sd,
               calder.cat %>% calc.0_rate(),
               metadata.range = NULL, legend.title = "0 rate", l10 = TRUE) +
  xlab("mean (log)") + ylab("standard deviation (log)") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE)

calder.cat.TL.per_0.params.gg <- calder.cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>% bind_rows(.id = "n") %>%
  ggplot(aes(V, beta, color = n, shape = data)) + geom_point()

#calder.cat.wilcox.V <- wilcox.test(calder.cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
#                                     bind_rows(.id = "n") %>% filter(data == "gradual") %>% dplyr::select(V) %>% unlist(),
#                                   calder.cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
#                                     bind_rows(.id = "n") %>% filter(data == "sudden") %>% dplyr::select(V) %>% unlist(),
#                                   paired = TRUE)
#
#calder.cat.wilcox.beta <- wilcox.test(calder.cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
#                                        bind_rows(.id = "n") %>% filter(data == "gradual") %>% dplyr::select(beta) %>% unlist(),
#                                      calder.cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
#                                        bind_rows(.id = "n") %>% filter(data == "sudden") %>% dplyr::select(beta) %>% unlist(),
#                                      paired = TRUE)

calder.cat.TL.per_0.tb <- calder.cat.TL.per_0 %>% lapply(function(x) x$TL$params) %>%
  bind_rows(.id = "n") %>% mutate(n := as.numeric(n))

calder.cat.V.lm <- lm(V ~ data * n * beta,
                      calder.cat.TL.per_0.tb)

calder.cat.beta.lm <- lm(beta ~ data * n,
                         calder.cat.TL.per_0.tb)

# manova

calder.cat.manova <- manova(cbind(V, beta) ~ data * sparsity, calder.cat.TL.per_0.tb)


