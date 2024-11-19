source("source_me.R")

# get positions all mutations

mut.pos <- list(gradual = gradual.mtx,
                sudden = sudden.mtx) %>% lapply(lapply, function(x) {
                  apply(x, 2, function(y) {
                    y <- y[!grepl("\\+|\\-", names(y))]
                    pos <- names(y)[ifelse(y > 0, TRUE, FALSE)] %>% gsub("\\D+", "", .) %>% as.numeric()
                    genome <- rep(0, 11703)
                    genome[pos] <- 1
                    genome
                  }) %>% t() %>% as.data.frame()
                })

h.exp <- list(gradual = fit.hurst(mut.pos$gradual, d = 2000) %>% bind_rows(.id = "time"),
              sudden = fit.hurst(mut.pos$sudden, d = 2000) %>% bind_rows(.id = "time")) %>% bind_rows(.id = "Treatment") %>%
  dplyr::mutate(Treatment := Treatment %>%
                  gsub("^g", "G", .) %>%
                  gsub("^s", "S", .))

h.exp.gg <- ggplot(h.exp, aes(time %>% factor(levels = time %>% unique()), He, color = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_dodge(width = .75)) +
  xlab("Passage") + ylab( expression(italic("H"))) +
  geom_hline(aes(yintercept = .7), linetype = "dashed")

# test effect of time + sample + data

h.exp.lm <- lm(He ~ Treatment + time + data, h.exp %>% mutate(time := time %>% as.numeric()))

# only below a threshold

mut.pos.filter <- list(gradual = gradual.mtx %>% lapply(apply.filter, freq.thresh = .3, rev = TRUE),
                       sudden = sudden.mtx %>% lapply(apply.filter, freq.thresh = .3, rev = TRUE)) %>% lapply(lapply, function(x) {
                         apply(x, 2, function(y) {
                           y <- y[!grepl("\\+|\\-", names(y))]
                           pos <- names(y)[ifelse(y > 0, TRUE, FALSE)] %>% gsub("\\D+", "", .) %>% as.numeric()
                           genome <- rep(0, 11703)
                           genome[pos] <- 1
                           genome
                         }) %>% t() %>% as.data.frame()
                       })

h.exp.filter <- list(gradual = fit.hurst(mut.pos.filter$gradual, d = 2000) %>% bind_rows(.id = "time"),
                     sudden = fit.hurst(mut.pos.filter$sudden, d = 2000) %>% bind_rows(.id = "time")) %>% bind_rows(.id = "treatment")

# NA where He could not be calculated

h.exp.filter[h.exp.filter == 0] <- NA

h.exp.filter.gg <- ggplot(h.exp.filter, aes(time %>% factor(levels = time %>% unique()), Hrs, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_dodge(width = .75)) +
  xlab("passage") + ylab("Hrs") +
  geom_hline(aes(yintercept = .7), linetype = "dashed")

h.exp.filter.lm <- lm(Hrs ~ treatment * time + data, h.exp.filter %>% mutate(time := time %>% as.numeric()))

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
  xlab("passage") + ylab("He") +
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

h.exp.pop <- fit.hurst(mut.pos.pop, d = 550) %>% bind_rows(.id = "population") %>% filter(!is.na(He))

h.exp.pop.gg <- ggplot(h.exp.pop, aes(data, He, color = data, group = data)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  xlab("") + ylab("He") +
  geom_hline(aes(yintercept = .7), linetype = "dashed")

h.exp.pop.wilcox <- wilcox.test(He ~ data,
                                h.exp.pop)

h.exp.pop.lm <- lm(He ~ data, h.exp.pop)

# pull all mutations by treatment and time

mut.cat.pos <- list(gradual = gradual.mtx,
                    sudden = sudden.mtx) %>% lapply(lapply, function(x) {
                      apply(x, 2, function(y) {
                        y <- y[!grepl("\\+|\\-", names(y))]
                        pos <- names(y)[ifelse(y > 0, TRUE, FALSE)] %>% gsub("\\D+", "", .) %>% as.numeric()
                        genome <- rep(0, 11703)
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
                           genome <- rep(0, 11703)
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

# aggregation from alleles distance



# make a figure

ggsave("6_hurst1.pdf", h.exp.gg, width = 6.85, height = 4)
ggsave("6_hurst1.png", h.exp.gg, width = 6.85, height = 4, dpi = 600)

hurst.gg <- ggarrange(h.exp.gg + labs(tag = "a"),
                      h.exp.pop.gg + theme(axis.ticks.x = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.title.y = element_blank()) + labs(tag = "b"),
                      common.legend = TRUE, widths = c(1, .5))

ggsave("hurst2.pdf", hurst.gg, width = 6.85, height = 4)
