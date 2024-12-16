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

# aggregation from alleles distance

mut.dist.diff <- list(gradual = gradual.mtx,
                      sudden = sudden.mtx) %>% lapply(function(x) {
                        
                        lapply(x, function(y) {
                          
                          alleles <- rownames(y)
                          alleles <- alleles[!grepl("\\+|\\-", alleles)]
                          
                          apply(y, 2, function(z) {
                            
                            z <- z[!grepl("\\+|\\-", names(z))]
                            dists <- alleles[z[!grepl("\\+|\\-", names(z))] > 0] %>% gsub("\\D+", "", .) %>% as.numeric() %>% diff()
                            list(Mean = mean(dists),
                                 `Standard deviation` = sd(dists))
                          }) %>% bind_rows(.id = "Passage")
                        }) %>% bind_rows(.id = "Population")
                      }) %>% bind_rows(.id = "Treatment") %>% dplyr::mutate(`Mean (log)` = log10(Mean),
                                                                            `Standard deviation (log)` = log10(`Standard deviation`))

mut.dist.diff.TL <- mut.dist.diff$Population %>% unique() %>%
  lapply(function(x) {
    
    df <- mut.dist.diff %>% dplyr::filter(Population == x)
    lm <- lm.log_log.mean_sd <- lm(`Standard deviation (log)` ~ `Mean (log)`, data = df)
    fit.summary <- summary(lm.log_log.mean_sd)
    
    list(Treatment = df$Treatment %>% unique() %>%
           gsub("^g", "G", .) %>%
           gsub("^s", "S", .),
         V = fit.summary$coefficients[1] %>% exp(),
         beta = fit.summary$coefficients[2],
         r.squared = fit.summary$r.squared)
  }) %>% bind_rows(.id = "Population")

mut.dist.diff %>%
  ggplot(aes(`Mean (log)`, `Standard deviation (log)`, color = Treatment, group = Treatment)) +
  geom_point() +
  geom_smooth(aes(group = Treatment), method = "lm", se = FALSE) +
  #facet_wrap(~Passage %>% as.numeric()) +
  #stat_regline_equation() +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

mut.dist.diff %>%
  ggplot(aes(Passage %>% as.numeric(), `Mean`, color = Treatment)) +
  #geom_point() +
  stat_mean(aes(group = interaction(Passage, Treatment)))

mut.dist.diff.TL %>%
  ggplot(aes(Treatment, beta, color = Treatment)) + geom_boxplot() + geom_jitter()

# test if gradual is aggregating by looking at residuals

#rownames(mut.dist.diff) <- paste0(mut.dist.diff$Treatment, "_", mut.dist.diff$Population, "_", mut.dist.diff$Passage)

mut.dist.diff.lm <- lm(`Standard deviation (log)` ~ `Mean (log)`, mut.dist.diff)

mut.dist.diff$Residual <- mut.dist.diff.lm$residuals

mut.dist.diff %>%
  ggplot(aes(Treatment, Residual)) + geom_boxplot()

mut.dist.diff.residual.wilcox <- wilcox.test(mut.dist.diff %>% filter(Treatment == "gradual") %>% dplyr::select(Residual) %>% unlist(),
            mut.dist.diff %>% filter(Treatment == "sudden") %>% dplyr::select(Residual) %>% unlist(), alternative = "greater")

mut.dist.diff.residual.lm <- lm(Residual ~ Treatment + Passage, mut.dist.diff %>% dplyr::mutate(Passage := as.numeric(Passage)))

# make a figure

ggsave("6_hurst1.pdf", h.exp.gg, width = 6.85, height = 4)
ggsave("6_hurst1.png", h.exp.gg, width = 6.85, height = 4, dpi = 600)

hurst.gg <- ggarrange(h.exp.gg + labs(tag = "a"),
                      h.exp.pop.gg + theme(axis.ticks.x = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.title.y = element_blank()) + labs(tag = "b"),
                      common.legend = TRUE, widths = c(1, .5))

ggsave("hurst2.pdf", hurst.gg, width = 6.85, height = 4)
