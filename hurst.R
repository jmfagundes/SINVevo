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

mut.dist.diff.tb <-  list(gradual = gradual.mtx,
                          sudden = sudden.mtx) %>% lapply(function(x) {
                            
                            lapply(x, function(y) {
                              
                              alleles <- rownames(y)
                              alleles <- alleles[!grepl("\\+|\\-", alleles)]
                              
                              apply(y, 2, function(z) {
                                
                                z <- z[!grepl("\\+|\\-", names(z))]
                                dists <- alleles[z[!grepl("\\+|\\-", names(z))] > 0] %>% gsub("\\D+", "", .) %>% as.numeric() %>% diff()
                                data.frame(Distances = dists)
                              }) %>% bind_rows(.id = "Passage")
                            }) %>% bind_rows(.id = "Population")
                          }) %>% bind_rows(.id = "Treatment")

mut.dist.diff.tb

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
                                                                            `Standard deviation (log)` = log10(`Standard deviation`),
                                                                            Passage := Passage %>% factor(levels = c(Passage %>% unique())))

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

mut.dist.diff.TL_by_t <- mut.dist.diff$Passage %>% unique() %>%
  lapply(function(x) {
    
    df <- mut.dist.diff %>% dplyr::filter(Passage == x)
    
    lapply(setNames(nm = c("gradual", "sudden")), function(y) {
      
      df.treat <- df %>% dplyr::filter(Treatment == y)
      
      lm <- lm.log_log.mean_sd <- lm(`Standard deviation (log)` ~ `Mean (log)`, data = df.treat)
      fit.summary <- summary(lm.log_log.mean_sd)
      
      list(Treatment = y %>%
             gsub("^g", "G", .) %>%
             gsub("^s", "S", .),
           V = fit.summary$coefficients[1] %>% exp(),
           beta = fit.summary$coefficients[2],
           r.squared = fit.summary$r.squared)
    }) %>% bind_rows()
    
  }) %>% bind_rows(.id = "Passage")

mut.dist.diff %>%
  ggplot(aes(`Mean (log)`, `Standard deviation (log)`, color = Treatment, group = Treatment)) +
  geom_point() +
  geom_smooth(aes(group = Treatment), method = "lm", se = FALSE) +
  #facet_wrap(~Passage %>% as.numeric()) +
  #stat_regline_equation() +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

mut.dist.diff %>%
  ggplot(aes(`Mean (log)`, `Standard deviation (log)`, group = Passage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Treatment)

mut.dist.diff %>%
  ggplot(aes(`Mean (log)`, `Standard deviation (log)`, group = Population)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Treatment)

mut.dist.diff %>%
  ggplot(aes(Passage %>% as.numeric(), `Mean (mean)`, color = Treatment)) +
  #geom_point() +
  stat_mean(aes(group = interaction(Passage, Treatment)))

mut.dist.diff %>%
  ggplot(aes(Passage %>% as.numeric(), `Standard deviation (mean)`, color = Treatment)) +
  #geom_point() +
  stat_mean(aes(group = interaction(Passage, Treatment)))

mut.dist.diff.TL %>%
  ggplot(aes(Treatment, beta, color = Treatment)) + geom_boxplot() + geom_jitter()

mut.dist.diff.TL_by_t %>%
  ggplot(aes(Treatment, beta, color = Treatment)) + geom_boxplot() + geom_jitter()

# SD / mean

mut.dist.diff %>% mutate(`SD / mean` = `Standard deviation` / Mean) %>%
  ggplot(aes(Passage, `SD / mean`, color = Treatment)) + geom_boxplot()

sd_mean.norm.lm <- mut.dist.diff %>% mutate(`SD / mean` = `Standard deviation` / Mean) %>%
  lm(`SD / mean` ~ Treatment + Population + Passage, .)

# test if gradual is aggregating by looking at residuals

#rownames(mut.dist.diff) <- paste0(mut.dist.diff$Treatment, "_", mut.dist.diff$Population, "_", mut.dist.diff$Passage)

mut.dist.diff.lm <- lm(`Standard deviation (log)` ~ `Mean (log)`, mut.dist.diff)

mut.dist.diff$Residual <- mut.dist.diff.lm$residuals

mut.dist.diff %>%
  ggplot(aes(Treatment, Residual)) + geom_boxplot()

mut.dist.diff.residual.wilcox <- wilcox.test(mut.dist.diff %>% filter(Treatment == "gradual") %>% dplyr::select(Residual) %>% unlist(),
            mut.dist.diff %>% filter(Treatment == "sudden") %>% dplyr::select(Residual) %>% unlist(), alternative = "greater")

mut.dist.diff.residual.lm <- lm(Residual ~ Treatment + Passage, mut.dist.diff %>% dplyr::mutate(Passage := as.numeric(Passage)))

# test fit by passage

wilcox.test(mut.dist.diff.TL_by_t %>% filter(Treatment == "Gradual") %>% dplyr::select(beta) %>% unlist(),
            mut.dist.diff.TL_by_t %>% filter(Treatment == "Sudden") %>% dplyr::select(beta) %>% unlist(), paired = TRUE)

# make a figure

ggsave("7_hurst1.pdf", h.exp.gg, width = 6.85, height = 4)
ggsave("7_hurst1.png", h.exp.gg, width = 6.85, height = 4, dpi = 600)
