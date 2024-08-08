source("source_me.R")

abs.diff.mtx <- list(gradual = gradual.mtx,
                     sudden = sudden.mtx) %>% lapply(function(treatment) {
                       
                       df <- lapply(treatment, function(x) {
                         lapply(1:(ncol(x) - 1), function(y) {
                           x <- x[y:(y + 1)]
                           x %>% calc.afd()
                         }) %>% unlist()
                       }) %>% bind_rows() %>% t() %>% as.data.frame()
                       
                       colnames(df) <- c("4-7", "7-10", "10-13", "13-16", "16-19", "19-22", "22-25")
                       df$rep <- row.names(df)
                       df %>% pivot_longer(-rep, names_to = "passage", values_to = "AFD")
                      
                     }) %>% bind_rows(.id = "treatment")

abs.diff.mtx$passage <- factor(abs.diff.mtx$passage, levels = c("4-7", "7-10", "10-13", "13-16", "16-19", "19-22", "22-25"))

afd.gg <- abs.diff.mtx %>% ggplot() +
  geom_boxplot(aes(passage, AFD), outlier.size = 0) +
  geom_line(aes(passage, AFD, group = rep, color = rep)) +
  facet_wrap(~treatment) + theme(legend.position = "none")

ggsave("afd.pdf", afd.gg, width = 6.85, height = 5)

# pairwise distances within passages



