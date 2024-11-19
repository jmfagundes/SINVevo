library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(readxl)
library(segmented)
library(ggpubr)
library(viridis)
library(ggrepel)
#library(randtests)
library(igraph)
library(minpack.lm)
library(pracma)
library(statcomp)
library(grid)
library(gridExtra)
library(ggbreak)
library(marginaleffects)

library(Hmisc)
library(cluster)
library(factoextra)

# load matrices

gradual.mtx <- lapply(setNames(1:9, nm = c(98, 101, 102,
                                           106, 107, 109,
                                           112, 115, 116) %>% as.character()), function(x) {
                                             y <- read_xlsx("data/Gradual.xlsx", sheet = x) %>% as.data.frame()
                                             rownames(y) <- y$'...1'
                                             y[-c(1,2)]
                                           })

sudden.mtx <- lapply(setNames(1:9, nm = c(2, 5, 6,
                                          10, 11, 13,
                                          16, 19, 20) %>% as.character()), function(x) {
                                            y <- read_xlsx("data/Sudden.xlsx", sheet = x) %>% as.data.frame()
                                            rownames(y) <- y$'...1'
                                            y[-c(1,2)]
                                          })

# get approxwf results

approxwf.res <- list(gradual = lapply(setNames(nm = c(98, 101, 102,
                                                      106, 107, 109,
                                                      112, 115, 116) %>% as.character()), function(x) {
                                                        y <- read.table(paste0("SE_results/Results/Gradual.", x, "_res.csv"),
                                                                        sep = ",", header = TRUE)
                                                        y$allele <- rownames(gradual.mtx[[as.character(x)]])
                                                        y
                                                      }) %>% bind_rows(.id = "population"),
                     sudden = lapply(setNames(nm = c(2, 5, 6,
                                                     10, 11, 13,
                                                     16, 19, 20) %>% as.character()), function(x) {
                                                       y <- read.table(paste0("SE_results/Results/Sudden.", x, "_res.csv"),
                                                                       sep = ",", header = TRUE)
                                                       y$allele <- rownames(sudden.mtx[[as.character(x)]])
                                                       y
                                                     }) %>% bind_rows(.id = "population")) %>%
  bind_rows(.id = "treatment")

# some ggplot themes

blank.theme <- theme(strip.text = element_text(hjust = 0, size = 12))

# functions

# return allele abundances

#' normalize.cells
#' 
#' Calculate allele abundances
#' 
#' @param mtx A matrix where rows are alleles and columns are cells
#' 
#' @return output allele abundances matrix
#' @export
#' 
#' @examples
normalize.cells <- function(mtx) {
  mtx <- mtx[colSums(mtx) > 0]
  t(t(mtx)/colSums(mtx))
}

#' fit.TL
#' 
#' Fit matrices in a list to Taylor's law
#' 
#' @param matrices List of matrices
#' @param prefix.xlsx Prefix of .xlsx files for cmplxcruncher. If NULL, don't save xlsx files
#' @param zero.rate.threshold Filter out alleles with percentage of zeros higher than zero.rate.threshold
#' @param normalize Boolean. Whether to transform the matrices to abundance matrix by performing cell size normalization
#' @param remove.zeros Do not include zeros when calculating mean and standard deviation
#' @param analyze.fluctuation.only Only analyze the fluctuation of mutations that reach fixation (removes all zeros and ones then add one zero and one one)
#' @param min.rows Skip matrix with number of elements (rows) < min.rows
#' @param sd.from.max Calculate max ~ sd from max
#' 
#' @return A list of the Taylor's parameters, the linear models and the log-transformed mean and standard deviation data
#' @export
#' 
#' @examples
fit.TL <- function(matrices,
                   prefix.xlsx = NULL,
                   zero.rate.threshold = .95,
                   normalize = TRUE,
                   remove.zeros = FALSE,
                   analyze.fluctuation.only = FALSE,
                   min.rows = NULL,
                   sd.from.max = FALSE) {
  
  params.tb <- data.frame(data = character(),
                          V = numeric(),
                          beta = numeric(),
                          R.squared = numeric(),
                          n.cells = numeric(),
                          n.alleles = numeric(),
                          sparsity = numeric(),
                          model = character())
  log_log.mean_sd.lst <- list()
  lm.lst <- list()

  for (m in matrices %>% names()) {
    
    allele.mtx <- matrices[[m]]
    
    if (is.null(ncol(allele.mtx))) next
    if (!is.null(min.rows)) if (nrow(allele.mtx) < min.rows) next
    
    if (!is.null(prefix.xlsx)) {
      
      # if too many cells, save three different subset files of 5000 cells
      
      if (ncol(allele.mtx) > 10000) {
        
        set.seed(12345)
        rep <- allele.mtx[sample(ncol(allele.mtx), 5000)]
        write_xlsx(cbind(" " = rownames(rep), rep), paste0("cmplxcruncher_analysis/", prefix.xlsx, "_", m, "_rep_1", ".xlsx"))
        
        set.seed(1000)
        rep <- allele.mtx[sample(ncol(allele.mtx), 5000)]
        write_xlsx(cbind(" " = rownames(rep), rep), paste0("cmplxcruncher_analysis/", prefix.xlsx, "_", m, "_rep_2", ".xlsx"))
        
        set.seed(255)
        rep <- allele.mtx[sample(ncol(allele.mtx), 5000)]
        write_xlsx(cbind(" " = rownames(rep), rep), paste0("cmplxcruncher_analysis/", prefix.xlsx, "_", m, "_rep_3", ".xlsx"))
        
      } else {
        write_xlsx(cbind(" " = rownames(allele.mtx), allele.mtx), paste0("cmplxcruncher_analysis/", prefix.xlsx, "_", m, ".xlsx"))
      }
    }
    if (!is.null(zero.rate.threshold)) {
      
      zero.rate.alleles <- calc.0_rate(list(mtx = allele.mtx))$mtx
      zero.rate.alleles <- zero.rate.alleles[zero.rate.alleles > zero.rate.threshold] %>% names()
      
      allele.mtx <- allele.mtx[!rownames(allele.mtx) %in% zero.rate.alleles,]
    }
    if (normalize) allele.mtx <- allele.mtx %>% normalize.cells()
    
    if (remove.zeros) {
      mean_sd <- data.frame(mean = apply(allele.mtx, 1, function(x) {
        if (sd.from.max) max(x[x > 0])
        else mean(x[x > 0])
      }),
      sd = apply(allele.mtx, 1, function(x) {
        if (sd.from.max) sqrt(sum((x[x > 0] - max(x))**2)/(length(x) - 1))
        else sd(x[x > 0])
      }))
      mean_sd <- mean_sd[mean_sd$mean > 0 & mean_sd$sd > 0,]
      
    } else if (analyze.fluctuation.only){
      mean_sd <- data.frame(mean = apply(allele.mtx, 1, function(x) {
        mean(c(x[x > 0 & x < 1], 0 ,1))
      }),
      sd = apply(allele.mtx, 1, function(x) {
        sd(c(x[x > 0 & x < 1], 0 ,1))
      }))
      mean_sd <- mean_sd[mean_sd$mean > 0 & mean_sd$sd > 0,]
      
    } else {
      mean_sd <- data.frame(mean = apply(allele.mtx, 1, function(x) {
        if (sd.from.max) max(x)
        else mean(x)
      }),
      sd = apply(allele.mtx, 1, function(x) {
        if (sd.from.max) sqrt(sum((x - max(x))**2)/(length(x) - 1))
        else sd(x)
      }))
      mean_sd <- mean_sd[mean_sd$mean > 0 & mean_sd$sd > 0,]
    }
    
    log_log.mean_sd <- log(mean_sd)
    lm.log_log.mean_sd <- lm(sd ~ mean, data = log_log.mean_sd)
    fit.summary <- summary(lm.log_log.mean_sd)
    
    params.tb[nrow(params.tb) + 1,] <- list(m, fit.summary$coefficients[1] %>% exp(),
                                            fit.summary$coefficients[2],
                                            fit.summary$r.squared,
                                            allele.mtx %>% ncol(),
                                            log_log.mean_sd %>% nrow(),
                                            coop::sparsity(allele.mtx %>% as.matrix()),
                                            "LLR")
    log_log.mean_sd.lst[[m]] <- log_log.mean_sd
    lm.lst[[m]] <- lm.log_log.mean_sd
  }
  return(list(params = params.tb,
              log_log.mean_sd = log_log.mean_sd.lst,
              lm.log_log.mean_sd = lm.lst))
}

#' calc.0_rate
#' 
#' Calculate the zero rate of alleles
#' 
#' @param matrices List of allele matrices
#' 
#' @return A named list with the zero rate of each allele
#' @export
#' 
#' @examples
calc.0_rate <- function(matrices) {
  
  zero_rate.lst <- list()
  
  for (m in names(matrices)) {
    
    allele.mtx <- matrices[[m]]
    zero_rate <- rowSums(allele.mtx == 0) / ncol(allele.mtx)
    
    zero_rate.lst[[m]] <- zero_rate
  }
  return(zero_rate.lst)
}

#' plot.power_law
#' 
#' Plot Taylor's plot with metadata 
#' 
#' @param log_log.mean_sd.matrices List of means and standard deviations in log-log
#' @param alleles.metadata.lst List containing, for each log-log mean-standard deviation matrix, a named list where names are the alleles
#' @param legend.title Legend title
#' @param metadata.range Range on legend. If NULL, assume values on alleles.metadata.lst are discrete
#' @param l10 If input is natural log, transform to log10
#' 
#' @return A ggplot
#' @export
#' 
#' @examples
plot.power_law <- function(log_log.mean_sd.matrices,
                           alleles.metadata.lst = NULL,
                           legend.title = NULL,
                           metadata.range = c(0, 1),
                           l10 = TRUE) {
  
  log_log.df <- log_log.mean_sd.matrices %>% bind_rows(.id = "title")
  log_log.df$alleles <- lapply(1:length(log_log.mean_sd.matrices),
                             function(x) log_log.mean_sd.matrices[[x]] %>% rownames()) %>% unlist()
  log_log.df$title <- factor(log_log.df$title, levels = names(log_log.mean_sd.matrices))
  
  if (l10) log_log.df <- log_log.df %>% mutate(mean = log10(exp(mean)), sd = log10(exp(sd)))
  
  if (is.null(alleles.metadata.lst)) {
    
    gg <- ggplot(log_log.df, aes(mean, sd))
    
  } else {
    
    value <- stack(alleles.metadata.lst)
    value$alleles <-  alleles.metadata.lst %>% lapply(function(x) names(x)) %>% unlist()
    
    log_log.df$value <- value[match(paste0(log_log.df$title, log_log.df$alleles),
                                    paste0(value$ind, value$alleles)), "values"]
    
    if (is.null(metadata.range)) log_log.df$value <- as.factor(log_log.df$value)
    
    gg <- ggplot(log_log.df, aes(mean, sd, color = value))
  }
  
  gg <- gg + facet_wrap(~title) +
    geom_point(shape = 20, stroke = 0, size = 1.2)
  
  if (!is.null(legend.title)) gg <- gg + labs(color = legend.title)
  else gg <- gg + theme(legend.title = element_blank())
  
  if (!is.null(metadata.range)) gg <- gg + scale_color_gradient(limits = metadata.range)
  else gg <- gg + scale_color_discrete(na.translate = FALSE)
  return(gg)
}

#' fit.hurst
#' 
#' Estimate H exponent for a list of matrices
#' 
#' @param matrices List of matrices
#' @param gene.lst Only estimate H for these genes
#' @param d Minimum window size for calculating H
#' 
#' @return A table with the results
fit.hurst <- function(matrices,
                      gene.lst = NULL,
                      d = 50) {
  
  hurst.fit <- list()
  
  if (is.null(gene.lst)) {
    gene.lst <- lapply(matrices, rownames) %>% unlist() %>% unique() %>% as.character()
  }
  
  for (gene in gene.lst) {
    H.lst <- mapply(function(x, y) {
      
      if (!gene %in% rownames(x)) {
        message(paste0("gene ", gene, " not found in matrix ", y, ", returning NA"))
        return(list(data = y, Hrs = NA, Ht = NA, Hal = NA, He = NA))
        
      } else c(data = y, hurstexp(x[gene,] %>% unlist(), display = FALSE, d))
    }, matrices, names(matrices), SIMPLIFY = FALSE) %>% bind_rows()
    hurst.fit[[gene]] <- H.lst
  }
  return(hurst.fit)
}

# filters for reconstructing haplotypes

# only include SNPs present in at least n.time points
n.time <- 1

# filter SNPs that never reach freq.thresh
freq.thresh <- .03

#' apply.filter
#' 
#' Filters a frequency matrix
#' 
#' @param x Input matrix
#' @param n.time Only include SNPs present in at least n.time times
#' @param freq.thresh Filter SNPs that never reach freq.thresh
#' @param rev Return SNPs below freq.thresh. Ignores n.time
#' @param exact Boolean. Only include SNPs present exactly n.time times
#' 
#' @return A filtered matrix
#' @export
#' 
#' @examples
apply.filter <- function(x,
                         n.time = 1,
                         freq.thresh = 0,
                         rev = FALSE,
                         exact = FALSE) {
  
  if (!exact) above_n.time <- apply(x, 1, function(y) length(y[y > 0]) >= n.time)
  else above_n.time <- apply(x, 1, function(y) length(y[y > 0]) == n.time)
  
  above_freq.thresh <- apply(x, 1, function(y) length(y[y >= freq.thresh]) > 0)
  
  if (!rev) return(x[above_n.time & above_freq.thresh,])
  else return(x[!above_freq.thresh,])
}

#' shannon.entropy
shannon.entropy <- function(freq.lst) -sum(lapply(freq.lst, function(x) x * log(x, 2)) %>% unlist())

shannon.entropy.lst <- function(obs) {
  freq <- table(obs) / length(obs)
  vec <- as.data.frame(freq)[,2]
  vec < -vec[vec > 0]
  -sum(vec * log2(vec))
}

#' calc.entropy
#' 
#' Calculate entropies for each allele frequency
#' 
#' @param alt.freq.list List of alternative allele frequencies
#' @param return.sum Return the sum of the entropies
#' 
#' @return Shannon entropy sum
#' @export
#' 
#' @examples
calc.entropy <- function(alt.freq.lst,
                         return.sum = TRUE) {
  
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
  
  if (return.sum) sum(entropies %>% unlist()) %>% return()
  else return(entropies)
}

#' calc.afd
#' 
#' Calculate AFD between two populations
#' 
#' @param AF.tb Two-column data frame with allele frequencies
#' 
#' @return AFD
#' @export
#' 
#' @examples
calc.afd <- function(AF.tb) {
  
  # remove indels
  
  AF.tb <- AF.tb[!grepl("\\+|\\-", row.names(AF.tb)),]

  # get snps at same site
  
  sites <- row.names(AF.tb) %>%
    gsub("[A-Z+-]*", "", .) %>% gsub("\\*", "", .) %>% unique()
  
  # calcute afd
  
  afd <- lapply(setNames(nm = sites), function(x) {
    
    alleles <- AF.tb[grepl(paste0("^", x, "[ACTG]|^", x, "\\*"), row.names(AF.tb)),]
    alleles["ref",] <- 1 - colSums(alleles)
    alleles.diff <- abs(alleles[1] - alleles[2])
    colSums(alleles.diff) / 2
    
  }) %>% unlist() %>% sum()
  
  return(afd)
}

#' fit.time.sd
#' 
#' Calculate deviation from next/previous point to investigate systemic deterministic behavior
#' 
#' @param matrices Input matrices
#' @param next.point Use deviation from next point
#' @param return.mean Returns the mean value of each allele
#' @param not.fix Do not analyse differences that lead to the loss/fixation of an allele
#' @param fix Only analyse alleles that are fixed
#' @param use.entropy Use the entropy of an allele (summing the frequencies of all other alternative allele)
#' @param use.quasi.entropy Use -freq * log(freq, 2)
#' @param remove.indels Remove indels. Set to TRUE if use.entropy = TRUE
#' @param logt Power-law relationship
#' 
#' @return A list with the results with same names as fit.TL for compatibility
#' @export
#' 
#' @examples
fit.time.sd <- function(matrices,
                        next.point = TRUE,
                        return.mean = TRUE,
                        not.fix = FALSE,
                        fix = FALSE,
                        use.entropy = FALSE,
                        use.quasi.entropy = FALSE,
                        remove.indels = FALSE,
                        logt = TRUE) {
  
  params.tb <- data.frame(data = character(),
                          V = numeric(),
                          beta = numeric(),
                          R.squared = numeric(),
                          n.cells = numeric(),
                          n.alleles = numeric(),
                          sparsity = numeric(),
                          model = character())
  log_log.mean_sd.lst <- list()
  lm.lst <- list()
  
  for (m in names(matrices)) {
    
    x <- matrices[[m]]
    
    if (fix) x <- x %>% apply.filter(freq.thresh = 1)

    if (nrow(x) == 0) next
    
    if (use.entropy | remove.indels) x <- x[!grepl("\\+|\\-", rownames(x)),]
    
    if (not.fix) x[x == 0 | x == 1] <- NA 

    d.x <- apply(x, 1, diff) %>% abs() %>% t() %>% as.data.frame()
    
    if (use.entropy) x <- x %>% mutate_all(function(y) {
      
      lapply(y, function(z) {
        if (is.na(z)) return(z)
        shannon.entropy(c(z, 1 - z))
      })
    })
    
    if (use.quasi.entropy) x <- x %>% mutate_all(function(y) {
      
      lapply(y, function(z) {
        if (is.na(z)) return(z)
        -z * log(z, 2)
      })
    })

    if (return.mean) {
      
      mean_sd <- data.frame(mean = rowMeans(x),
                            sd = rowMeans(d.x ** 2))
        
      mean_sd <- mean_sd[mean_sd$mean > 0 & mean_sd$sd > 0,]
      
      log_log.mean_sd <- log(mean_sd)
      lm.log_log.mean_sd <- lm(sd ~ mean, data = log_log.mean_sd)
      fit.summary <- summary(lm.log_log.mean_sd)
      
      
      params.tb[nrow(params.tb) + 1,] <- list(m, fit.summary$coefficients[1] %>% exp(),
                                             fit.summary$coefficients[2],
                                             fit.summary$r.squared,
                                             x %>% ncol(),
                                             log_log.mean_sd %>% nrow(),
                                             coop::sparsity(x %>% as.matrix()),
                                             "LLR")
      log_log.mean_sd.lst[[m]] <- log_log.mean_sd
      lm.lst[[m]] <- lm.log_log.mean_sd
    
    } else {
      
      # keep mean as value for compatibility with other functions
      
      if (next.point) x <- x[1:(length(x) - 1)]
      else x <- x[2:length(x)]

      v_sd <- data.frame(mean = stack(x)$values,
                         sd = stack(d.x)$values,
                         name = rownames(x))
      
      v_sd <- v_sd[v_sd$mean > 0 & v_sd$sd > 0 & !is.na(v_sd$sd),]
      
      if (logt) v_sd <- log(v_sd[c("mean", "sd")]) %>% cbind(v_sd["name"])
      
      lm.fit <- lm(sd ~ mean, data = v_sd)
      fit.summary <- summary(lm.fit)
      
      params.tb[nrow(params.tb) + 1,] <- list(m, fit.summary$coefficients[1] %>% exp(),
                                              fit.summary$coefficients[2],
                                              fit.summary$r.squared,
                                              x %>% ncol(),
                                              v_sd %>% nrow(),
                                              0,#coop::sparsity(x %>% as.matrix()),
                                              "LLR")
      log_log.mean_sd.lst[[m]] <- v_sd
      lm.lst[[m]] <- lm.fit
      
    }  
  }
  return(list(params = params.tb,
              log_log.mean_sd = log_log.mean_sd.lst,
              lm.log_log.mean_sd = lm.lst))
}
