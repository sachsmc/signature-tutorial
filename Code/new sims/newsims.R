
library(ggplot2)

parse_sims <- function(file, path = ""){

  colnms <- c('holdout.5.AUC', 'holdout.5.OR', 'holdout.3.AUC', 'holdout.3.OR', 'cv.10.AUC', 'cv.10.OR', 'cv.preval.AUC', 'cv.preval.OR', 'cv.100.AUC', 'cv.100.OR', 'boot.AUC', 'boot.OR', 'zhu.hold.AUC', 'zhu.hold.OR', 'zhu.hold2.AUC', 'zhu.hold2.OR', 'zhu.cv.AUC', 'zhu.cv.OR', 'naive.AUC', 'naive.OR')


  titles <- c("N = 1000; p = 10",
              "N = 1000; p = 100",
              "N = 1000; p = 500",
              "N = 200; p = 10",
              "N = 200; p = 100",
              "N = 200; p = 500",
              "N = 200; p = 5000"
              )
  names(titles) <- c("results100010.txt" , "results1000100.txt",
                     "results1000500.txt", "results20010.txt" ,
                     "results200100.txt",  "results200500.txt", "results2005000.txt")
  tite <- titles[file]

  cvres <- read.table(paste0(path, file))
  colnames(cvres) <- colnms
  cvlong <- do.call("rbind", lapply(1:ncol(cvres), function(i){

    cl <- cvres[, i]
    splow <- strsplit(colnames(cvres)[i], ".", fixed = TRUE)
    nm <- paste(unlist(sapply(splow, function(s) rev(rev(s)[-1]))), collapse = ".")
    cls <- sapply(splow, function(s) rev(s)[1])
    data.frame(value = cvres[, i], stat = cls, scen = nm, stringsAsFactors = FALSE)

  }))
  cvlong$value[cvlong$stat == "OR"] <- exp(cvlong$value[cvlong$stat == "OR"])

  cvlong$scen <- factor(cvlong$scen, levels = c("naive", "zhu.cv", "zhu.hold", "zhu.hold2",
                                                "cv.preval", "cv.10", "cv.100", "holdout.3",
                                                "holdout.5", "boot"), ordered = TRUE)
  cvlong$N <- as.numeric(strsplit(tite, "( = |;)")[[1]][2])
  cvlong$p <- as.numeric(strsplit(tite, "( = |;)")[[1]][4])

  labs <- c("Resubstitution", "Partial CV",
            "Partial Holdout", "Partial Resubstitution",
            "Pre-validation", "Leave 10 out CV", "Leave 100 out CV",
            "30% Holdout", "50% Holdout", "Bootstrap")

  names(labs) <- levels(cvlong$scen)

  plot.auc <- ggplot(subset(cvlong, stat == "AUC"), aes(y = value, x = scen)) + geom_violin(fill = "grey60") +
    theme_bw(base_size = 13, base_family = "serif") + coord_flip() +
    scale_x_discrete("Estimation Approach", labels = labs) +
    ylab("Area Under the ROC Curve") + ggtitle(tite)

  plot.or <- ggplot(subset(cvlong, stat == "OR"), aes(y = value, x = scen)) + geom_violin(fill = "grey60") +
    theme_bw(base_size = 13, base_family = "serif") + coord_flip() +
    scale_x_discrete("Estimation Approach", labels = labs) + scale_y_log10(breaks = c(.5, 1, 2, 4)) +
    ylab("Odds Ratio (log scale)") + ggtitle(tite)


  cvlong$value[cvlong$stat == "OR"] <- log(cvlong$value[cvlong$stat == "OR"])

  cvlong$Approach <- factor(labs[cvlong$scen], levels = labs, ordered = TRUE)
  cvlong$Truth <- ifelse(cvlong$stat == "AUC", .5, 0.0)
  cvlong$value[!is.finite(cvlong$value)] <- NA
  cvwide <- cvlong %>% group_by(stat, Approach) %>%
    summarize(`mean` = mean(value, na.rm =TRUE), `std.dev` = sd(value, na.rm = TRUE),
              `Bias` = mean(value - Truth[1], na.rm = TRUE))

  cvwide2 <- left_join(cvwide[1:10, -1], cvwide[11:20, -1], by = c("Approach"))
  colnames(cvwide2) <- gsub(".x", " AUC", colnames(cvwide2))
  colnames(cvwide2) <- gsub(".y", " OR", colnames(cvwide2))

  table <- cvwide2 %>%
    kable(digits = 2, caption = paste("Comparison of different approaches to estimating the Area Under the ROC Curve (AUC) and the log odds ratio (OR) in the setting where a dataset is used to both develop the signature and evaluate its performance. The true value of the AUC is 0.5 and the true value of the Log OR is 0.0. Estimates are based on 1000 replicates of the numerical experiment. In each replicate, ", tite, ". CV = Cross validation. "), format = "latex")

  list(plot.auc = plot.auc, plot.or = plot.or, table = table, datalong = cvlong)

}

