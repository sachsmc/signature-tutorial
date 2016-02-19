library(survival)
library(tidyr)
library(dplyr)

load("../Data/jbl10-data-2016-02-09.RData")

## supervised pca fit for heldout
##    - select x most correlated genes
##    - subset those, form 3 pcas
##    - Fit multivariate model
## cutpoint at median
## predict for held in dataset

scl.cen <- function(x){

  (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)

}

## dichotomize survival
gdat$DSS.status[gdat$DSS.status == "dead"] <- "Dead"
gdat$surv2.5 <- with(gdat, ifelse(DSS.time >= 5, 1, -1))
psdex <- 2:(grep("tissue", colnames(gdat)) - 1)
gdat[, psdex] <- lapply(gdat[, psdex], scl.cen)
ldat <- gather_(gdat, "gene", "expression", colnames(gdat[, psdex]))
ldat <- ldat %>% group_by(gene) %>% mutate(zexpression = (expression - mean(expression, na.rm = TRUE))/sd(expression, na.rm = TRUE))
cldat <- gdat[, -psdex]

fit.superpc <- function(dat0, cldat){

  cors <- dat0 %>% group_by(gene) %>%
    do(data.frame(
      t(summary(coxph(Surv(DSS.time, DSS.status == "Dead") ~ zexpression,
                           data = .,
                           eps = 1e-4))$coefficients[c(1,5)]))
      )

  ## select number of genes

  cors <- cors[unlist(cors$X2) < .005, ]
  n.gen.cands <- 5:nrow(cors)

  p.ngenes <- sapply(n.gen.cands, function(n.genes){

    sel.genes <- cors$gene[rank(-abs(cors$X1)) <= n.genes]
    coeffs <- cors$X1[rank(-abs(cors$X1)) <= n.genes]

    pcas <- dat0 %>% select(ID, gene, zexpression) %>%
      filter(gene %in% sel.genes) %>% spread(gene, zexpression)
    pcamat <- as.matrix(pcas[, -1]) * matrix(rep(coeffs, nrow(pcas)), nrow = nrow(pcas), byrow = TRUE)
    rownames(pcamat) <- pcas$ID

    pc.fit <- prcomp(pcamat, center = TRUE, scale. = TRUE)
    pc.dat <- as.data.frame(pc.fit$x[, 1:5])
    pc.dat$ID <- rownames(pc.dat)

    pdonk <- merge(cldat, pc.dat, by = "ID")
    fit.cox <- coxph(Surv(DSS.time, DSS.status == "Dead") ~ PC1 + PC2 + PC3 + PC4 + PC5, data = pdonk)
    summary(fit.cox)$concordance[1]

  })

  n.gene.sel <- n.gen.cands[which(order(-abs(p.ngenes)) == 1)]
  sel.genes <- cors$gene[rank(-abs(cors$X1)) <= n.gene.sel]

  pcas <- dat0 %>% select(ID, gene, zexpression) %>%
    filter(gene %in% sel.genes) %>% spread(gene, zexpression)
  pcamat <- as.matrix(pcas[, -1])
  rownames(pcamat) <- pcas$ID

  pc.fit <- prcomp(pcamat, center = FALSE, scale. = FALSE)
  pc.dat <- as.data.frame(pc.fit$x[, 1:5])
  pc.dat$ID <- rownames(pc.dat)

  pdonk <- merge(cldat, pc.dat, by = "ID")
  fit.cox <- coxph(Surv(DSS.time, DSS.status == "Dead") ~ PC1 + PC2 + PC3 + PC4 + PC5, data = pdonk)

  lps <- predict(fit.cox, type = 'lp')

  ## middle 75%
  c.cands <- sort(lps)[floor(.25 * length(lps)):floor(.75 * length(lps))]
  p.cands <- sapply(c.cands, function(c.cand){
    ltest <- survdiff(Surv(pdonk$DSS.time, pdonk$DSS.status == "Dead") ~ I(lps > c.cand))
    pchisq(ltest$chisq, df = 1, lower.tail = FALSE)
  })
  cutoff <- c.cands[which(order(p.cands) == 1)]
  list(pc.fit = pc.fit, cox.fit = fit.cox, cutoff = cutoff)

}


predict.superpc <- function(dat0, las.fit){

  pcests <- predict(las.fit$pc.fit, newdata = dat0)[, 1:5]
  dat0 <- cbind(dat0, pcests)
  coxests <- predict(las.fit$cox.fit, newdata = dat0, type = 'lp')
  list(lps = coxests, riskgrp = coxests > las.fit$cutoff)


}

eval.predict <- function(dat0, las.fit){

  lpreds <- predict.superpc(dat0, las.fit)
  sc <- survConcordance(Surv(dat0$DSS.time, dat0$DSS.status == "Dead") ~ lpreds$lps)

  disc <- coxph(Surv(dat0$DSS.time, dat0$DSS.status == "Dead") ~ lpreds$riskgrp)$coefficients

  c(calib = sc$concordance, discrim = disc)

}


## substitution estimate using full dataset

fit.all <- fit.superpc(subset(ldat, Post.Surgical.Treatment == "OBS"), cldat)
#fit.obs <- fit.superpc(ldat, cldat)
fitted <- predict.superpc(subset(gdat, Post.Surgical.Treatment == "OBS"), fit.all)
eval.predict(subset(gdat, Post.Surgical.Treatment == "OBS"), fit.all)

## partial resubstitution






