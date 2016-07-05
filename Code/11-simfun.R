
auc <- function(lll, ddd){

  ROCR::performance(ROCR::prediction(lll, ddd), "auc")@y.values[[1]]

}

or <- function(sc, resp){

  sc01 <- sc > median(sc, na.rm = TRUE)
  tab <- table(sc01, resp) + .5
  log(tab[1 ,1]) + log(tab[2, 2]) - (log(tab[1, 2]) + log(tab[2, 1]))

}

createnosignaldata <- function(n, d){

  data.frame(Y = rbinom(n, 1, p = .3), X = matrix(rnorm(n * d), nrow = n))

}

selectvars <- function(data0){

  pees <- sapply(data0[, -1], function(v1){

    fit0 <- glm(data0$Y ~ v1, family = "binomial", model = FALSE, y = FALSE)
    summary(fit0)$coefficients[2, 4]

  })

  order(pees)[1:25]

}

runclassifier <- function(data0, selected){

  dat2 <- data0[, c(1, selected + 1)]
  fit1 <- glm(Y ~ ., data = dat2, family = "binomial", model = FALSE, y = FALSE)

  fit1

}

fitclassifier <- function(fit1, data1){

  scores <- predict(fit1, newdata = data1, type = "response")
  response <- as.factor(data1$Y)

  c(AUC = auc(scores, response), OR = or(scores, response))

}

naiveest <- function(n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)

  train <- dat0

  selecttr <- selectvars(train)
  fit <- runclassifier(train, selecttr)
  fitclassifier(fit, train)

}

zhutest.hold <- function(trratio = .5, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  npart.tr <- floor(trratio * n)
  tr.dex <- sample(1:n, npart.tr)
  ho.dex <- setdiff(1:n, tr.dex)

  train <- dat0[tr.dex, ]
  hold <- dat0[ho.dex, ]

  selecttr <- selectvars(dat0) ## selection performed on entire dataset
  fit <- runclassifier(train, selecttr)
  fitclassifier(fit, hold)

}

zhutest.hold2 <- function(trratio = .5, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  npart.tr <- floor(trratio * n)
  tr.dex <- sample(1:n, npart.tr)

  train <- dat0[tr.dex, ]

  selecttr <- selectvars(train) ## selection and fitt performed on holdout
  fit <- runclassifier(train, selecttr)
  fitclassifier(fit, dat0)  ## evaluated on entire dataset

}


holdoutest <- function(trratio = .5, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  npart.tr <- floor(trratio * n)
  tr.dex <- sample(1:n, npart.tr)
  ho.dex <- setdiff(1:n, tr.dex)

  train <- dat0[tr.dex, ]
  hold <- dat0[ho.dex, ]

  selecttr <- selectvars(train)
  fit <- runclassifier(train, selecttr)
  fitclassifier(fit, hold)

}


partitiondex <- function(dat0, K){

  n0 <- sum(dat0$Y == 0)
  n1 <- sum(dat0$Y == 1)

  k0 <- floor(n0 / K)
  k1 <- floor(n1 / K)

  dex.0 <- sample(which(dat0$Y == 0), n0, replace = FALSE)
  dex.1 <- sample(which(dat0$Y == 1), n1, replace = FALSE)

  i0 <- i1 <- 1
  ret <- vector("list", length = K)
  for(p in 1:K){

    ret[[p]] <- c(dex.0[i0:(i0 + k0 - 1)], dex.1[i1:(i1 + k1 - 1)])

    i0 <- i0 + k0
    i1 <- i1 + k1

  }

  left <- c(dex.0[i0:length(dex.0)], dex.1[i1:length(dex.1)])
  left <- left[sample(1:length(left), length(left), replace = FALSE)]
  Ki <- 1
  for(j in left){
    ret[[Ki]] <- c(ret[[Ki]], j)
    if(Ki < 50){
      Ki <- Ki + 1
    } else {
      Ki <- 1
    }
  }

  ret

}

zhu.cv <- function(k = 10, K = 50, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  selectin <- selectvars(dat0)

  over <- partitiondex(dat0, K = K)

  cvests <- sapply(over, function(oot.dex){

    oot <- dat0[oot.dex, ]
    estin <- dat0[setdiff(1:n, oot.dex), ]

    fit <- runclassifier(estin, selectin)
    fitclassifier(fit, oot)

  })
  rowMeans(cvests, na.rm = TRUE)

}


cvest <- function(k = 10, K = 50, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  over <- partitiondex(dat0, K = K)

  cvests <- sapply(over, function(oot.dex){

    oot <- dat0[oot.dex, ]
    estin <- dat0[setdiff(1:n, oot.dex), ]

    selectin <- selectvars(estin)
    fit <- runclassifier(estin, selectin)

    fitclassifier(fit, oot)

  })
  rowMeans(cvests, na.rm = TRUE)

}

cvpreval <- function(k = 10, K = 50, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  over <- partitiondex(dat0, K = K)

  cvests <- lapply(over, function(oot.dex){

    oot <- dat0[oot.dex, ]
    estin <- dat0[setdiff(1:n, oot.dex), ]

    selectin <- selectvars(estin)
    fit <- runclassifier(estin, selectin)

    preval <- predict(fit, newdata = oot, type = "response")
    data.frame(scores = preval, response = as.factor(oot$Y))

  })

  pval <- do.call("rbind", cvests)
  c(AUC = auc(pval$scores, pval$response), OR = or(pval$scores, pval$response))

}


bootest <- function(B = 50, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  bootests <- replicate(B, {

    boot.dex <- sample(1:n, n, replace = TRUE)
    bin <- dat0[boot.dex, ]
    notbin <- setdiff(1:n, unique(boot.dex))

    selectin <- selectvars(bin)
    fit <- runclassifier(bin, selectin)
    fitclassifier(fit, dat0[notbin, ])

  })

  rowMeans(bootests)

}

runsim <- function(n = 1000, d = 500, B = 200){

  c(holdout.5 = holdoutest(trratio = .5, n = n, d = d),
    holdout.3 = holdoutest(trratio = .666, n = n, d = d),
    cv.10 = cvest(k = 50, K = 50, n = n, d = d),
    cv.preval = cvpreval(k = 50, K = 200, n = n, d = d),
    cv.100 = cvest(k = 100, K = 10, n = n, d = d),
    boot = bootest(B = 200, n = n, d = d),
    zhu.hold = zhutest.hold(trratio = .5, n = n, d = d),
    zhu.hold2 = zhutest.hold2(trratio = .5, n = n, d = d),
    zhu.cv = zhu.cv(k = 50, K = 50, n = n, d = d),
    naive = naiveest(n = n, d = d))

}

sim.fun <- function(n.samp, d.k, n.sim, output, seed){

  res <- NULL
  for(i in 1:n.sim){

    res <- rbind(res, runsim(n = n.samp, d = d.k, B = 200))

  }

  write.table(res, output, append = TRUE, row.names = FALSE, col.names = FALSE)

}

