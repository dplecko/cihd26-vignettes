
infer_search <- function(dat, spc, solver, gt = FALSE, se = FALSE, 
                         detect_viol = FALSE, mc_cores = NULL, nbreaks = 5,
                         tau_xw = NULL) {
  
  pattern <- spc$pattern
  if (spc$type == "bound") {
    
    spc$opt_params <- parallel::mclapply(
      seq_len(nrow(spc$params)),
      function(i) bounded_optim(pattern, spc$opt_params[[i]], spc$params[i,],
                                dat, spc$method, solver),
      mc.cores = n_cores()
    )
  }
  
  # if replicating with known ground truth, fit models
  if (gt) {
    
    if (is.element("W", names(dat))) {
      
      cor_mat <- t(cbind(dat$R, dat$W)) %*% cbind(dat$R, dat$W) / nrow(dat$R)
      mu_w <- colMeans(dat$W)
      k <- ncol(dat$R)
      d <- ncol(dat$W)
      d_index <- k + seq_len(d)
      theta <- tail(binsensate:::cormat_to_Sigma_cts(cor_mat, mu_w, k, d), n = 1L)[[1]]
      
      theta_gt <- binsensate:::fit_lmb_rw(dat$X, dat$Y, dat$R, dat$W, 
                                          tau_xw = tau_xw)
      theta_gt[["theta"]] <- theta
      attr(gt, "theta_gt") <- theta_gt
    } else {
      
      Sigma0 <- binsensate:::infer_Sigma_IF(dat$X, dat$Y, dat$R, 
                                            fi = list(list(0, 0), list(0, 0)))
      lambda <- glm(X ~ ., data = data.frame(X = dat$X, R = dat$R), 
                    family = "binomial")$coef
      mu_mod <- glm(Y ~ ., data = data.frame(Y = dat$Y, X = dat$X, R = dat$R), 
                    family = "binomial")
      mu <- mu_mod$coef[-2]
      theta_gt <- list(Sigma0 = Sigma0, lambda = lambda, mu_mod = mu_mod, mu = mu)
      spc$params$ATE_gt <- NA
      attr(gt, "theta_gt") <- theta_gt
    }
    
  } else theta_gt <- NULL
  
  params <- lapply(
    seq_len(nrow(spc$params)),
    function(i) {
      # cat(i)
      spc_optim(spc$type, spc$opt_params[[i]], spc$params[i,],
                dat, spc$method, solver, gt, se, detect_viol, mc_cores,
                tau_xw)
    }#, mc.cores = n_cores()
  )
  cat("\n")
  
  spc$params <- do.call(rbind, params)
  spc
}

search_space <- function(pattern = c("agn", "x", "y"),
                         type = c("range", "bound"),
                         method = c("IF", "ZINF"),
                         fi, fi2 = NULL, k) {
  
  pattern <- match.arg(pattern, c("agn", "x", "y"))
  type <- match.arg(type, c("range", "bound"))
  method <- match.arg(method, c("IF", "ZINF"))
  
  params <- rng_to_df(pattern, fi_seq = fi, fi_seq2 = fi2)
  params$ATE <- NA
  
  if (method == "ZINF") k <- 2^k
  
  if (type == "bound") {
    
    opt_params <- lapply(
      seq_len(nrow(params)),
      function(i) {
        
        fi <- list(
          list(rep(params[i,]$fi_x0y0, k), rep(params[i,]$fi_x0y1, k)),
          list(rep(params[i,]$fi_x1y0, k), rep(params[i,]$fi_x1y1, k))
        )
      }
    )
  } else opt_params <- NULL
  
  list(
    pattern = pattern, fi = fi, type = type, method = method,
    params = params, opt_params = opt_params
  )
}

rng_to_df <- function(pattern, fi_seq, fi_seq2 = NULL) {
  
  if (is.null(fi_seq2)) fi_seq2 <- fi_seq
  
  pattern <- match.arg(pattern, c("agn", "x", "y"))
  if (pattern == "agn") {
    
    df <- data.frame(pattern = rep(pattern, length(fi_seq)))
    df$fi_x0y0 <- df$fi_x1y0 <- df$fi_x0y1 <- df$fi_x1y1 <- fi_seq
  } else if (pattern == "x") {
    
    fi_grid <- expand.grid(a = fi_seq, b = fi_seq2)
    
    df <- data.frame(pattern = rep(pattern, nrow(fi_grid)))
    df$fi_x0y0 <- df$fi_x0y1 <- fi_grid$a
    df$fi_x1y0 <- df$fi_x1y1 <- fi_grid$b
  } else if (pattern == "y") {
    
    fi_grid <- expand.grid(a = fi_seq, b = fi_seq2)
    
    df <- data.frame(pattern = rep(pattern, nrow(fi_grid)))
    df$fi_x0y0 <- df$fi_x1y0 <- fi_grid$a
    df$fi_x0y1 <- df$fi_x1y1 <- fi_grid$b
  }
  
  df
}

spc_optim <- function(type, opt_params, params, dat, method, solver, 
                      gt, se, detect_viol, mc_cores, tau_xw) {
  
  if (type == "bound") {
    
    fi <- opt_params
  } else {
    
    fi <- list(
      list(params$fi_x0y0, params$fi_x0y1),
      list(params$fi_x1y0, params$fi_x1y1)
    )
  }
  
  if (gt) {
    
    assert_that(!is.na(params$ATE[1]))
    if (is.element("W", names(dat))) {
      
      theta_gt <- attr(gt, "theta_gt")
      seed <- sample.int(10^5, size = 1)
      dati <- synth_data_mid(
        n = nrow(dat$R), k = ncol(dat$R),
        seed = seed, fi = fi, method = method,
        class = "expfam-2d-cts",
        Sigma = theta_gt$Sigma, 
        theta = theta_gt$theta,
        lam = c(theta_gt$lambda_z, theta_gt$lambda_w),
        lam_tauw = theta_gt$lambda_tauw,
        tau_xw = theta_gt$tau_xw,
        mu = c(theta_gt$mu_z, theta_gt$mu_w),
        icept_x = theta_gt$lambda_icept, icept_y = theta_gt$mu_icept,
        beta = ate_to_or_cts(
          params$ATE[1], 
          c(theta_gt$mu_icept, theta_gt$mu_z, theta_gt$mu_w), 
          dat$R, dat$W
        )
      )
    } else {
      
      # infer coefficients
      theta_gt <- attr(gt, "theta_gt")
      Sigma0 <- theta_gt$Sigma0
      lambda <- theta_gt$lambda
      mu <- theta_gt$mu
      mu_mod <- theta_gt$mu_mod
      
      seed <- sample.int(10^5, size = 1)
      dati <- synth_data_mid(
        n = nrow(dat$R), k = nrow(Sigma0),
        seed = seed, fi = fi, method = method,
        class = "expfam-2d",
        Sigma = Sigma0, lam = lambda[-1], mu = mu[-1],
        icept_x = lambda[1], icept_y = mu[1],
        beta = ate_to_or(params$ATE[1], mu_mod, dat$R)
      )
    }
  } else dati <- dat
  
  if (is.element("W", names(dat))) {
    res <- binsensate(dati$X, dati$Y, dati$R, W = dati$W, fi = fi, method = method,
                      solver = solver, se = se, detect_viol = detect_viol,
                      mc_cores = mc_cores, tau_xw = tau_xw)
  } else {
    
    res <- binsensate(dati$X, dati$Y, dati$R, W = dati$W, fi = fi, method = method,
                      solver = solver, se = se, detect_viol = detect_viol)
  }
  
  # need to run bootstrap for backward solver
  if (solver == "backward" & se) {
    
    nboot <- 500
    ate_boot <- c()
    for (nb in seq_len(nboot)) {
      
      b_idx <- sample.int(length(dati$X), replace = TRUE)
      ate_b <- binsensate(dati$X[b_idx], dati$Y[b_idx], dati$R[b_idx, ], 
                          fi = fi, method = method, solver = solver)$ATE
      ate_boot <- c(ate_boot, ate_b)
    }
    ate_sd <- sd(ate_boot)
    res$ATE_lwr <- res$ATE - 1.96 * ate_sd
    res$ATE_upr <- res$ATE + 1.96 * ate_sd
  }
  
  if (gt) {
    
    params$ATE_gt[1] <- res$ATE
    params$fi_change[1] <- res$fi_change 
    if (!is.null(res$ATE_lwr)) {
      
      params$ATE_gt_lwr[1] <- res$ATE_lwr
      params$ATE_gt_upr[1] <- res$ATE_upr
    }
  } else {
    
    params$ATE[1] <- res$ATE
    params$fi_change[1] <- res$fi_change
    if (!is.null(res$ATE_lwr)) {
      
      params$ATE_lwr[1] <- res$ATE_lwr
      params$ATE_upr[1] <- res$ATE_upr
    }
  }
  
  params
}

spc_plot <- function(spc, spc2 = NULL) {
  
  pattern <- spc$pattern
  method <- spc$method
  params <- spc$params
  params$method <- method
  if (!is.null(spc2)) {
    
    params <- rbind(params, cbind(spc2$params, method = spc2$method))
    pattern <- "multi"
  }
  
  grid_to_plt(params, pattern, method)
}

grid_to_plt <- function(lgrid, pattern, method) {
  
  is_const <- function(x) all(x == mean(x))
  
  to_model <- function(method) 
    if (method == "IF") "Independent Fidelity" else "Zero-Inflation"
  
  if (pattern == "multi") {
    
    p_multi <- ggplot(lgrid, aes(x = fi_x0y0, y = ATE, color = method)) +
      geom_point() + geom_line() +
      theme_bw() +
      xlab(latex2exp::TeX("$\\phi$")) +
      ylab("ATE Estimate") +
      scale_y_continuous(labels = scales::percent) +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.4),
        legend.box.background = element_rect()
      ) +
      scale_color_discrete(name = "Pattern", labels = c("Indep. Fidelity",
                                                        "Zero-Inflation"))
    
    return(p_multi)
  }
  
  model <- to_model(method)
  
  df <- subset(lgrid, pattern == pattern)
  txt_clr <- ifelse(df$fi_change == TRUE, "darkorange", "black")
  
  if (pattern == "agn" || 
      (pattern == "x" & (is_const(df$fi_x0y0) || is_const(df$fi_x1y0))) ||
      (pattern == "y" & (is_const(df$fi_x0y0) || is_const(df$fi_x0y1)))
  ) {
    
    if (pattern == "agn") phi <- "$\\phi$" else {
      
      if (pattern == "x" & is_const(df$fi_x0y0)) {
        
        phi <- "$\\phi_{x_1}$"
      } else phi <- "$\\phi_{x_0}$"
      
      if (pattern == "y" & is_const(df$fi_x0y0)) {
        
        phi <- "$\\phi_{y_1}$"
      } else if (pattern == "y") phi <- "$\\phi_{y_0}$"
    }
    
    if (!is.null(df$ATE_gt)) {
      
      tmp <- df$ATE
      df$ATE <- df$ATE_gt
      df$ATE_gt <- tmp
      
      if (!is.null(df$ATE_lwr)) {
        
        tmp <- df$ATE_lwr
        df$ATE_lwr <- df$ATE_gt_lwr
        df$ATE_gt_lwr <- tmp
        
        tmp <- df$ATE_upr
        df$ATE_upr <- df$ATE_gt_upr
        df$ATE_gt_upr <- tmp
      }
    }
    
    plt <- ggplot(df, aes(x = fi_x0y0, y = ATE, color = "Estimand")) + 
      geom_point() + geom_line() +
      theme_bw() +
      xlab(latex2exp::TeX(paste0(phi, " (", model, ")"))) + 
      ylab(latex2exp::TeX(
        "ATE Estimate"
        # "ATE$_{x_0, x_1}(y ; \\Phi = (\\phi, \\phi, \\phi, \\phi))$"
      )) +
      scale_y_continuous(labels = scales::percent)
    
    if (!is.null(df$ATE_gt)) {
      
      plt <- plt + 
        geom_line(aes(y = ATE_gt, color = "Ground Truth")) +
        geom_point(aes(y = ATE_gt, color = "Ground Truth")) +
        scale_color_manual(
          name = "Quantity",
          values = c("black", "red"), 
          labels = c("Estimated", "Ground Truth")
        ) +
        theme(
          legend.position = "inside", legend.position.inside = c(0.8, 0.16),
          legend.box.background = element_rect()
        )
    } else plt <- plt + scale_color_manual(values = "black") + 
      guides(color = "none")
    
    if (!is.null(df$ATE_lwr)) {
      
      plt <- plt +
        geom_ribbon(aes(ymin = ATE_lwr, ymax = ATE_upr), alpha = 0.4)
    }
    
  } else {
    
    if (pattern == "y") df$fi_x1y0 <- df$fi_x0y1
    plt <- ggplot(df, aes(x = fi_x0y0, y = fi_x1y0)) + 
      geom_tile(aes(fill = symmetric_trim(ATE, 0.02)), color = 'white') +
      scale_fill_gradient2(
        low="blue", high="red", mid="white", midpoint=0,
        limits = c(-0.02, 0.02),
        labels = scales::percent
      ) +
      theme_bw() + 
      labs(fill="ATE") +
      geom_text(aes(label=sprintf("%.2f%%", 100 * ATE), color = txt_clr), vjust=1) +
      scale_color_identity() +
      xlab(latex2exp::TeX(paste0("$\\phi_{", pattern, "_0}$", " (", model, ")"))) + 
      ylab(latex2exp::TeX(paste0("$\\phi_{", pattern, "_1}$", " (", model, ")")))
  }
  
  plt +
    theme(axis.title = element_text(size = 12), 
          axis.text = element_text(size = 10))
}

symmetric_trim <- function(x, eps) {
  
  x[x > eps] <- eps
  x[x < -eps] <- -eps
  x
}
