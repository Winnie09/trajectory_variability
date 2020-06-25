monocle3_group <- function(expr,branch) {
  res <- mclapply(row.names(expr),function(i) {
    summary(speedglm::speedglm(expr[i,]~branch,family = stats::gaussian(),acc=1e-3, model=FALSE,y=FALSE))$coefficients[2,3:4]
  },mc.cores=detectCores())
  names(res) <- row.names(expr)
  res <- do.call(rbind,res)
  colnames(res)[2] <- 'pval'
  res$fdr <- p.adjust(res$pval,method='fdr')
  res
}

monocle3_time <- function(expr,cell_coords) {
  neighbor_graph = "knn"
  reduction_method = "UMAP"
  k = 25
  method = c('Moran_I')
  alternative = 'greater'
  cores=1
  verbose=FALSE
  expression_family='uninormal'
  
  exprs_mat <- expr
  
  my.moran.test <- function (x, listw, wc, alternative = "greater",
                             randomisation = TRUE) {
    zero.policy = TRUE
    adjust.n = TRUE
    na.action = stats::na.fail
    drop.EI2 = FALSE
    xname <- deparse(substitute(x))
    wname <- deparse(substitute(listw))
    NAOK <- deparse(substitute(na.action)) == "na.pass"
    x <- na.action(x)
    na.act <- attr(x, "na.action")
    if (!is.null(na.act)) {
      subset <- !(1:length(listw$neighbours) %in% na.act)
      listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    n <- length(listw$neighbours)
    if (n != length(x))
      stop("objects of different length")
    
    S02 <- wc$S0 * wc$S0
    res <- spdep::moran(x, listw, wc$n, wc$S0, zero.policy = zero.policy,
                        NAOK = NAOK)
    I <- res$I
    K <- res$K
    
    EI <- (-1)/wc$n1
    if (randomisation) {
      VI <- wc$n * (wc$S1 * (wc$nn - 3 * wc$n + 3) - wc$n *
                      wc$S2 + 3 * S02)
      tmp <- K * (wc$S1 * (wc$nn - wc$n) - 2 * wc$n * wc$S2 +
                    6 * S02)
      if (tmp > VI)
        warning(paste0("Kurtosis overflow,\ndistribution of variable does ",
                       "not meet test assumptions"))
      VI <- (VI - tmp)/(wc$n1 * wc$n2 * wc$n3 * S02)
      if (!drop.EI2)
        VI <- (VI - EI^2)
      if (VI < 0)
        warning(paste0("Negative variance,\ndistribution of variable does ",
                       "not meet test assumptions"))
    }
    else {
      VI <- (wc$nn * wc$S1 - wc$n * wc$S2 + 3 * S02)/(S02 *
                                                        (wc$nn - 1))
      if (!drop.EI2)
        VI <- (VI - EI^2)
      if (VI < 0)
        warning(paste0("Negative variance,\ndistribution of variable does ",
                       "not meet test assumptions"))
    }
    ZI <- (I - EI)/sqrt(VI)
    statistic <- ZI
    names(statistic) <- "Moran I statistic standard deviate"
    if (alternative == "two.sided")
      PrI <- 2 * stats::pnorm(abs(ZI), lower.tail = FALSE)
    else if (alternative == "greater")
      PrI <- stats::pnorm(ZI, lower.tail = FALSE)
    else PrI <- stats::pnorm(ZI)
    if (!is.finite(PrI) || PrI < 0 || PrI > 1)
      warning("Out-of-range p-value: reconsider test arguments")
    vec <- c(I, EI, VI)
    names(vec) <- c("Moran I statistic", "Expectation", "Variance")
    method <- paste("Moran I test under", ifelse(randomisation,
                                                 "randomisation", "normality"))
    
    res <- list(statistic = statistic, p.value = PrI, estimate = vec)
    if (!is.null(na.act))
      attr(res, "na.action") <- na.act
    class(res) <- "htest"
    res
  }
  
  my.geary.test <- function (x, listw, wc, randomisation = TRUE,
                             alternative = "greater")
  {
    zero.policy = TRUE
    adjust.n = TRUE
    spChk = NULL
    alternative <- match.arg(alternative, c("less", "greater",
                                            "two.sided"))
    if (!inherits(listw, "listw"))
      stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (!is.numeric(x))
      stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (any(is.na(x)))
      stop("NA in X")
    n <- length(listw$neighbours)
    if (n != length(x))
      stop("objects of different length")
    if (is.null(spChk))
      spChk <- spdep::get.spChkOption()
    if (spChk && !spdep::chkIDs(x, listw))
      stop("Check of data and weights ID integrity failed")
    S02 <- wc$S0 * wc$S0
    res <- spdep::geary(x, listw, wc$n, wc$n1, wc$S0, zero.policy)
    C <- res$C
    if (is.na(C))
      stop("NAs generated in geary - check zero.policy")
    K <- res$K
    EC <- 1
    if (randomisation) {
      VC <- (wc$n1 * wc$S1 * (wc$nn - 3 * n + 3 - K * wc$n1))
      VC <- VC - ((1/4) * (wc$n1 * wc$S2 * (wc$nn + 3 * n -
                                              6 - K * (wc$nn - n + 2))))
      VC <- VC + (S02 * (wc$nn - 3 - K * (wc$n1^2)))
      VC <- VC/(n * wc$n2 * wc$n3 * S02)
    }
    else {
      VC <- ((2 * wc$S1 + wc$S2) * wc$n1 - 4 * S02)/(2 * (n +
                                                            1) * S02)
    }
    ZC <- (EC - C)/sqrt(VC)
    statistic <- ZC
    names(statistic) <- "Geary C statistic standard deviate"
    PrC <- NA
    if (is.finite(ZC)) {
      if (alternative == "two.sided")
        PrC <- 2 * stats::pnorm(abs(ZC), lower.tail = FALSE)
      else if (alternative == "greater")
        PrC <- stats::pnorm(ZC, lower.tail = FALSE)
      else PrC <- stats::pnorm(ZC)
      if (!is.finite(PrC) || PrC < 0 || PrC > 1)
        warning("Out-of-range p-value: reconsider test arguments")
    }
    vec <- c(C, EC, VC)
    names(vec) <- c("Geary C statistic", "Expectation", "Variance")
    method <- paste("Geary C test under", ifelse(randomisation,
                                                 "randomisation", "normality"))
    data.name <- paste(deparse(substitute(x)), "\nweights:",
                       deparse(substitute(listw)), "\n")
    res <- list(statistic = statistic, p.value = PrC, estimate = vec,
                alternative = ifelse(alternative == "two.sided", alternative,
                                     paste("Expectation", alternative,
                                           "than statistic")),
                method = method, data.name = data.name)
    class(res) <- "htest"
    res
  }
  
  calculateLW <- function() {
    if(verbose) {
      message("retrieve the matrices for Moran's I test...")
    }
    knn_res <- NULL
    principal_g <- NULL
    
    
    if (neighbor_graph == "knn") {
      knn_res <- RANN::nn2(cell_coords, cell_coords,
                           min(k + 1, nrow(cell_coords)),
                           searchtype = "standard")[[1]]
    } else if(neighbor_graph == "principal_graph") {
      
    }
    
    jaccard_coeff <- function(idx,weight){
      weights <- matrix(0,nrow=nrow(idx)*ncol(idx), ncol=3);
      r <- 1
      for(i in 1:nrow(idx)) {
        for(j in 1:ncol(idx)) {
          k = idx[i,j] - 1
          weights[r, 1] = i + 1
          weights[r, 2] = k + 1
          weights[r, 3] = 1
          r <- r + 1
        }
      }
      weights[,2] = weights[,2] / max(weights[,2])
      weights
    }
    
    exprs_mat <- expr
    if(neighbor_graph == "knn") {
      if(is.null(knn_res)) {
        knn_res <- RANN::nn2(cell_coords, cell_coords,
                             min(k + 1, nrow(cell_coords)),
                             searchtype = "standard")[[1]]
      }
      links <- jaccard_coeff(knn_res[, -1], F)
      links <- links[links[, 1] > 0, ]
      relations <- as.data.frame(links)
      colnames(relations) <- c("from", "to", "weight")
      knn_res_graph <- igraph::graph.data.frame(relations, directed = T)
      
      knn_list <- lapply(1:nrow(knn_res), function(x) knn_res[x, -1])
      region_id_names <- colnames(expr)
      
      id_map <- 1:ncol(expr)
      names(id_map) <- id_map
      
      points_selected <- 1:nrow(knn_res)
      
      knn_list <- lapply(points_selected,
                         function(x) id_map[as.character(knn_res[x, -1])])
    } else if (neighbor_graph == "principal_graph") {
      # mapping from each cell to the principal points
      cell2pp_map <-
        cds@principal_graph_aux[[
          reduction_method]]$pr_graph_cell_proj_closest_vertex
      if(is.null(cell2pp_map)) {
        stop(paste("Error: projection matrix for each cell to principal",
                   "points doesn't exist, you may need to rerun learn_graph"))
      }
      
      # This cds object might be a subset of the one on which ordering was
      # performed, so we may need to subset the nearest vertex and low-dim
      # coordinate matrices:
      cell2pp_map <-  cell2pp_map[row.names(cell2pp_map) %in%
                                    row.names(colData(cds)),, drop=FALSE]
      cell2pp_map <- cell2pp_map[colnames(cds), ]
      
      if(verbose) {
        message("Identify connecting principal point pairs ...")
      }
      # an alternative approach to make the kNN graph based on the principal
      # graph
      knn_res <- RANN::nn2(cell_coords, cell_coords,
                           min(k + 1, nrow(cell_coords)),
                           searchtype = "standard")[[1]]
      # convert the matrix of knn graph from the cell IDs into a matrix of
      # principal points IDs
      # kNN_res_pp_map <- matrix(cell2pp_map[knn_res], ncol = k + 1, byrow = F)
      
      # kNN can be built within group of cells corresponding to each principal
      # points
      principal_g_tmp <- principal_g
      diag(principal_g_tmp) <- 1 # so set diagnol as 1
      cell_membership <- as.factor(cell2pp_map)
      uniq_member <- sort(unique(cell_membership))
      
      membership_matrix <- Matrix::sparse.model.matrix( ~ cell_membership + 0)
      colnames(membership_matrix) <- levels(uniq_member)
      # sparse matrix multiplication for calculating the feasible space
      feasible_space <- membership_matrix %*%
        Matrix::tcrossprod(principal_g_tmp[as.numeric(levels(uniq_member)),
                                           as.numeric(levels(uniq_member))],
                           membership_matrix)
      
      links <- jaccard_coeff(knn_res[, -1], F)
      links <- links[links[, 1] > 0, ]
      relations <- as.data.frame(links)
      colnames(relations) <- c("from", "to", "weight")
      knn_res_graph <- igraph::graph.data.frame(relations, directed = T)
      
      # remove edges across cells belong to two disconnected principal points
      tmp_a <- igraph::get.adjacency(knn_res_graph)
      block_size <- 10000
      num_blocks = ceiling(nrow(tmp_a) / block_size)
      if(verbose) {
        message('start calculating valid kNN graph ...')
      }
      
      tmp <- NULL
      
      for (j in 1:num_blocks){
        if (j < num_blocks){
          block_a <- tmp_a[((((j-1) * block_size)+1):(j*block_size)), ]
          block_b <- feasible_space[((((j-1) * block_size)+1):(j*block_size)), ]
        }else{
          block_a <- tmp_a[((((j-1) * block_size)+1):(nrow(tmp_a))), ]
          block_b <- feasible_space[((((j-1) * block_size)+1):(nrow(tmp_a))), ]
        }
        
        cur_tmp <- block_a * block_b
        
        if(is.null(tmp)) {
          tmp <- cur_tmp
        } else {
          tmp <- Matrix::rBind(tmp, cur_tmp)
        }
      }
      
      #close(pb_feasible_knn)
      if(verbose) {
        message('Calculating valid kNN graph, done ...')
      }
      
      region_id_names <- colnames(expr)
      
      id_map <- 1:ncol(expr)
      names(id_map) <- id_map
      
      knn_list <-
        slam::rowapply_simple_triplet_matrix(slam::as.simple_triplet_matrix(tmp),
                                             function(x) {
                                               res <- which(as.numeric(x) > 0)
                                               if(length(res) == 0)
                                                 res <- 0L
                                               res
                                             })
    } else {
      stop("Error: unrecognized neighbor_graph option")
    }
    # create the lw list for moran.test
    names(knn_list) <- id_map[names(knn_list)]
    class(knn_list) <- "nb"
    attr(knn_list, "region.id") <- region_id_names
    attr(knn_list, "call") <- match.call()
    # attr(knn_list, "type") <- "queen"
    lw <- spdep::nb2listw(knn_list, zero.policy = TRUE)
    lw
  }
  
  lw <- calculateLW()
  wc <- spdep::spweights.constants(lw, zero.policy = TRUE, adjust.n = TRUE)
  test_res <- pbmcapply::pbmclapply(row.names(exprs_mat),
                                    FUN = function(x, sz, alternative,
                                                   method, expression_family) {
                                      exprs_val <- exprs_mat[x, ]
                                      
                                      if (expression_family %in% c("uninormal", "binomialff")){
                                        exprs_val <- exprs_val
                                      }else{
                                        exprs_val <- log10(exprs_val / sz + 0.1)
                                      }
                                      
                                      test_res <- tryCatch({
                                        if(method == "Moran_I") {
                                          mt <- suppressWarnings(my.moran.test(exprs_val, lw, wc, alternative = alternative))
                                          data.frame(status = 'OK', p_value = mt$p.value,
                                                     morans_test_statistic = mt$statistic,
                                                     morans_I = mt$estimate[["Moran I statistic"]])
                                        } else if(method == 'Geary_C') {
                                          gt <- suppressWarnings(my.geary.test(exprs_val, lw, wc, alternative = alternative))
                                          data.frame(status = 'OK', p_value = gt$p.value,
                                                     geary_test_statistic = gt$statistic,
                                                     geary_C = gt$estimate[["Geary C statistic"]])
                                        }
                                      },
                                      error = function(e) {
                                        data.frame(status = 'FAIL', p_value = NA, morans_test_statistic = NA,
                                                   morans_I = NA)
                                      })
                                    }, sz = sz, alternative = alternative, method = method,
                                    expression_family = expression_family, mc.cores=cores,
                                    ignore.interactive = TRUE)
  
  
  test_res <- do.call(rbind.data.frame, test_res)
  row.names(test_res) <- row.names(expr)
  test_res$fdr <- p.adjust(test_res$p_value,method='fdr')
  test_res
}

