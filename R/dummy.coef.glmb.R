#' Extract Coefficients in Original Coding
#'
#' This extracts coefficients in terms of the original levels of the coefficients
#' rather than the coded variable
#' 
#' @param object a \code{glmb} model fit
#' @param use.na logical flag for coefficients in a singular model. If use.na is true, 
#' undetermined coefficients will be missing; if false they will get one possible value.
#' @param x object to be printed
#' @param \ldots arguments passed to or from other methods
#' @details A fitted linear model has coefficients for the contrasts of the factor terms, 
#' usually one less in number than the number of levels. This function re-expresses the 
#' coefficients in the original coding; as the coefficients will have been fitted in the 
#' reduced basis, any implied constraints (e.g., zero sum for contr.helmert or contr.sum) 
#' will be respected. There will be little point in using dummy.coef for contr.treatment 
#' contrasts, as the missing coefficients are by definition zero. 
#' @return A list giving for each term the draws for the coefficients. 
#' @example inst/examples/Ex_confint.glmb.R
#' @exportClass  dummy.coef.glmb
#' @method dummy.coef glmb
#' @export


dummy.coef.glmb<-function (object, use.na = FALSE, ...) 
{
  ## Top part added
  
  object_means=object
  object_mode=object
  object_means$coefficients=object$coef.means
  object_mode$coefficients=object$coef.mode
  
  coef.means=dummy.coef.lm(object_means)
  coef.mode=dummy.coef.lm(object_mode)
  
    ## Note only one change from dummy.coef.lm
  ## The matrix with coefficients in glmb
  ## has coefficients across and draws coming down
  ## while dummy.coef.glm has ability to handle reverse
  
  xl <- object$xlevels
  if (!length(xl)) 
    return(as.list(coef(object)))
  Terms <- terms(object)
  tl <- attr(Terms, "term.labels")
  int <- attr(Terms, "intercept")
  facs <- attr(Terms, "factors")[-1, , drop = FALSE]
  Terms <- delete.response(Terms)
  mf <- object$model
  if (is.null(mf)) 
    mf <- model.frame(object)
  vars <- dimnames(facs)[[1]]
  xtlv <- lapply(mf[, vars, drop = FALSE], levels)
  nxl <- pmax(lengths(xtlv), 1L)
  lterms <- apply(facs, 2L, function(x) prod(nxl[x > 0]))
  nl <- sum(lterms)
  args <- sapply(vars, function(i) if (nxl[i] == 1) 
    rep.int(1, nl)
    else factor(rep.int(xtlv[[i]][1L], nl), levels = xtlv[[i]]), 
    simplify = FALSE)
  dummy <- do.call(data.frame, args)
  names(dummy) <- vars
  pos <- 0L
  rn <- rep.int(tl, lterms)
  rnn <- character(nl)
  for (j in tl) {
    i <- vars[facs[, j] > 0]
    ifac <- i[nxl[i] > 1]
    lt.j <- lterms[[j]]
    if (length(ifac) == 0L) {
      rnn[pos + 1L] <- j
    }
    else {
      p.j <- pos + seq_len(lt.j)
      if (length(ifac) == 1L) {
        dummy[p.j, ifac] <- x.i <- xtlv[[ifac]]
        rnn[p.j] <- as.character(x.i)
      }
      else {
        tmp <- expand.grid(xtlv[ifac], KEEP.OUT.ATTRS = FALSE)
        dummy[p.j, ifac] <- tmp
        rnn[p.j] <- apply(as.matrix(tmp), 1L, paste, 
                          collapse = ":")
      }
    }
    pos <- pos + lt.j
  }
  attr(dummy, "terms") <- attr(mf, "terms")
  lcontr <- object$contrasts
  lci <- vapply(dummy, is.factor, NA)
  lcontr <- lcontr[names(lci)[lci]]
  mm <- model.matrix(Terms, dummy, lcontr, xl)
  if (anyNA(mm)) {
    warning("some terms will have NAs due to the limits of the method")
    mm[is.na(mm)] <- NA
  }
  
  ## KN Edited this to be transpose
  coef <- t(object$coefficients)
  if (!use.na) 
    coef[is.na(coef)] <- 0
  asgn <- attr(mm, "assign")
  res <- setNames(vector("list", length(tl)), tl)
  if (isM <- is.matrix(coef)) {
    for (j in seq_along(tl)) {
      keep <- which(asgn == j)
      cf <- coef[keep, , drop = FALSE]
      ij <- rn == tl[j]
      cf <- if (any(na <- is.na(cf))) {
        if (ncol(cf) >= 2) 
          stop("multivariate case with missing coefficients is not yet implemented")
        rj <- t(mm[ij, keep[!na], drop = FALSE] %*% cf[!na])
        rj[apply(mm[ij, keep[na], drop = FALSE] != 0, 
                 1L, any)] <- NA
        rj
      }
      else t(mm[ij, keep, drop = FALSE] %*% cf)
      dimnames(cf) <- list(colnames(coef), rnn[ij])
      res[[j]] <- cf
    }
  }
  else {
    for (j in seq_along(tl)) {
      keep <- which(asgn == j)
      cf <- coef[keep]
      ij <- rn == tl[j]
      res[[j]] <- if (any(na <- is.na(cf))) {
        rj <- setNames(drop(mm[ij, keep[!na], drop = FALSE] %*% 
                              cf[!na]), rnn[ij])
        rj[apply(mm[ij, keep[na], drop = FALSE] != 0, 
                 1L, any)] <- NA
        rj
      }
      else setNames(drop(mm[ij, keep, drop = FALSE] %*% 
                           cf), rnn[ij])
    }
  }
  if (int > 0)     res <- c(list(`(Intercept)` = if (isM) coef[int, 
                                                ] else coef[int]), res)
  
    coefficients=res
    structure(coefficients, class = "dummy_coef",  matrix = isM)
  
    res2=list(coefficients=coefficients,coef.means=coef.means,
              coef.mode=coef.mode)
    
    structure(res2, class = "dummy_coef.glmb",  matrix = isM)
}


#' @rdname dummy.coef.glmb
#' @method print dummy_coef.glmb
#' @importFrom utils head
#' @export

print.dummy_coef.glmb<-function(x, ...){

cat("\nPrinting coefficients (first 5 records)\n")  
  cn=names(x$coefficients)
  for(i in 1:length(x$coefficients)){
    cat(cn[i],"\n")
    print(head(x$coefficients[[i]],5))
  }
  
  
cat("coef.means:","\n")
print(x$coef.means)
cat("coef.mode:","\n")
print(x$coef.mode)

}
