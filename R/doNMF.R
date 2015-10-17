.updateH <-
  function(W,V,H) {
    Hn <- H * ((t(W) %*% V) / ((t(W) %*% W %*% H)))
    Hn[Hn<1e-16] <- 1e-16
    return(Hn)
  }
.updateW <-
  function(W,V,H) {
    Wn <- W * ((V %*% t(H)) / ((W %*% H %*% t(H))))
    Wn[Wn<1e-16] <- 1e-16
    return(Wn)
  }
### Do a non-negative matrix factorization of V
### using rank k
nnmf <- function(V,k=2,niter=2000,verbose=FALSE) {
  m <- ncol(V)
  n <- nrow(V)
  W <- matrix(abs(rnorm(n*k)),ncol=k)
  H <- matrix(abs(rnorm(k*m)),ncol=m)
  h <- H
  w <- W
  for(i in 1:niter) {
    w1 <- .updateW(w,V,h)
    h1 <- .updateH(w,V,h)
    w <- w1
    h <- h1
  }
  return(list(h=h,w=w))
}

.preprocessnnmf <- function(mat,multfactor=10000) {
  return(return(apply(mat,2,function(x) 10000*((rank(x)-1)/(nrow(mat)-1)))))
}

setGeneric("preprocessnnmf",function(object,multfactor=10000,...) {
  standardGeneric("preprocessnnmf")
})

setMethod("preprocessnnmf",c("matrix"),
          function(object,multfactor=10000,...) {
            return(.preprocessnnmf(mat))
          })

setMethod("preprocessnnmf",c("ExpressionSet"),
          function(object,multfactor=10000,...) {
            exprs(object) <- .preprocessnnmf(exprs(object))
          })
