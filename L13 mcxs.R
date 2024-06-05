
############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 14: SVARs Bayesian Estimation I
############################################################

rm(list=ls())

orthogonal.complement.matrix.TW = function(x){
  # x is a mxn matrix and m>n
  # the function returns a mx(m-n) matrix, out, that is an orthogonal complement of x, i.e.:
  # t(x)%*%out = 0 and det(cbind(x,out))!=0
  N     = dim(x)
  tmp   = qr.Q(qr(x, tol = 1e-10),complete=TRUE)
  out   = as.matrix(tmp[,(N[2]+1):N[1]])
  return(out)
}

# Estimating models with exclusion restrictions: example
# Use Algorithm 1
############################################################

R1            = matrix(0,2,6)
R1[1,1]       = 1
R1[2,4]       = 1
R2            = matrix(c(0,0,0,1,0,0),nrow=1)
A             = t(matrix(c(.5,.5,0,-1.25,.25,0,-1,0,.5),3,3))
Sigma         = matrix(c(1,.5,1,.5,4.25,2.5,1,2.5,3),3,3)
B0.tilde      = t(solve(chol(Sigma)))
B1.tilde      = B0.tilde%*%A
IR.0.tilde    = solve(B0.tilde)
IR.infinity.tilde = solve(B0.tilde-B1.tilde)

fAA           = rbind(IR.0.tilde,IR.infinity.tilde)
p1            = orthogonal.complement.matrix.TW(t(R1 %*% fAA))
p2            = orthogonal.complement.matrix.TW(cbind(t(R2 %*% fAA),p1))
p3            = orthogonal.complement.matrix.TW(cbind(p1,p2))
Q             = t(cbind(p1,p2,p3))
Q[2,]         = -Q[2,]
B0            = Q%*%B0.tilde          # the resulting contemporaneous relationship matrix
B1            = Q%*%B1.tilde          # the resulting autoregressive matrix
IR.0          = solve(B0)             # this matrix has zero in the [1.1] element
IR.infinity   = solve(B0-B1)          # this matrix has zeros in the [1.1] and [1.2] elements 

R1%*%rbind(IR.0,IR.infinity)%*%diag(3)[,1]    # Verifying that the restrictions hold: OK!
R2%*%rbind(IR.0,IR.infinity)%*%diag(3)[,2]    # Verifying that the restrictions hold: OK! (this is a numerical zero)



# Estimating models with sign restrictions: example
# Use Algorithm 2
############################################################
rm(list=ls())
set.seed(123456)
sign.restrictions = c(-1,-1,1,-1,-1,1)
R1            = diag(sign.restrictions)
A             = t(matrix(c(.5,.5,0,-1.25,.25,0,-1,0,.5),3,3))
Sigma         = matrix(c(1,.5,1,.5,4.25,2.5,1,2.5,3),3,3)
B0.tilde      = t(solve(chol(Sigma)))
B1.tilde      = B0.tilde%*%A
IR.0.tilde    = solve(B0.tilde)
IR.1.tilde    = solve(B0.tilde)%*%B1.tilde%*%solve(B0.tilde)

sign.restrictions.do.not.hold = TRUE
i=1
while (sign.restrictions.do.not.hold){
  X           = matrix(rnorm(9),3,3)
  QR          = qr(X, tol = 1e-10)
  Q           = qr.Q(QR,complete=TRUE)
  R           = qr.R(QR,complete=TRUE)
  Q           = t(Q %*% diag(sign(diag(R))))
  B0          = Q%*%B0.tilde
  B1          = Q%*%B1.tilde
  B0.inv      = solve(B0)
  check       = prod(R1 %*% rbind(B0.inv,B0.inv%*%B1%*%B0.inv) %*% diag(3)[,1] > 0)
  if (check==1){sign.restrictions.do.not.hold=FALSE}
  i=i+1
}
i
Q
B0
B1
IR.0        = B0.inv
IR.1        = B0.inv%*%B1%*%B0.inv
IR.0
IR.1
R1 %*% rbind(IR.0,IR.1) %*% diag(3)[,1]
