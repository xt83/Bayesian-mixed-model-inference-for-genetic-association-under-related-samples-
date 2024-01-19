library(gdata)
library(R.matlab)
library(statmod)
library(Rlab)
library(Matrix)
library(MASS)
library(matrixcalc)
library(lqmm)
library(mvtnorm)

# Data prerpocessing
## Kinship matrix 
prep.kin<-read.table("HCP_pedigree.grm.sp",header = F)
k <- nrow(prep.kin)
# n: number of subjects
n <- length(unique(prep.kin$V1))
kinship.ori <- array(rep(0, n*n), dim = c(n,n))
for (m in 1:k) {
  i <- prep.kin[m,1]+1
  j <- prep.kin[m,2]+1
  kinship.ori[i,j]<-prep.kin[m,3]
  kinship.ori[j,i]<-prep.kin[m,3]
}
# Extract ID for kinship 
ped <- read.csv("pedigree_HPC.csv")
ped <- ped[complete.cases(ped),]
id_ped <- ped$Subject
rownames(kinship.ori) <- id_ped
colnames(kinship.ori) <- id_ped

## A - connectivity matrix
hcp.matl <- readMat("HCP_networkdata.mat")
id_A <- hcp.matl$all.id

## Covariate: X0 -- 1st column: ID
load("X0.rdata")
id_X0 <- X0$id
X0  <- X0[, c(1,3,2,4:ncol(X0))]

## SNP 
load("sigSNP.rdata")
id_SNP <- rownames(myf)

## Overlap id selection (covariate & kinship & SNP & A)
id.temp1 <- Reduce(intersect,list(id_ped, id_X0))
id.temp2 <- Reduce(intersect,list(id.temp1, id_SNP))
id <- Reduce(intersect,list(id.temp2, id_A))
N_subj <- length(id)    ## Number of subjects: 1010

km.temp <- kinship.ori[rownames(kinship.ori) %in% id, ]
kinship <- km.temp[ ,colnames(km.temp) %in% id ] #Kinship matrix generation

A.ori <- hcp.matl$all.network 
m.id_A <- as.data.frame(id_A)
colnames(m.id_A) <- 'id'
m.id_A$ind <- "A"
m.id <- as.data.frame(id)
colnames(m.id) <- 'id'
m.id$ind <- "ID"
test <- merge(x=m.id_A,y=m.id,by="id",all=TRUE)
aleft <- which(is.na(test$ind.y))
A <- A.ori[ , ,-aleft]

sigSNP <- myf[rownames(myf) %in% id, ]

X0 <- X0[X0$id %in% id, ]

## Extract the unique edges
# N_subj: number of subjects
V <- dim(A)[1]
n.up <- V*(V-1)/2
y <- array(rep(NA, N_subj*n.up), dim = c(N_subj,n.up)) 

for (n in 1:N_subj) {
  up_mat <- upperTriangle(A[,,n], diag=FALSE, byrow=TRUE)
  y[n, ] <- up_mat
}

############## Preliminary Step ##############
## 1. Covariates adjustment: age/gender/10 genetic PCs
X <- as.matrix(X0[,-1])
q <- dim(X)[2]
I <- diag(N_subj)
P0 <- I-X%*%solve(t(X)%*%X)%*%t(X)
U <- svd(P0)$u[,1:(N_subj-q)]
Y <- t(U)%*%y
# Assign Y back into upper triangular mat
N_new_subj <- dim(Y)[1]
A_new <- array(rep(0, V*V*N_new_subj), dim = c(V,V,N_new_subj))
for (n in 1:N_new_subj) {
  temp <- Y[n,]
  mat <- A_new[,,n]
  mat[lower.tri(mat, diag=FALSE)] <- temp
  A_new[,,n] <- t(mat)
  A_new[,,n] <- A_new[,,n]+t(A_new[,,n])
}

## 2. select one SNP from 1860 SNP
sigSNP <- myf[rownames(myf) %in% id,]
args <- commandArgs(trailingOnly = TRUE)
id.SNP <- args[1]
X_snp = sigSNP[, id.SNP] 

############## Model fitting ##############
## Algorithm
## update \alpha
update_alpha<-function(h,H,p,sigma_a,kinship,sigma_e,alpha_h,eta,Atrain,xtrain,sigma_prim){
  n=length(xtrain)
  matrix_I <- diag(n)
  v=dim(alpha_h)[2]
  A_tild<-array(rep(0, v*v*n), dim = c(v,v,n))
  sub <- array(rep(0, v*v), dim = c(v,v))
  term_1 <- 0
  temp <- 0
  int<-1:H
  int1<-int[-h]
  int2<-1:v
  int3<-int2[-p]
  mid <- solve(sigma_a*kinship+sigma_e*matrix_I)
  for (j in int1) {
    m1<-eta[j]*outer(alpha_h[j,],alpha_h[j,],"*")
    sub <-sub+m1
  }
  A_tild <- Atrain-outer(sub,xtrain,"*")[,,,1]
  for (g in int3) {
    alpha_hat <- alpha_h[h,g]*eta[h]*xtrain
    A_tild_pg <- A_tild[g,p,]
    term_1 <- term_1 + crossprod(alpha_hat,mid)%*%alpha_hat
    temp <- temp + crossprod(A_tild_pg,mid)%*%alpha_hat
  }
  var_alpha <- 1/(term_1 + sigma_prim[h,p])
  mean_alpha <- temp*var_alpha
  alpha <- rnorm(n=1, mean = mean_alpha, sd = sqrt(var_alpha))
  return(alpha)
}

# update \eta
# temp function for eta updates
form_Q<-function(l,p,h,xtrain,alpha_h){
  n<-length(xtrain)
  Q<-matrix(0,nrow=n,ncol=h)
  for(j in 1:h){
    Q[,j]<-alpha_h[j,l]*alpha_h[j,p]*xtrain
  }
  return(Q)
}


update_eta <- function(h,xtrain,kinship,Atrain,alpha_h,sigma_a,sigma_e,w_vec){
  n <- dim(kinship)[1]
  matrix_I <- diag(n)
  term_1 <- matrix(0,nrow=h,ncol=h)
  temp <- matrix(0,nrow=1,ncol=h)
  D=diag(x=w_vec)
  v=dim(alpha_h)[2]
  mid<-solve(sigma_a*kinship+sigma_e*matrix_I)
  for (p in 2:v) {
    for (l in 1:(p-1)) {
      Q_lp <- form_Q(l,p,h,xtrain,alpha_h)
      A_tild_lp <- Atrain[l,p,]
      term_1 <- term_1 + (crossprod(Q_lp,mid)%*%Q_lp)
      temp <- temp + crossprod(A_tild_lp,mid)%*%Q_lp
    }
  }
  var_eta <- solve(term_1+solve(D))
  mean_eta <- temp%*%var_eta
  eta<-rmvnorm(1, mean = mean_eta, sigma = var_eta)
  return(eta)
}

# update \omega
update_omega_h <- function(h,eta,v1,v2){
  # tau=1
  p_1<- 1/(sqrt(v2))*exp(-0.5*(eta[h]/v2)^2)
  # tau=0
  p_0<- 1/(sqrt(v1))*exp(-0.5*(eta[h]/v1)^2)
  if((p_0==0)||(p_1==0)){
    tau=1
  }
  else{
    prob_1=p_1/(p_1+p_0)
    tau=rbern(1, prob=prob_1)
  }
  if (tau==1){
    omega_h=v2
  } else{
    omega_h=v1
  }
  return(list(omega_h=omega_h,indicator=tau))
}

#update \sigma_prim
update_sigma_prim <- function(v,alpha_hp){
  mean <- v/abs(alpha_hp)
  sigma_prim<-rinvgauss(n = 1,mean = mean,dispersion = v^2)
  return(sigma_prim)
}

# temp function for sigma updates
#' @param opt 1 or 2 -- 1: sigma_a; 2: sigma_e
form_pi <- function(H,eta,opt,sigma_a,sigma_e,alpha_h,alpha=0.1,beta=0.1,kinship,xtrain,Atrain){
  sigma_list<- c(sigma_a,sigma_e)
  n <- dim(kinship)[1]
  matrix_I <- diag(n)
  sigma_sum <- sigma_a*kinship+sigma_e*matrix_I
  mu <- array(rep(0, v*v*n), dim = c(v,v,n)) 
  mu_0 <-array(rep(0, v*v), dim = c(v,v))  
  temp <- 0
  for (h in 1:H) {
    mu_0<- mu_0 + eta[h]*tcrossprod(alpha_h[h,],alpha_h[h,])
  }
  mu <- outer(mu_0,xtrain,"*")[,,,1]
  inv.sigsum<-solve(sigma_sum)
  par1 <- log((det(2*pi*sigma_sum))^(-0.5))
  ele <- rep(NA,v*(v-1)/2)
  ii <-0
  for (p in 2:v) {
    for (l in 1:(p-1)) {
      ii <- ii+1
      sub <- Atrain[l,p,] - mu[l,p,]
      ele[ii] <- par1 + log(exp(-0.5*(crossprod(sub,inv.sigsum)%*%sub)))
    }
  }
  llh <- sum(ele)
  sigma <- sigma_list[opt]
  term2 <- sigma^(-alpha-1)*exp(-beta/sigma)
  pi <- llh+log(term2)
  return(pi)
}

# update sigma_e
update_sigma_e <- function(rho,H,eta,sigma_a,sigma_e,alpha_h,kinship,xtrain,Atrain,d_e){
  sigma_p <- abs(rnorm(1, sigma_e,rho))
  pi_num <- form_pi(H=H,eta=eta,opt=2,sigma_a=sigma_a,sigma_e=sigma_p,alpha_h=alpha_h,
                    alpha=0.1,beta=0.1,kinship=kinship,xtrain=xtrain,Atrain=Atrain)
  pi_den <- d_e
  ind<-pi_num-pi_den
  if (is.nan(ind)){
    new_sigma_e = sample(c(sigma_p, sigma_e), size=1, prob=c(0.5,0.5))
    new_de <- ifelse(new_sigma_e==sigma_e,d_e,pi_num)
  }
  else if ( ind>=0){
    new_sigma_e=sigma_p
    new_de <- pi_num
  }
  else {
    R <- exp(ind)
    new_sigma_e = sample(c(sigma_p, sigma_e), size=1, prob=c(R,1-R))
    new_de <- ifelse(new_sigma_e==sigma_e,d_e,pi_num)
  }
  return(list(sigma_e=new_sigma_e,d_e=new_de))
}

# update sigma_a
update_sigma_a <- function(rho,H,eta,sigma_a,sigma_e,alpha_h,kinship,xtrain,Atrain,d_a){
  sigma_p <- abs(rnorm(1, sigma_a,rho))
  pi_num <- form_pi(H=H,eta=eta,opt=1,sigma_a=sigma_p,sigma_e=sigma_e,alpha_h=alpha_h,
                    alpha=0.1,beta=0.1,kinship=kinship,xtrain=xtrain,Atrain=Atrain)
  pi_den <- d_a
  ind<-pi_num-pi_den
  if (is.nan(ind)){
    new_sigma_a = sample(c(sigma_p, sigma_a), size=1, prob=c(0.5,0.5))
    new_da <- ifelse(new_sigma_a==sigma_e,d_a,pi_num)
  }
  else if ( ind>=0){
    new_sigma_a=sigma_p
    new_da <- pi_num
  }
  else {
    R <- exp(ind)
    new_sigma_a = sample(c(sigma_p, sigma_a), size=1, prob=c(R,1-R))
    new_da <- ifelse(new_sigma_a==sigma_e,d_a,pi_num)
  }
  return(list(sigma_a=new_sigma_a, d_a=new_da))
}

### Initialization
delta <- t(U)%*%kinship%*%U

#' @param kinship: delta = transformed kinship matrix: t(U)*kinship*U
#' @param A: A_new = t(U)*A
#' @param H: 3 (1st try)
#' @param X: t(U)*X_snp
#' @param alpha: H = 3, v = 87
#' @param v: V = 87
#' 

Atrain <- A_new
xtrain <- t(U)%*%X_snp
kinship <- delta
v <- V 
H <- 3
h <- H
alpha_h <- array(rnorm(H*V, 0, 1), dim = c(H,V))
eta <- c(1,1,1)
sigma_e<-0.05
sigma_a<-0.05
v1<-0.5
v2<-2
omega_h <- rep(2,h)
sigma_prim <- array(rep(0.5, h*v), dim = c(h,v))
rho=0.5
vp=0.8

len <- 1000
sigmae_new<-vector("list", len)
sigmaa_new<-vector("list", len)
alpha_new<-vector("list", len)
eta_new<-vector("list", len)
omega_new<-vector("list", len)
sigmaprim_new<-vector("list", len)
tau_new<-vector("list", len)
next.d_e <- form_pi(H=H,eta=eta,opt=2,sigma_a=sigma_a,sigma_e=sigma_e,alpha_h=alpha_h,
                    alpha=0.1,beta=0.1,kinship=kinship,xtrain=xtrain,Atrain=Atrain)
next.d_a <- form_pi(H=H,eta=eta,opt=1,sigma_a=sigma_a,sigma_e=sigma_e,alpha_h=alpha_h,
                    alpha=0.1,beta=0.1,kinship=kinship,xtrain=xtrain,Atrain=Atrain)

#### MCMC implementation ####
tau1<-rep(NA,3)
for(i in 1:1000){
  eta<-update_eta(h,xtrain,kinship,Atrain,alpha_h,sigma_a,sigma_e,omega_h)
  for(j in 1:h){
    upd.ome <-update_omega_h(j,eta,v1,v2)
    omega_h[j]<-upd.ome$omega_h
    tau1[j]<-upd.ome$indicator
  }
  set_e<-update_sigma_e(rho,h,eta,sigma_a,sigma_e,alpha_h,kinship,xtrain,Atrain,next.d_e)
  set_a<-update_sigma_a(rho,h,eta,sigma_a,sigma_e,alpha_h,kinship,xtrain,Atrain,next.d_a)
  sigma_e <- set_e$sigma_e
  sigma_a <- set_a$sigma_a
  next.d_e <- set_e$d_e
  next.d_a <- set_a$d_a
  for(f in 1:h){
    for(p in 1:v){
      alpha_h[f,p]<-update_alpha(f,3,p,sigma_a,kinship,sigma_e,alpha_h,eta,Atrain,xtrain,sigma_prim)
    }
  }
  for(f in 1:h){
    for(p in 1:v){
      sigma_prim[f,p]<-update_sigma_prim(vp,alpha_h[f,p])
    }
  }
  print(i)
  alpha_new[[i]]<-alpha_h
  eta_new[[i]]<-eta
  omega_new[[i]]<-omega_h
  sigmaa_new[[i]]<-sigma_a
  sigmae_new[[i]]<-sigma_e
  tau_new[[i]]<-tau1
  sigmaprim_new[[i]] <- sigma_prim
}
