# LASSO_Compare_Coordinate_Descent_and_Minorize_Maximize

Gaussian:

```
set.seed(0)
nobs<-1000
ndim<-20
X<-matrix(runif(nobs*ndim),nobs,ndim)
beta_true<-c(1,0.8,0.6,0.4,0.2,rep(0,15))
Y<-X%*%beta_true+rnorm(nobs)*0.5
lambda<-0.1

## Coordinate Descend

beta_est<-solve(t(X)%*%X)%*%t(X)%*%Y
for(iter in 1:1000){
  beta_old<-beta_est
  for(pp in 1:ndim){
    temp_beta<-beta_est
    temp_beta[pp]<-0
    R<-Y-X%*%temp_beta
    ARG1<-X[,pp]%*%R
    ARG2<-lambda*nobs
    if(ARG1>ARG2){
      S<-ARG1-ARG2
    }else if(ARG1< -ARG2){
      S<-ARG1+ARG2
    }else{
      S<-0
    }
    beta_est[pp]<-S/sum(X[,pp]^2)
  }
  eps<-sum(abs(beta_est-beta_old))
  if(eps<1e-10){
    print(iter)
    break
  }
}

(beta_CD<-c(beta_est))

## Majorization

epsilon<-1e-10
beta_est<-solve(t(X)%*%X)%*%t(X)%*%Y
for(iter in 1:10000){
  beta_old<-beta_est
  dl<- -t(X)%*%(Y-X%*%beta_est)/nobs
  ddl<- t(X)%*%X/nobs
  dpen<- lambda*beta_est/(epsilon+abs(beta_est))
  ddpen<- lambda/(epsilon+abs(beta_est))
  dl<-dl+dpen
  diag(ddl)<-diag(ddl)+ddpen
  beta_est<-beta_est-solve(ddl,dl)
  eps<-sum(abs(beta_est-beta_old))
  if(eps<1e-10){
    print(iter)
    break
  }
}

(beta_MM<-round(beta_est,9))

cbind(beta_MM,beta_CD)

sprintf("%.10f",t(Y-X%*%beta_MM)%*%(Y-X%*%beta_MM)/nobs+lambda*sum(abs(beta_MM)))
sprintf("%.10f",t(Y-X%*%beta_CD)%*%(Y-X%*%beta_CD)/nobs+lambda*sum(abs(beta_CD)))
```

Binomial (logistic)

```
################################
# weighted logistic regression #
################################

set.seed(100)
nobs<-1000
ndim<-20
X<-matrix(rnorm(nobs*ndim),nobs,ndim)
beta_true<-c(1,0.8,0.6,0.4,0.2,rep(0,15))
eta<-X%*%beta_true
ppi<-exp(eta)/(1+exp(eta))
Y<-rbinom(nobs,1,ppi)

## No penalty

## Coordinate Descend

beta_est<-rep(1,ndim)
for(iter in 1:1000){
  beta_old<-beta_est
  eta<-X%*%beta_est
  ppi<-exp(eta)/(1+exp(eta))
  wwi<-ppi*(1-ppi)
  for(pp in 1:ndim){
    beta_est[pp]<-sum(X[,pp]*(X[,pp]*beta_est[pp]+(Y-ppi)/ppi/(1-ppi))*ppi*(1-ppi))/sum(X[,pp]^2*ppi*(1-ppi))
  }
  if(sum(abs(beta_est-beta_old))<1e-10)break
}

cbind(beta_est,coef(glm(Y~X-1,family=binomial(),control=list(epsilon=1e-10))))

## Coordinate Descend

lambda<-0.01
beta_est<-rep(0,ndim)
for(iter in 1:1000){
  beta_old<-beta_est
  eta<-X%*%beta_est
  ppi<-exp(eta)/(1+exp(eta))
  wwi<-ppi*(1-ppi)
  for(pp in 1:ndim){
    ARG1<-sum(X[,pp]*(X[,pp]*beta_est[pp]+(Y-ppi)/ppi/(1-ppi))*ppi*(1-ppi))
    ARG2<-lambda*nobs
    if(ARG1>ARG2){
      S<-ARG1-ARG2
    }else if(ARG1< -ARG2){
      S<-ARG1+ARG2
    }else{
      S<-0
    }
    beta_est[pp]<-S/sum(X[,pp]^2*ppi*(1-ppi))
  }
  print(beta_est)
  if(sum(abs(beta_est-beta_old))<1e-10)break
}
cbind(beta_est,coef(glm(Y~X-1,family=binomial())))
```

Survival (Cox)

```
##########################
# Cox Coordinate Descent #
##########################

set.seed(100)
nobs<-1000
ndim<-20
X<-matrix(rnorm(nobs*ndim),nobs,ndim)
beta_true<-c(1,0.8,0.6,0.4,0.2,rep(0,15))
eta<-X%*%beta_true
DDi<-rexp(nobs,exp(eta))/2
CCi<-rep(1,nobs)
TTi<-pmin(DDi,CCi)
deltai<-TTi<CCi

####################
# CD without LASSO #
####################

beta_est<-rep(0,ndim)
for(iter in 1:1000){
  beta_old<-beta_est
  wwi<-exp(X%*%beta_est)
  for(pp in 1:ndim){
    A1<-0
    A2<-0
    B<-0
    for(ii in 1:nobs){
      if(!deltai[ii])next
      denominator<-0
      numerator_1<-0
      numerator_2<-0
      for(kk in 1:nobs){
        if(TTi[kk]<=TTi[ii])next
        denominator<-denominator+wwi[kk]
        numerator_1<-numerator_1+wwi[kk]*X[kk,pp]
        numerator_2<-numerator_2+wwi[kk]*X[kk,pp]*X[kk,pp]
      }
      B<-B+(X[ii,pp]-numerator_1/denominator)
      A1<-A1+numerator_2/denominator
      A2<-A2+(numerator_1/denominator)^2
    }
    beta_est[pp]<-beta_est[pp]+B/(A1-A2)
  }
  eps<-sum(abs(beta_est-beta_old))
  print(eps)
  if(eps<1e-15){
    print(iter)
    break
  }
}

beta_CD<-beta_est

####################
# NR without LASSO #
####################

beta_est<-rep(0,ndim)
for(iter in 1:1000){
  beta_old<-beta_est
  wwi<-exp(X%*%beta_est)
  
  U<-rep(0,ndim)
  I<-matrix(0,ndim,ndim)
  for(ii in 1:nobs){
    if(!deltai[ii])next
    denominator<-0
    numerator_vec<-rep(0,ndim)
    numerator_mat<-matrix(0,ndim,ndim)
    for(kk in 1:nobs){
      if(TTi[kk]<=TTi[ii])next
      denominator<-denominator+wwi[kk]
      numerator_vec<-numerator_vec+wwi[kk]*X[kk,]
      numerator_mat<-numerator_mat+outer(wwi[kk]*X[kk,],X[kk,])
    }
    U<-U+(X[ii,]-numerator_vec/denominator)
    I<-I+numerator_mat/denominator-outer(numerator_vec/denominator,numerator_vec/denominator)
  }
  beta_est<-beta_est+solve(I,U)
  eps<-sum(abs(beta_est-beta_old))
  print(eps)
  if(eps<1e-15){
    print(iter)
    break
  }
}
beta_NR<-beta_est

library(survival)
cbind(beta_CD,beta_NR,coef(coxph(Surv(TTi,deltai)~X,control=coxph.control(eps=1e-15,iter.max=50,toler.chol=1e-16))))
```
