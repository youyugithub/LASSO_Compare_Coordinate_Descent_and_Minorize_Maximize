# LASSO_Compare_Coordinate_Descent_and_Minorize_Maximize
LASSO_Compare_Coordinate_Descent_and_Minorize_Maximize
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
