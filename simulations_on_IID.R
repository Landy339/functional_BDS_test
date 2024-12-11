library(far)
library('Rfast')
#install.packages('svMisc')
require(svMisc)
#max norm BDS test random normality check
M=seq(2,10)
KS_test_table_maxnorm=matrix(,nrow=8,ncol=9)
for(m in 1:9){
  for(r in 6:8){
  print(paste('m=',M[m]))
  print(paste('r=',r))
  random_bds_statistic_maxnorm=matrix(,nrow=200,ncol=1)
  for(i in 1:200){
  Sim=simul.far(m=100,
                n=500,
                base=base.simul.far(24, 5),
                d.rho=diag(c(0.001, 0.01, 0.001, 0.001)),
                alpha=diag(c(0.01, 0.02, 0.001)),
                cst1=0.05)
  X=Sim[["var"]]
  S=fda_BDS_test_maxnorm_quick(X,M[m],R_MaxNorm[r])
  random_bds_statistic_maxnorm[i]=S[[1]]
  print(i)
}
KS_test_table_maxnorm[r,m]=ks.test(random_bds_statistic_maxnorm,'pnorm')$p.value
  }
}

X_mean=rowMeans(X)
MaxNorm_distance=matrix(,nrow=500,ncol=1)
for(i in 1:ncol(X)){
  MaxNorm_distance[i]=max(abs(X[,i]-X_mean))
}
sd=sqrt(mean(MaxNorm_distance^2))
R_MaxNorm=seq(0.25,2,by=0.25)*sd


#max norm BDS test power check
M=seq(2,10)
accuracy_maxnorm=matrix(,nrow=5,ncol=9)
for(m in 1:9){
  for(r in 1:5){
    print(paste('m=',M[m]))
    print(paste('r=',r))
    random_bds_statistic_maxnorm=matrix(,nrow=200,ncol=1)
    for(i in 1:200){
      Sim=simul.far(m=100,
                    n=500,
                    base=base.simul.far(100, 30),
                    d.rho=diag(c(0.1, 0.9, 0.34, 0.2)),
                    alpha=diag(c(0.5, 0.23, 0.018)),
                    cst1=0.2)
      X=Sim[["var"]]
      S=fda_BDS_test_maxnorm_quick(X,M[m],R_MaxNorm[r])
      random_bds_statistic_maxnorm[i]=S[[1]]
      print(i)
    }
    accuracy_maxnorm[r,m]=sum(abs(random_bds_statistic_maxnorm)>1.96)/200
  }
}

