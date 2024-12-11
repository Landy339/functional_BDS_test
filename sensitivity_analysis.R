accuracy_L1norm=matrix(,nrow=1,ncol=5)
accuracy_L2norm=matrix(,nrow=1,ncol=5)
accuracy_Maxnorm=matrix(,nrow=1,ncol=5)
N=c(100,250,500,750,1000)
for(K in 1:5){
    random_bds_statistic_L1norm=matrix(,nrow=200,ncol=1)
    random_bds_statistic_L2norm=matrix(,nrow=200,ncol=1)
    random_bds_statistic_Maxnorm=matrix(,nrow=200,ncol=1)
    for(i in 1:200){
      Sim=fgarch_sim(100,N[K])$garch_mat
      X=Sim
      S1=fda_BDS_test_L1norm_quick(X,3,R_L1norm[4])
      random_bds_statistic_L1norm[i]=S1[[1]]
      S2=fda_BDS_test_L2norm_quick(X,3,R_L2norm[4])
      random_bds_statistic_L2norm[i]=S2[[1]]
      S3=fda_BDS_test_maxnorm_quick(X,3,R_MaxNorm[4])
      random_bds_statistic_maxnorm[i]=S3[[1]]
      print(i)
    }
    accuracy_L1norm[1,K]=sum(abs(random_bds_statistic_L1norm)>1.96)/200
    accuracy_L2norm[1,K]=sum(abs(random_bds_statistic_L2norm)>1.96)/200
    accuracy_maxnorm[1,K]=sum(abs(random_bds_statistic_maxnorm)>1.96)/200
  }

M=seq(2,10)
accuracy_L2norm=matrix(,nrow=8,ncol=9)
for(m in 1:9){
  for(r in 1:8){    
    print(paste('m=',M[m]))
    print(paste('r=',r))
    random_bds_statistic_L2norm=matrix(,nrow=200,ncol=1)
    for(i in 1:200){
      Sim=fgarch_sim(100,500)$garch_mat
      X=Sim
      S=fda_BDS_test_L2norm_quick(X,M[m],R_L2norm[r])
      random_bds_statistic_L2norm[i]=S[[1]]
      print(i)
    }
    accuracy_L2norm[r,m]=sum(abs(random_bds_statistic_L2norm)>1.96)/200
  }
}

X_mean=rowMeans(X)
L2_distance=matrix(,nrow=500,ncol=1)
for(i in 1:ncol(X)){
  L2_distance[i]=sqrt(sum((X[,i]-X_mean)^2))
}
sd=sqrt(mean(L2_distance^2))
R_L2norm=seq(0.25,2,by=0.25)*sd

Hbase=c(1,10)
pbase=c(3,4,5,8,10)
power_GK=matrix(,nrow=2,ncol=5)
for(m in 1:2){
  for(r in 1:5){    
    print(paste('H=',Hbase[m]))
    print(paste('p=',pbase[r]))
    p_GK_trial=matrix(,nrow=200,ncol=1)
    for(i in 1:200){
      X=fgarch(days_to_simulate = 500,no_grid=100, seed_number = i)
      S=GK_test(X,K=pbase[r],H=Hbase[m])
      p_GK_trial[i]=S[[1]]
      print(i)
    }
    power_GK[m,r]=sum(p_GK_trial<0.05)/200
  }
}

#15 * t * (1 - t) * s * (1 - s)

for(m in 2){
  for(r in 1:5){    
    print(paste('H=',Hbase[m]))
    print(paste('p=',pbase[r]))
    p_GK_trial=matrix(,nrow=200,ncol=1)
    for(i in 1:200){
      X=fgarch_sim(100,500)$garch_mat
      S=GK_test(X,K=pbase[r],H=Hbase[m])
      p_GK_trial[i]=S[[1]]
      print(i)
    }
    power_GK[m,r]=sum(p_GK_trial<0.05)/200
  }
}
