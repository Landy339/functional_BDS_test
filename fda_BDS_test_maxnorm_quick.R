################################################################################################################################
# fda_obs: functional data set of dimension u * N where u is the number of discerte points and N is the number of observations.
# m: embedding dimension hyperparameter. Common choices are 2, 3 or 4.
# r: distance threshold value. Common choices of r are 0.75s.d.(fda_obs), s.d.(fda_obs) and 1.25s.d.(fda_obs)
# If p-value is smaller than alpha = 0.05, we reject null hypothesis and conclude rejection of independence.
################################################################################################################################

fda_BDS_test_maxnorm_quick <- function(fda_obs, m, r)
{
    n = ncol(fda_obs)
    indicator = matrix(NA, nrow=n, ncol=n)
    for(i in 1:n)
    {
        for(j in 1:n)
        {
            m_history1 = fda_obs[,i]
            m_history2 = fda_obs[,j]
            indicator[i,j] = (max(colMaxs(abs(m_history1-m_history2), value = TRUE))<r)
        }
    }
    bds_c_T = (sum(rowSums(indicator)^2)-3*sum(indicator)+2*n)/(n*(n-1)*(n-2))

    counts = (n-m+1)*(n-m)/2
    M = indicator[1:(n-m+1),1:(n-m+1)]
    for(m_ind in 2:m)
    {
        M = M*indicator[m_ind:(n+m_ind-m),m_ind:(n+m_ind-m)]
    }
    counts_bds_m = sum(M[upper.tri(M, diag = FALSE)])
    bds_c_m = counts_bds_m/counts

    M = indicator[m:n,m:n]
    counts_bds_m1 = sum(M[upper.tri(M, diag = FALSE)])
    bds_c_m1 = counts_bds_m1/counts

    counts = n*(n-1)/2
    M = indicator[1:n,1:n]
    counts_bds_m1 = sum(M[upper.tri(M, diag = FALSE)])
    C = counts_bds_m1/counts

    A = 0
    for(j in 1:(m-1))
    {
        A = A+(bds_c_T^(m-j))*(C^(2*j))
    }
    sigma = 4*(bds_c_T^m+2*A+((m-1)^2)*(C^(2*m))-(m^2)*bds_c_T*(C^((2*m)-2)))
    BDS_statistic = (bds_c_m-(bds_c_m1^m))*sqrt(n-m+1)/sqrt(sigma)
    p_value = min(c(pnorm(BDS_statistic),1-pnorm(BDS_statistic)))
    return(list(BDS_statistic, p_value))
}
