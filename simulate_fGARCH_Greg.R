##############################
# Functional GARCH(1,1) model
##############################
library('Sim.DiffProc')
fgarch_sim <- function(point_grid, size_sample)
{
    times = 1:point_grid/point_grid
    error_sim <- function(grid_point, samplesize)
    {
        T = 2^(400)/(log(2))
        epsilon = matrix(NA, grid_point, samplesize)
        for(ik in 1:samplesize)
        {
            epsilon[,ik] = BM(x = 0, t0 = 0, T = (2^(400 * 1)/(log(2))), N = (grid_point - 1)) * (1/sqrt( (T * (1:(grid_point))/(grid_point) )))
        }
        return(epsilon)
    }    

    int_approx <- function(x)
    {
        temp_n = NROW(x)
        return((1/(temp_n)) * sum(x))
    }

    no_proc = 1000 + size_sample
    sim_sigma2_matrix = sim_garch_matrix = matrix(NA, point_grid, no_proc)
    error_matrix = error_sim(grid_point = point_grid, samplesize = no_proc)
    
    #fGARCH parameters:

    beta_par = alpha_par = function(t,s)
    {
        return(6 * t * (1-t) * s * (1-s))
    }
    delta_par = 0.01

    #fill in process value for first day and sig2 values for first day:

    sim_sigma2_matrix[,1] = rep(delta_par, point_grid)
    sim_garch_matrix[,1] = delta_par * error_matrix[,1]

    #fill in the rest of sigma2 and functional garch process matrices:
    for(j in 2:no_proc)
    {
        #first fill in sigma2 column:
        for(i in 1:point_grid)
        {
            vector_for_alpha_op = alpha_par(times[i], times) * ((sim_garch_matrix[,(j-1)])^2)
            vector_for_beta_op = beta_par(times[i], times) * ((sim_sigma2_matrix[,j-1]))
            sim_sigma2_matrix[i,j] = delta_par + int_approx(vector_for_alpha_op) + int_approx(vector_for_beta_op)
        }		
        #fill in garch process values for column:
        sim_garch_matrix[,j] = sqrt(sim_sigma2_matrix[,j]) * error_matrix[,j]
    }
    return(list(garch_mat = sim_garch_matrix[,1001:(1000 + size_sample)], 
                sigma_mat = sim_sigma2_matrix[,1001:(1000 + size_sample)]))
}

