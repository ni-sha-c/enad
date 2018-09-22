from pylab import *
from numpy import *
from scipy.special import *
from scipy.stats import linregress
n_plots = 6
c_b = linspace(10.0,50.0,n_plots)
c_var = linspace(1,500,n_plots)
lambda_1 = linspace(1.0, 1.0, n_plots)
gamma_1 = linspace(1.0, 1.0, n_plots)
n_points = 100
slope_lmsevsT = ones(n_points)
T = logspace(1.0,10.0, n_points)
errT = zeros(n_points)
tau = zeros(n_points)
lmse = zeros(n_points)
def mintau(g1,l1,cb,cvar,T):
    c = 2.0*cb*exp(1.0 + g1/l1)*g1*(l1 + g1)*T/cvar/l1
    tau =  -0.5/l1 + 0.5*lambertw(c)/(l1 + g1)
    return tau

def mse(tau,cb,cvar,g1,l1,T):
    err = tau/T*cvar*exp(2.0*l1*tau) + cb*exp(-2.0*g1*tau)
    return err


#if __name__=="__main__":
def plot_lmse_vs_k():
    n_plots = 6
    n_rvalues = 3
    fig, ax = subplots(1,n_rvalues,figsize=(15,6.5))
    fig1, ax1 = subplots(1,n_rvalues,figsize=(15,6.5))
    # k is defined as cb/cvar*gamma1
    # but here I am keeping gamma1 fixed at 1
    k = linspace(0.1,10.,n_plots) 
    c_var = linspace(10.,500.,n_plots)
    c_b = k*c_var
    tau = zeros(n_points)
    tau_grid = zeros(n_points)
    mse_grid = zeros(n_points)
    lmse = zeros(n_points)

    colors = cm.rainbow(linspace(0,1,n_plots))
    rvalues = array([0.1,1.0,10.0])
    gamma_1 = 1.0
    for r_iter,r in enumerate(rvalues):
        ax[r_iter].tick_params(labelsize=14)
        ax[r_iter].set_title("r = %2.1f" % r)
        ax1[r_iter].tick_params(labelsize=14)
        ax1[r_iter].set_title("r = %2.1f" % r)

        for i in range(n_plots):
            for j in range(n_points):
                tau[j] = mintau(gamma_1, gamma_1/r, c_b[i], c_var[i],\
                    T[j])
                lmse[j] = mse(tau[j],c_b[i],c_var[i],gamma_1,\
                    gamma_1/r,T[j])
            ax[r_iter].loglog(T, lmse, label=r"$ k = %2.2f$" % k[i], \
                color=colors[i])
            tau_grid = linspace(tau[-1]/2.0, 3.0*tau[-1]/2.0, n_points)
            for j in range(n_points):
                mse_grid[j] = mse(tau_grid[j],c_b[i],c_var[i],gamma_1,\
                        gamma_1/r,T[-1])
            ax1[r_iter].semilogy(tau_grid,mse_grid,label=r"$ k = %2.2f$" % k[i], \
                color=colors[i])
            ax1[r_iter].semilogy(tau[-1],mse_grid[n_points//2],"o",ms=10,color=colors[i])

    ax[1].legend()
    ax[1].set_xlabel(r"$T$",fontsize=14)
    ax[0].set_ylabel(r"$\tilde{e}_{\rm min}(T)$",fontsize=14,\
            rotation="horizontal",labelpad=30)
    fig.subplots_adjust(wspace=0.3)
    #savefig("paper/eminvsk_python.png",dpi=500)

    ax1[1].legend()
    ax1[1].set_xlabel(r"$\tau$",fontsize=14)
    ax1[0].set_ylabel(r"$\tilde{e}(\tau,T = 10^{10})$",fontsize=14,\
            rotation="horizontal",labelpad=40)
    fig1.subplots_adjust(wspace=0.3)
    savefig("paper/evstau_python.png",dpi=500)
#if __name__ == "__main__":
def plot_lmse_vs_r():
    n_plots = 6
    n_kvalues = 3
    fig, ax = subplots(1,n_kvalues,figsize=(15,6.5))
    # r is defined as gamma_1/lambda_1
    # but here I am keeping gamma_1 fixed at 1
    r = linspace(0.1,10.,n_plots) 
    lambda_1 = linspace(1.,50.,n_plots)
    gamma_1 = r*lambda_1
    c_var = 100
    tau = zeros(n_points)
    colors = cm.rainbow(linspace(0,1,n_plots))
    kvalues = array([0.1,1.0,10.0])
    for k_iter, k in enumerate(kvalues):
        ax[k_iter].tick_params(labelsize=14)
        ax[k_iter].set_title("k = %2.1f" % k)
        for i in range(n_plots):
            c_b = c_var*k/gamma_1[i]
            for j in range(n_points):
                tau[j] = mintau(gamma_1[i], lambda_1[i], c_b , c_var,\
                    T[j])
                lmse[j] = mse(tau[j], c_b, c_var, gamma_1[i],\
                    lambda_1[i], T[j])
            ax[k_iter].loglog(T, lmse, label=r"$ r = %2.2f$" % r[i], \
                color=colors[i])
    ax[1].legend()
    ax[1].set_xlabel(r"$T$",fontsize=14)
    ax[0].set_ylabel(r"$\tilde{e}_{\rm min}(T)$",fontsize=14,\
            rotation="horizontal",labelpad=30)
    fig.subplots_adjust(wspace=0.3)
    savefig("paper/eminvsr_python.png",dpi=500)


#if __name__=="__main__":
def plot_convergence_vs_ratio():
    n_kvalues = 6
    n_rvalues = 100
    n_points = 100
    fig, ax = subplots(figsize=(8,8))
    kvalues = linspace(0.1,10.,n_kvalues) 
    c_var = 50
    tau = zeros(n_points)
    lmse = zeros(n_points)
    colors = cm.rainbow(linspace(0,1,n_kvalues))
    gamma_1 = 1.0
    rvalues = linspace(0.1,50,n_rvalues)
    slope = zeros(n_rvalues)
    T = logspace(0.5, 10, n_points)
    intercept = zeros(n_rvalues)
    for k_iter, k in enumerate(kvalues):
        c_b = k*c_var/gamma_1
        for r_iter, r in enumerate(rvalues):
            for j in range(n_points):
                tau[j] = mintau(gamma_1, gamma_1/r, c_b, c_var,\
                    T[j])
                lmse[j] = mse(tau[j], c_b, c_var, gamma_1,\
                    gamma_1/r, T[j])
            slope[r_iter], intercept[r_iter], _, _, _ = linregress(log(T), log(lmse))
        ax.loglog(rvalues, -slope, label=r"$ k = %2.2f$" % k, \
                color=colors[k_iter])
        
           
    ax.legend()
    ax.set_xlabel(r"$r$",fontsize=14)
    ax.set_ylabel(r"$\beta$",fontsize=14,\
            rotation="horizontal",labelpad=20)
    ax.grid(True,which='both',axis='both')
    ax.tick_params(labelsize=14)
    #xaxis_formatter = LogFormatter(base=10.0,labelOnlyBase=False,\
    #        minor_thresholds=(inf,inf))
    #ax.xaxis.set_minor_formatter(xaxis_formatter)
    savefig("paper/slopeeminvsT_vs_ratiog1l1.png",dpi=500)


if __name__=="__main__":
#def plot_lorenz63_bias():
    exp_tangent = loadtxt("data/hypersonic_new/expected_value_tangent.txt")
    exp_adjoint = loadtxt("data/hypersonic_new/expected_value_adjoint.txt")
    exp_finite_difference = loadtxt("data/hypersonic_new/expected_value_finite_difference.txt")
    exp_tangent_part2 = loadtxt("data/voyager_new/expected_value_tangent.txt")
    exp_adjoint_part2 = loadtxt("data/voyager_new/expected_value_adjoint.txt")
    exp_finite_difference_part2 = loadtxt("data/voyager_new/expected_value_finite_difference.txt")

    exp_tangent = append(exp_tangent,exp_tangent_part2)
    exp_adjoint = append(exp_adjoint,exp_adjoint_part2)
    exp_finite_difference = append(exp_finite_difference,exp_finite_difference_part2)

    n_tau = exp_tangent.shape[0]
    
    tau_beg = 30
    tau_step_size = 10
    tau_arr = array(range(tau_beg, tau_beg + tau_step_size*n_tau, \
            tau_step_size)) 
    dt = 0.005

    sq_bias_tangent = (exp_tangent - 1.0)**2.0
    sq_bias_adjoint = (exp_adjoint - 1.0)**2.0
    sq_bias_finite_difference = (exp_finite_difference - 1.0)**2.0
        
    slope_tangent, intercept_tangent, _, _, _ = linregress(2.0*tau_arr[50:]*dt, \
            log(sq_bias_tangent[50:]))  

    
    slope_adjoint, intercept_adjoint, _, _, _ = linregress(2.0*tau_arr*dt, \
            log(sq_bias_adjoint))  

    slope_finite_difference, intercept_finite_difference, _, _, _ = \
            linregress(2.0*tau_arr*dt, \
            log(sq_bias_finite_difference))  
    
    fig, ax = subplots(1,1,figsize=(15,15))   
    ax.semilogy(tau_arr*dt, sq_bias_tangent,label="tangent")
    ax.semilogy(tau_arr*dt, sq_bias_adjoint,label="adjoint")
    ax.semilogy(tau_arr*dt, sq_bias_finite_difference, label="FD")
    ax.legend()

    ax.set_xlabel(r"$\tau$",fontsize=14)
    ax.set_ylabel(r"${\rm b}^2$",fontsize=14,\
            rotation="horizontal",labelpad=20)
    ax.grid(True,which='both',axis='both')
    ax.tick_params(labelsize=14)
    #xaxis_formatter = LogFormatter(base=10.0,labelOnlyBase=False,\
    #        minor_thresholds=(inf,inf))
    #ax.xaxis.set_minor_formatter(xaxis_formatter)
    savefig("paper/sqbias_vs_tau.png",dpi=500)


#if __name__ == "__main__":
def plot_lorenz63_variance():
    variance_tangent = loadtxt("data/hypersonic_new/variance_tangent.txt")
    variance_adjoint = loadtxt("data/hypersonic_new/variance_adjoint.txt")
    variance_finite_difference = loadtxt("data/hypersonic_new/variance_finite_difference.txt")
    variance_tangent_part2 = loadtxt("data/voyager_new/variance_tangent.txt")
    variance_adjoint_part2 = loadtxt("data/voyager_new/variance_adjoint.txt")
    variance_finite_difference_part2 = loadtxt("data/voyager_new/variance_finite_difference.txt")

    variance_tangent = append(variance_tangent,variance_tangent_part2)

    variance_adjoint = append(variance_adjoint,variance_adjoint_part2)

    variance_finite_difference = append(variance_finite_difference,variance_finite_difference_part2)

    n_tau = variance_tangent.shape[0]
    
    tau_beg = 30
    tau_step_size = 10
    tau_arr = array(range(tau_beg, tau_beg + tau_step_size*n_tau, \
            tau_step_size)) 
    dt = 0.005

    T = 1000000

    sane_indices_tangent = variance_tangent > 0.0
    sane_indices_adjoint = variance_adjoint > 0.0
    sane_indices_finite_difference = variance_finite_difference > 0.0
    variance_tangent = variance_tangent[sane_indices_tangent]
    variance_adjoint = variance_adjoint[sane_indices_adjoint]
    variance_finite_difference = variance_finite_difference[sane_indices_finite_difference]
    
    
    tau_arr_tangent = copy(tau_arr[sane_indices_tangent])
    tau_arr_adjoint = copy(tau_arr[sane_indices_adjoint])
    tau_arr_finite_difference = copy(tau_arr[sane_indices_finite_difference])
    
    N_times_variance_tangent = (T/tau_arr_tangent)*variance_tangent
    N_times_variance_adjoint = (T/tau_arr_adjoint)*variance_adjoint
    N_times_variance_finite_difference = (T/tau_arr_finite_difference)*\
            variance_finite_difference

    
    slope_tangent, intercept_tangent, _, _, _ = linregress(2.0*tau_arr_tangent[50:]*dt, \
            log(N_times_variance_tangent[50:]))  

    
    slope_adjoint, intercept_adjoint, _, _, _ = linregress(2.0*tau_arr_adjoint*dt, \
            log(N_times_variance_adjoint))  

    slope_finite_difference, intercept_finite_difference, _, _, _ = \
            linregress(2.0*tau_arr_finite_difference*dt, \
            log(N_times_variance_finite_difference))  
    
    fig, ax = subplots(1,1,figsize=(15,15))   
    ax.semilogy(tau_arr_tangent*dt, N_times_variance_tangent,label="tangent")
    ax.semilogy(tau_arr_adjoint*dt, N_times_variance_adjoint,label="adjoint")
    ax.semilogy(tau_arr_finite_difference*dt, N_times_variance_finite_difference, label="FD")
    ax.legend()

    ax.set_xlabel(r"$\tau$",fontsize=14)
    ax.set_ylabel(r"$(T/\tau){\rm var}$",fontsize=14,\
            rotation="horizontal",labelpad=20)
    ax.grid(True,which='both',axis='both')
    ax.tick_params(labelsize=14)
    #xaxis_formatter = LogFormatter(base=10.0,labelOnlyBase=False,\
    #        minor_thresholds=(inf,inf))
    #ax.xaxis.set_minor_formatter(xaxis_formatter)
    savefig("paper/Ntimesvariance_vs_tau.png",dpi=500)
    



