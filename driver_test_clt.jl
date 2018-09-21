@everywhere include("test_ensemble_sensitivity.jl")
using DelimitedFiles
n_tau = 5
tau_arr = [800, 1000, 1200, 1800, 3000]
n = 10000
N = 200000
theta_tauN_tangent = zeros(n,n_tau) 
theta_tauN_adjoint = zeros(n,n_tau)
theta_tauN_fd = zeros(n,n_tau)

for (i, tau) in enumerate(tau_arr)
	print(i)
	theta_tauN_tangent[:,i] = test_variance_at_tau(tau, N, "tangent")
	theta_tauN_adjoint[:,i] = test_variance_at_tau(tau, N, "adjoint")
	theta_tauN_fd[:,i] = test_variance_at_tau(tau, N, "fd")
end

open("data/test_variance_estimator_tangent.txt","w") do io
		writedlm(io, theta_tauN_tangent)
end
open("data/test_variance_estimator_adjoint.txt","w") do io
		writedlm(io, theta_tauN_adjoint)
end
open("data/test_variance_estimator_fd.txt","w") do io
		writedlm(io, theta_tauN_fd)
end


