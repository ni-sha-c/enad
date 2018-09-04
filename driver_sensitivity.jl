@everywhere include("Lorenz63.jl")
@everywhere include("ensemble_sensitivity.jl")
using DelimitedFiles
TT = 5000000
n_estimators = 10000
T = 1000000
ds = [0.,1.0,0.]
s = copy(s0)
n_tau = 50 
tau_arr = range(30, step=10, length=n_tau)

exp_theta_tauN = zeros(n_tau)
var_theta_tauN = zeros(n_tau)

for (tau_ind, tau) in enumerate(tau_arr)
	N = floor(Int,(T/tau))
	exp_theta_tauN[tau_ind] = @sync @distributed (+) for i = 1:TT
			compute_tangent_sensitivity(rand(3),
										s,tau,ds)/TT
	end
		
	var_theta_tauN[tau_ind] = @sync @distributed (+) for i = 1:n_estimators
		theta = 0.0
		for j = 1:N
			theta += compute_tangent_sensitivity(rand(3),
														s, tau, ds)/N
		end
		theta*theta/n_estimators
	end
		
	var_theta_tauN[tau_ind] -= exp_theta_tauN[tau_ind]*exp_theta_tauN[tau_ind]
		
end

open("expected_value_tangent.txt","w") do io
		writedlm(io, exp_theta_tauN)
end
open("variance_tangent.txt","w") do io
		writedlm(io, var_theta_tauN)
end


