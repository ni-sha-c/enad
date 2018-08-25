@everywhere include("Lorenz63.jl")
@everywhere include("ensemble_sensitivity.jl")
tau = 85
T = 5000
N = 100
n_estimators = Int(T/N)

ds = [0.,1.0,0.]
s = copy(s0)
n_tau = 4
tau_arr = range(60, step=10, length=n_tau)

exp_theta_tauN = zeros(n_tau)
var_theta_tauN = zeros(n_tau)

for (tau_ind, tau) in enumerate(tau_arr)
		exp_theta_tauN[tau_ind] = @sync @distributed (+) for i = 1:T
			compute_tangent_sensitivity(rand(3),
										s,tau,ds)/T
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




