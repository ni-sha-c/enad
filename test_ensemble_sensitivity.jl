@everywhere include("Lorenz63.jl")
@everywhere include("ensemble_sensitivity.jl")
@everywhere using LinearAlgebra
using PyPlot
ds = [0.,1.0,0.]
function test_growth_tangent()
	u = rand(3)
	tau_len = 40
	tau_arr = range(30.0, stop=2000, length=tau_len)
	tan_sens = zeros(tau_len)
	for (tau_ind, tau) in enumerate(tau_arr)
		tan_sens[tau_ind] = 
		compute_tangent_sensitivity(rand(3), s0, 
				floor(Int,tau), ds)
	end	
	semilogy(tau_arr, abs.(tan_sens))

end

#function test_growth()
	u0 = rand(3)
	u0 = Step(u0,s0,τ_runup)
	n_dim = length(u0)
	du = rand(n_dim)
	du0 = copy(du)
	tau = 2000
	u = zeros(n_dim,tau)
	u[:,1] = copy(u0)
	for n = 2:tau
		u[:,n] = Step(u[:,n-1],s0,1)
    	end
	theta = 0.
	lyap_exp = 0.
	du_np1 = copy(du)
    	for n = 1:tau-1
		println(n)
		du_np1 = tangent_step(du, u[:,n], s0, zeros(3))
		#lyap_exp += log(norm(du_np1)/norm(du))/tau
		du = copy(du_np1)
    	end
	lyap_exp_alt = (1.0/tau)*log(norm(du_np1)/norm(du0))




#end
