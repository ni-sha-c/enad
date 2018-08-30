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

function test_lyapunov_exponent_tangent()
	u0 = rand(3)
	u0 = Step(u0,s0,τ_runup)
	n_dim = length(u0)
	du = rand(n_dim)
	du0 = copy(du)
	tau = 20000
	u = zeros(n_dim,tau)
	u[:,1] = copy(u0)
	for n = 2:tau
		u[:,n] = Step(u[:,n-1],s0,1)
    	end
	theta = 0.
	lyap_exp_tangent = 0.
	du_np1 = copy(du)
	
    	for n = 1:tau-1
		du_np1 = tangent_step(du, u[:,n], s0, zeros(3))
		lyap_exp_tangent += log(norm(du_np1)/norm(du))/(tau*dt)
		du[:] = du_np1
    	end

	println(lyap_exp_tangent)
	@assert abs(lyap_exp_tangent - 0.96) < 0.1  

end

function test_lyapunov_exponent_adjoint()
	u0 = rand(3)
	u0 = Step(u0,s0,τ_runup)
	n_dim = length(u0)
	dw = rand(n_dim)
	dw0 = copy(dw)
	tau = 20000
	u = zeros(n_dim,tau)
	u[:,1] = copy(u0)
	for n = 2:tau
		u[:,n] = Step(u[:,n-1],s0,1)
    	end
	theta = 0.
	lyap_exp_adjoint = 0.
	dw_nm1 = copy(dw)
	
    	for n = tau:-1:2
		global lyap_exp_adjoint
		global dw_np1, dw
		dw_nm1 = adjoint_step(dw, u[:,n-1], s0)
		lyap_exp_adjoint += log(norm(dw_nm1)/norm(dw))/(tau*dt)
		dw[:] = dw_nm1
    	end
	println(lyap_exp_adjoint)
	@assert abs(lyap_exp_adjoint - 0.96) < 0.1  

end
