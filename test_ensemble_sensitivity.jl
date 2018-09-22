@everywhere include("Lorenz63.jl")
@everywhere include("ensemble_sensitivity.jl")
@everywhere using LinearAlgebra
@everywhere using SharedArrays
ds = [0.,1.0,0.]
function fitaline(x, y)
	n = length(y)
	@assert length(x)==n
	xbar = sum(x)/n
	ybar = sum(y)/n
	sumxy = dot(x,y)
	sumxsq = sum(x.^2)
	varx = sumxsq - n*xbar*xbar
	slope = (sumxy - n*xbar*ybar)/varx
	intercept = (ybar*sumxsq - xbar*sumxy)/varx
	return slope, intercept

end
function test_growth_tangent()
	u = rand(3)
	tau_len = 1000
	tau_arr = range(30.0, stop=20000, length=tau_len)
	tan_sens = zeros(tau_len)
    adj_sens = zeros(tau_len)
	for (tau_ind, tau) in enumerate(tau_arr)
		tan_sens[tau_ind] = 
		compute_tangent_sensitivity(rand(3), s0, 
				floor(Int,tau), ds)
		adj_sens[tau_ind] = 
		compute_adjoint_sensitivity(rand(3), s0, 
				floor(Int,tau), ds)

	end	
    slope_tangent, intercept_tangent = fitaline(tau_arr[20:end], log.(abs.(tan_sens[20:end])))
    slope_adjoint, intercept_adjoint = fitaline(tau_arr[20:end], log.(abs.(adj_sens[20:end])))
    @assert abs(slope_tangent - slope_adjoint) < 1.e-1
    @assert abs(slope_tangent/dt - 0.96) < 0.05
	
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
		dw_nm1 = adjoint_step(dw, u[:,n-1], s0)
		lyap_exp_adjoint += log(norm(dw_nm1)/norm(dw))/(tau*dt)
		dw[:] = dw_nm1
    	end
	println(lyap_exp_adjoint)
	@assert abs(lyap_exp_adjoint - 0.96) < 0.1  

end

function test_expectation_at_tau(tau::Int64, N_arr::Array{Int64,1}, method::String="tangent")
		if(method=="tangent")
				ens_sens_method = compute_tangent_sensitivity
		elseif(method=="adjoint")
				ens_sens_method = compute_adjoint_sensitivity
		elseif(method=="fd")
				ens_sens_method = compute_finite_difference_sensitivity
		end
		ds0 = [0.,1.,0.]
		n_N = length(N_arr)
		exp_theta_tauN = zeros(n_N)
		for (i, N) in enumerate(N_arr)
				exp_theta_tauN[i] = @sync @distributed (+) for n = 1:N	
						ens_sens_method(rand(3), s0, tau, ds0)/N
				end
		end
		cumN_arr = cumsum(N_arr)
		cumexp_theta_tauN = cumsum(N_arr.*exp_theta_tauN)./cumN_arr
		return cumexp_theta_tauN			
end
function test_variance_at_tau(tau::Int64, N::Int64, method::String="tangent")
		if(method=="tangent")
				ens_sens_method = compute_tangent_sensitivity
		elseif(method=="adjoint")
				ens_sens_method = compute_adjoint_sensitivity
		elseif(method=="fd")
				ens_sens_method = compute_finite_difference_sensitivity
		end
		ds0 = [0.,1.,0.]
		n_N = 10000
		theta_tauN = zeros(n_N)
		for i in 1:n_N
			theta_tauN[i] = @sync @distributed (+) for n = 1:N	
				ens_sens_method(rand(3), s0, tau, ds0)/N
			end
		end
		return theta_tauN		
end
function test_rare_event(tau::Int64, method::String="tangent")
	if(method=="tangent")
				ens_sens_method = compute_tangent_sensitivity
		elseif(method=="adjoint")
				ens_sens_method = compute_adjoint_sensitivity
		end
		ds0 = [0.,1.,0.]
		N = 100000
		exp_theta_tauN = 0.0
		u = rand(3)
		exp_theta_tauN = @sync @distributed (+) for n = 1:N	
			u = Step(u,s0,2000)
			ens_sens_method(u, s0, tau, ds0)/N
		end
		print(exp_theta_tauN)
		n_samples = 1000
		u0 = SharedArray(rand(3,n_samples))
		theta_tauN = SharedArray(zeros(n_samples))
		@sync @distributed for i = 1:n_samples
				u0[:,i] = Step(u0[:,i], s0, 2000)
		end

		@sync @distributed for i = 1:n_samples	
				theta_tauN[i] = ens_sens_method(u0[:,i], s0, tau, ds0)	
		end

		factor_exp_thetaN = abs.(theta_tauN./exp_theta_tauN)
	
		return u0, factor_exp_thetaN
end
		





		
		


	
		


