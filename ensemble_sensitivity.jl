τ_runup = 2000 
function compute_tangent_sensitivity(u0::Array{Float64,1},s::Array{Float64,1},
							 tau::Int64,ds::Array{Float64,1})
    u0 = Step(u0,s,τ_runup)
	n_dim = length(u0)
	du = zeros(n_dim)
	u = zeros(n_dim,tau)
	u[:,1] = copy(u0)
	J_weights = dt*ones(tau)
    	J_weights[1] /= 2.0
	J_weights[tau] /= 2.0

	for n = 2:tau
		u[:,n] = Step(u[:,n-1],s,1)
    	end
	theta = 0.
    	for n = 1:tau-1
		du = tangent_step(du, u[:,n], s, ds)
		theta += du'*nabla_objective(u[:,n+1],s)*
						J_weights[n+1]
    	end
	return theta
end
function compute_finite_difference_sensitivity(u0::Array{Float64,1},s::Array{Float64,1},
							 tau::Int64,ds::Array{Float64,1})
    u0 = Step(u0,s,τ_runup)
	n_dim = length(u0)

	uplus = zeros(n_dim,tau)
	uplus[:,1] = copy(u0)
	uminus = zeros(n_dim,tau)
	uminus[:,1] = copy(u0)

	J_weights = dt*ones(tau)
	J_weights[1] /= 2.0
	J_weights[tau] /= 2.0
	
	splus = copy(s)
	sminus = copy(s)
	epsi = 5.e-3
	splus[2] += epsi
	sminus[2] -= epsi

    theta = 0.0
	for n = 2:tau
		uplus[:,n] = Step(uplus[:,n-1],splus,1)
		uminus[:,n] = Step(uminus[:,n-1],sminus,1)
		theta += J_weights[n]*(objective(uplus[:,n-1],splus) -
							   objective(uminus[:,n-1],sminus))/
								(2*epsi)
    end
	

	return theta
end
function compute_adjoint_sensitivity(u0::Array{Float64,1},s::Array{Float64,1},
							 tau::Int64,ds::Array{Float64,1})
    u0 = Step(u0,s,τ_runup)
	n_dim = length(u0)
	du = zeros(n_dim)
	u = zeros(n_dim,tau)
	u[:,1] = copy(u0)
	J_weights = dt*ones(tau)
    J_weights[1] /= 2.0
	J_weights[tau] /= 2.0

	for n = 2:tau
		u[:,n] = Step(u[:,n-1],s,1)
    end
    
	du = adjoint_source(du, u[:,tau], s0, 
						J_weights[tau])
	theta = 0.0
    for n = tau-1:-1:1
		theta += du'*(dt*dfds(u[:,n],s0,ds))
		du = adjoint_step(du, u[:,n], s)
		du = adjoint_source(du, u[:,n], s0, 
							J_weights[n])	
    end
	return theta
end

