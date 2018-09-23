@everywhere using LinearAlgebra
@everywhere include("Lorenz63.jl")
@everywhere using SharedArrays
function compute_CLVs(u0::Array{Float64,1})
	n_dim = length(u0)
	du = rand(n_dim, n_dim)
	du /= norm(du)
	Q_init, R_init = qr(du)
	tau = 10000
	u = zeros(n_dim,tau)
	u[:,1] = copy(u0)
	for n = 2:tau
		u[:,n] = Step(u[:,n-1],s0,1)
    end
	Q = SharedArray(zeros(n_dim, n_dim, tau))
	Q[:,:,1] = Q_init
	R = SharedArray(zeros(n_dim, n_dim, tau))
	R[:,:,1] = R_init
	lyap_exps = zeros(3)
	C = Matrix{Float64}(I,3,3)
    for n = 2:tau
		@sync @distributed for k = 1:n_dim
					Q[:,k,n] = tangent_step(Q[:,k,n-1], 
								u[:,n-1], s0, zeros(3))
		end
		Q[:,:,n], R[:,:,n] = qr(Q[:,:,n])
		lyap_exps += log.(abs.(diag(R[:,:,n])))/(tau*dt)
    end
	for n = tau:-1:2
		C = \(R[:,:,n], C)
    end
	return lyap_exps, C
end
