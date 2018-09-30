@everywhere using LinearAlgebra
@everywhere include("Lorenz63.jl")
@everywhere using SharedArrays
function compute_CLVs(u0::Array{Float64,1})
    u0 = Step(rand(3),s0,1000)
    n_dim = length(u0)
	du = rand(n_dim, n_dim)
	du /= norm(du)
	Q_init, R_init = qr(du)
	t = 1000
	tau = 600
	u = zeros(n_dim,t)
	u[:,1] = copy(u0)
	for n = 2:t
		u[:,n] = Step(u[:,n-1],s0,1)
	end
	Q = SharedArray(zeros(n_dim, n_dim, t))
	Q[:,:,1] = Q_init
	R = SharedArray(zeros(n_dim, n_dim, t))
	R[:,:,1] = R_init
	lyap_exps = zeros(3)
	C = Matrix{Float64}(I,3,1)
    for n = 2:t
		@sync @distributed for k = 1:n_dim
					Q[:,k,n] = tangent_step(Q[:,k,n-1], 
								u[:,n-1], s0, zeros(3))
		end
		Q[:,:,n], R[:,:,n] = qr(Q[:,:,n])
		lyap_exps += log.(abs.(diag(R[:,:,n])))/(t*dt)
    end
	for n = t:-1:100
		C = \(R[:,:,n], C)
		C = C'.\diag(C)
		C = C'
		println(C[:,1]')
    end
	return lyap_exps, C, Q, R
end
