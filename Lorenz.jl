#Lorenz'63 System
module Lorenz

function f(X::Array{Float64,1},r::Float64)
	params = get_system_params()
	σ = params[1]
	b = params[2]
	x = X[1]
	y = X[2]
	z = X[3]
	dXdt = [-σ*x + σ*y, 
			-x*z + r*x - y,
			x*y - b*z]

	return dXdt 
end
function get_system_params()
	σ = 10.e0
	b = 8.0/3.0
	return [σ, b]

end
function get_r()
		r = 28.
end
function get_dt()
	dt = 0.005e0
	return dt
end
function get_dr()
	dr = 0.005e0
	return dr
end
function get_N()
#Long integration time
	N = 26272
	return N
end
function get_n(M::Int64)
#Number of Random Trajectories
	N = get_N()
	n = floor(Int64,N/M)
	return n
end
function get_zbar_direct(r::Float64,N::Int64=get_N(),
		X0::Array{Float64,1}=[-2.4e0,-3.7e0,14.98e0],
		flag::Int64=1)	
# zbar = (1/τ) ∫_0^τ z(t) dt

		dt = get_dt()
		dr = get_dr()
		x0 = X0[1]
		y0 = X0[2]
		z0 = X0[3]
		τ  = N*dt

		
		X0 = [x0, y0, z0]
		X = zeros(3,N)
		X[:,1] = X0
		for j = 2:N
				Xjm1 = X[:,j-1]
				k1 = dt*f(Xjm1,r)
				k2 = dt*f(Xjm1+0.5*k1,r)
				k3 = dt*f(Xjm1+0.5*k2,r)
				k4 = dt*f(Xjm1+k3,r)
				X[:,j] = Xjm1 + 1./6.*k1 + 
						 1./3.*k2 + 1./3.*k3 +
						 1./6.*k4 
		end
		
		zbar = (0.5*X[3,1] + sum(X[3,2:end-1]) + 0.5*X[3,end])/N

		if(flag==1)
			return zbar
		end
		return X
end
#=
		w1 = [0, x0*dt, 0]
		w2 = [0, 0, 1]
		θ  = zeros(m)
	 	for i = 1:m
			for k = 1:n
				
				θ[i] += w1'*J(k,X0)*w2       
	
			end
			θ[i] /= n
		end

		print(mean(θ))
=#

function Jcbn(Xnm1::Array{Float64,1},r)
		
	params = get_system_params()
	σ = params[1]
	b = params[2]

	dt = get_dt()
	Ja = zeros(3,3)

	x = Xnm1[1]
	y = Xnm1[2]
	z = Xnm1[3]

	Ja[1,1] = 1. - σ*dt
	Ja[1,2] = σ*dt
	Ja[1,3] = 0.0

	Ja[2,1] = (r-z)*dt
	Ja[2,2] =  1. - dt
	Ja[2,3] = -x*dt

	Ja[3,1] = y*dt
	Ja[3,2] = x*dt
	Ja[3,3] = 1. - b*dt

	return Ja
end
function get_dzbardr_adjoint(X1::Array{Float64,1},M::Int64)


		r = get_r()
		dt = get_dt()
		dzndX2 = zeros(3,M)
		X = get_zbar_direct(r,M,X1,0)

		for j=2:M
		    dzjdX2 = [0,0,1]
			for k=1:j-2
				dzjdX2 = (Jcbn(X[:,j-k],r))'*dzjdX2
			end
			dzndX2[:,j] = dzjdX2
		end

		dX2dr = [0.  X[1,1]*dt  0.]
		dzndr = dX2dr*dzndX2

		dzbardr = (0.5*dzndr[1] + 0.5*dzndr[end] +
					sum(dzndr[2:end-1]))/M

		return dzbardr

end
function get_dzbardr_EA(M::Int64)

	n = get_n(M)		
	x0 = -15.e0 + 30.e0*rand(1,n)
	y0 = -25.e0 + 50.e0*rand(1,n)
	z0 = 5.0 + 35.0*rand(1,n)
	
	dzbardr = 0.0
	for i = 1:n
			X0 = [x0[i],y0[i],z0[i]]
			dzbardr += get_dzbardr_adjoint(X0,M)
	end
	dzbardr /= n
	return dzbardr
end
end
