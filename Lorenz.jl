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
function dfdr(X::Array{Float64,1},r::Float64)

		return [0, X[1], 0]
end
function dfdX(X::Array{Float64,1},r::Float64)
	
	dfdX_res = zeros(3,3)
	params = get_system_params()
	σ = params[1]
	b = params[2]
	x = X[1]
	y = X[2]
	z = X[3]

	dfdX_res[1,1] = -σ
	dfdX_res[1,2] = σ
	dfdX_res[2,1] = -z + r
	dfdX_res[2,2] = -1.0
	dfdX_res[2,3] = -x
	dfdX_res[3,1] = y
	dfdX_res[3,2] = x
	dfdX_res[3,3] = -b

	return dfdX_res

end
function dvdt(X::Array{Float64,1},v::Array{Float64,1},r::Float64)
	
		return dfdr(X,r) + dfdX(X,r)*v
		
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

		λ = zeros(3,M)
		c = [0.0, 0.0, 1./M]
		b = zeros(3*M,1)
		λ[:,M] = [0,0,1./M]
		#A = eye(3*M,3*M)
		
		for i=1:M-1
			
			b[i*3+1:(i+1)*3] = dfdr(X[:,i+1],r)*dt	
				
			λ[:,M-i] = c + (eye(3,3) + dt*dfdX(X[:,M-i],r))'*λ[:,M-i+1]	

			#A[3*i + 1:(i+1)*3,(i-1)*3+1:i*3] = -eye(3,3) - dt*dfdX(X[:,i],r)
			
			
		end

		λ = λ[:]
		dzbardr = λ'*b
		#v = A\b


		return dzbardr

end
function get_dzbardr_tangent(X1::Array{Float64,1},M::Int64,flag::Int64=0)

	r = get_r()
	dt = get_dt()

	v = zeros(3,M)
	X = get_zbar_direct(r,M,X1,0)


	for i =2:M
		
		k1 = dt*dvdt(X[:,i-1], v[:,i-1], r)
		k2 = dt*dvdt(X[:,i-1], v[:,i-1] + 0.5*k1, r)
		k3 = dt*dvdt(X[:,i-1], v[:,i-1] + 0.5*k2, r)
		k4 = dt*dvdt(X[:,i-1], v[:,i-1] + k3, r)
		v[:,i] = v[:,i-1] + 1./6.*k1 + 1./3.*k2 + 1./3.*k3 +
				1./6.*k4


		
	end
	dzbardr_tangent = 1./M*sum([0., 0., 1.0]'*v)
	
	if(flag==0)
		return dzbardr_tangent
	end

	return v
	
end
function get_dzbardr_EA(M::Int64)

	n = get_n(M)		
	x0 = -15.e0 + 30.e0*rand(1,n)
	y0 = -25.e0 + 50.e0*rand(1,n)
	z0 = 5.0 + 35.0*rand(1,n)
	
	dzbardr = zeros(n)
	for i = 1:n
			X0 = [x0[i],y0[i],z0[i]]
			dzbardr[i] = get_dzbardr_adjoint(X0,M)[1,1]
	end
	dzbardr_mean = mean(dzbardr)
	dzbardr_var = var(dzbardr)
	return dzbardr_mean,dzbardr_var
end
end
