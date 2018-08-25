#Lorenz'63 System
module Lorenz

function f(X::Array{Float64,1},r::Float64,t::Float64=1.)
	params = get_system_params()
	σ = params[1]
	b = params[2]
	x = X[1]
	y = X[2]
	z = X[3]
	dXdt = [t*(-σ*x + σ*y), 
			t*(-x*z + r*x - y),
			t*(x*y - b*z)]

	return dXdt 
end
function dfdr(X::Array{Float64,1},r::Float64,t::Float64=1.)

		return [0, t*X[1], 0]
end
function dfdX(X::Array{Float64,1},r::Float64,t::Float64=1.)
	
	dfdX_res = zeros(3,3)
	params = get_system_params()
	σ = params[1]
	b = params[2]
	x = X[1]
	y = X[2]
	z = X[3]

	dfdX_res[1,1] = -t*σ
	dfdX_res[1,2] = t*σ
	dfdX_res[2,1] = t*(-z + r)
	dfdX_res[2,2] = -t*1.0
	dfdX_res[2,3] = -t*x
	dfdX_res[3,1] = t*y
	dfdX_res[3,2] = t*x
	dfdX_res[3,3] = -t*b

	return dfdX_res

end
function dvdt(X::Array{Float64,1},v::Array{Float64,1},r::Float64,t::Float64=1.)
	
		return dfdr(X,r,t) + dfdX(X,r,t)*v
		
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
	N = 1000000
#	N = 100000
	return N
end
function get_n(M::Int64,N::Int64=get_N())
#Number of Random Trajectories
#	if(M>100) 
		#	N = 10000000
#	end
	n= floor(Int64,N/M)
	
	return n
end
function get_zbar_direct(r::Float64,M::Int64,
		X0::Array{Float64,1}=[-2.4e0,-3.7e0,14.98e0],
		flag::Int64=1,t::Float64=1.,N::Int64=get_N())	
# zbar = (1/τ) ∫_0^τ z(t) dt

		dt = get_dt()
		dr = get_dr()
		x0 = X0[1]
		y0 = X0[2]
		z0 = X0[3]
		τ  = N*dt

		
		X0 = [x0, y0, z0]
		Xj = X0
	    if(flag!=2)	
			X = zeros(3,M)
			N = M
			X[:,1] = X0
		end
		for j = 2:N
				Xjm1 = Xj
				k1 = dt*f(Xjm1,r,t)
				k2 = dt*f(Xjm1+0.5*k1,r,t)
				k3 = dt*f(Xjm1+0.5*k2,r,t)
				k4 = dt*f(Xjm1+k3,r,t)
				Xj = Xjm1 + 1./6.*k1 + 
						 1./3.*k2 + 1./3.*k3 +
						 1./6.*k4
				if(flag!=2)	
					X[:,j] = Xj
				end
		end
		
		#zbar = (0.5*X[3,1] + sum(X[3,2:end-1]) + 0.5*X[3,end])/N

		#if(flag==1)
		#	return zbar
		#end
		if(flag==2)
			return Xj
		end
		return X
end
function get_dzbardr_adjoint(X1::Array{Float64,1},M::Int64,t::Float64=1.,flag::Int64=1)


		r = get_r()
		dt = get_dt()
		X = get_zbar_direct(r,M,X1,0,t)

		λ = zeros(3,M-1)
		c = [0.0, 0.0, 1./M]
		b = zeros(3*(M-1),1)
		λ[:,M-1] = [0,0,1./M]
		A = eye(3*(M-1),3*(M-1))
		b[1:3] = dfdr(X[:,1],r,t)*dt	
		for i=1:M-2
			
			b[i*3+1:(i+1)*3] = dfdr(X[:,i+1],r,t)*dt	
				
			λ[:,M-i-1] = c + (eye(3,3) + dt*dfdX(X[:,M-i],r,t))'*λ[:,M-i]	

			A[3*i + 1:(i+1)*3,(i-1)*3+1:i*3] = -eye(3,3) - dt*dfdX(X[:,i+1],r)
		
			
		end

		λ = λ[:]
		dzbardr = λ'*b
		v = A\b

		if(flag==1)
			return dzbardr
		end
		return v

end
function get_dzbardr_adjoint_2(X1::Array{Float64,1},M::Int64,t::Float64=1.,flag::Int64=1)


		r = get_r()
		dt = get_dt()
		X = get_zbar_direct(r,M,X1,0,t)

		λ = zeros(3,M-1)
		c = [0.0, 0.0, 1./M]
		b = zeros(3*(M-1),1)
		A = eye(3*(M-1),3*(M-1))
		b[1:3] = dfdr(X[:,1],r,t)*dt	
		λ[:,M-1] = (eye(3,3) - dt*dfdX(X[:,M],r,t))'\c
		for i=1:M-2
			
			b[i*3+1:(i+1)*3] = dfdr(X[:,i+1],r,t)*dt	
				
			λ[:,M-i-1] = (eye(3,3) - dt*dfdX(X[:,M-i],r,t))'\(c + λ[:,M-i])	

			A[3*i + 1:(i+1)*3,(i-1)*3+1:i*3] = -eye(3,3) - dt*dfdX(X[:,i+1],r)
		
			
		end

		λ = λ[:]
		dzbardr = λ'*b
		v = A\b

		if(flag==1)
			return dzbardr
		end
		return λ

end

function get_dzbardr_adjoint_3(X1::Array{Float64,1},M::Int64,t::Float64=1.,flag::Int64=1)


		r = get_r()
		dt = get_dt()
		X = get_zbar_direct(r,M,X1,0,t)

		λ = zeros(3,M-1)
		c = [0.0, 0.0, 1.0]
		b = zeros(3*(M-1),1)
	
		b[1:3] = dfdr(X[:,1],r,t)*dt	
		λ[:,M-1] = (eye(3,3) - dt*dfdX(X[:,M],r,t))'\c
		for i=1:M-2
			
			b[i*3+1:(i+1)*3] = dfdr(X[:,i+1],r,t)*dt	
				
			λ[:,M-i-1] = (eye(3,3) - dt*dfdX(X[:,M-i],r,t))'\(-dt*c + λ[:,M-i])	

		
		
			
		end

		λ = λ[:]
		dzbardr = -λ'*b
	

		if(flag==1)
			return dzbardr
		end
		return λ

end



function get_dzbardr_tangent(X1::Array{Float64,1},M::Int64,flag::Int64=0,
		t::Float64=1.)

	r = get_r()
	dt = get_dt()

	v = zeros(3,M)
	X = get_zbar_direct(r,M,X1,0,t)


	for i =2:M
		
		k1 = dt*dvdt(X[:,i-1], v[:,i-1], r,t)
		k2 = dt*dvdt(X[:,i-1], v[:,i-1] + 0.5*k1, r,t)
		k3 = dt*dvdt(X[:,i-1], v[:,i-1] + 0.5*k2, r,t)
		k4 = dt*dvdt(X[:,i-1], v[:,i-1] + k3, r,t)
		v[:,i] = v[:,i-1] + 1./6.*k1 + 1./3.*k2 + 1./3.*k3 +
				1./6.*k4


		
	end
	dzbardr_tangent = 1./M*sum([0., 0., 1.0]'*v)
	
	if(flag==0)
		return dzbardr_tangent
	end

	return v
	
end
function get_dzbardr_EA_fixedCost(M::Int64,T::Int64,t::Float64=1.)
	#number of trajectories in one experiment: n
	#total number of experiments: k
	n = floor(Int64,T/M)
	k = 100
	x0 = -15.e0 + 30.e0*rand(1,n*k)
	y0 = -25.e0 + 50.e0*rand(1,n*k)
	z0 = 5.0 + 35.0*rand(1,n*k)

	
	dzbardr = zeros(n*k)
	for i = 1:n*k
			X0 = [x0[i],y0[i],z0[i]]
			X = get_zbar_direct(get_r(),100,X0,2,1.,5000)	

			X0 = X
		
	#		dzbardr[i] = get_dzbardr_adjoint_3(X0,M,t)[1,1]
    		dzbardr[i] = get_dzbardr_tangent(X0,M,0,t)#[1,1]
	end
	dzbardr_mean = mean(dzbardr)
	dzbardr_var = 0.0
	for j=1:k
			dzbardr_var += (mean(dzbardr[(j-1)*n+1:j*n])-dzbardr_mean).^2.0

			#dzbardr_var += (mean(dzbardr[(j-1)*n +1:j*n])-dzbardr_mean)^2.0
	end
	dzbardr_var = dzbardr_var/k

	return dzbardr_mean,dzbardr_var
	return X0
end
function get_dzbardr_EA_variableCost(N::Int64,M::Int64,t::Float64=1.)

	n = get_n(M,N)		
	x0 = -15.e0 + 30.e0*rand(1,n)
	y0 = -25.e0 + 50.e0*rand(1,n)
	z0 = 5.0 + 35.0*rand(1,n)
	
	dzbardr = zeros(n)
	for i = 1:n
			X0 = [x0[i],y0[i],z0[i]]
			dzbardr[i] = get_dzbardr_adjoint(X0,M,t)[1,1]
	end
	#dzbardr_mean = mean(dzbardr)
	#dzbardr_var = var(dzbardr)
	return dzbardr
end

end
