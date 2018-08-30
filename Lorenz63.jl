#Lorenz'63 System
dt = 0.005e0
s0 = [10.,28.,8.0/3.0]
function Step(u0::Array{Float64,1},s::Array{Float64,1},N::Int64)
		u = copy(u0)
		for j = 1:N
			
			x = u[1]
			y = u[2]
			z = u[3]

			u[1] = x + dt*s[1]*(y-x)	
			u[2] = y + dt*(x*(s[2]-z) - y)
			u[3] = z + dt*(x*y - s[3]*z)	


		end

	return u
end

function tangent_step(v0::Array{Float64,1},u::Array{Float64,1},
					  s::Array{Float64,1},ds::Array{Float64,1})

	x = u[1]
	y = u[2]
	z = u[3]

    v = copy(v0)
	vx = v[1]
	vy = v[2]
	vz = v[3]

	v[1] = vx + dt*s[1]*(vy - vx) + dt*ds[1]*(y - x)
	v[2] = vy + dt*(vx*(s[2] - z) - vy - x*vz) + dt*ds[2]*x
	v[3] = vz + dt*(vx*y + vy*x - s[3]*vz) + dt*ds[3]*z
				
    return v
end

function adjoint_step(w0::Array{Float64,1},u::Array{Float64,1},
					  s::Array{Float64,1})
	
    w = copy(w0)
	wx = w[1]
	wy = w[2]
	wz = w[3]

	x = u[1]
	y = u[2]
	z = u[3]

	w[1] = wx*(1. - dt*s[1]) + 
	wy*dt*(s[2] - z) + 
	wz*dt*y 

	w[2] = wx*dt*s[1] + 
	wy*(1. - dt) + 
	wz*dt*x


	w[3] = -wy*x*dt + 
	wz*(1.0 - dt*s[3])		 

	return w


end

function objective(u::Array{Float64,1},s::Array{Float64,1})

	return u[3]	

end

function nabla_objective(u::Array{Float64,1},s::Array{Float64,1})

	return [0.,0.,1.0]	

end

function dfds(u::Array{Float64,1},
			s::Array{Float64,1},ds::Array{Float64,1})

	v = zeros(3)
	x = u[1]
	y = u[2]
	z = u[3]
	
	v[1] = ds[1]*(y-x)
	v[2] = ds[2]*x
	v[3] = ds[3]*z

	return v

end

function tangent_source(v0::Array{Float64,1},u::Array{Float64,1},
						s::Array{Float64,1},ds::Array{Float64,1})
	x = u[1]
	y = u[2]
	z = u[3]
	
    v = copy(v0)
	v[1] += dt*ds[1]*(y-x)
	v[2] += dt*ds[2]*x
	v[3] += dt*ds[3]*z


	return v

end

function adjoint_source(w0::Array{Float64,1},u::Array{Float64,1},
			s::Array{Float64,1},dm::Float64)

    w = copy(w0)
	dJ = nabla_objective(u,s)
	w += dm*dJ
    return w 

end

#=
function gradFs(u::Array{Float64,1},s::Array{Float64,1})

	x = u[1]
	y = u[2]
	z = u[3]

	A = zeros(3,3)
	A[1,1] = -s[1]
	A[1,2] = s[1]
	A[2,1] = s[2]-z
	A[2,2] = -1
	A[2,3] = -x
	A[3,1] = y 
	A[3,2] = x
	A[3,3] = -s[3]

	return A

end
=#

