include("Lorenz.jl")
nop = 50
N = linspace(500,40000,nop)
N = floor(Int64,N)
θ = zeros(nop,1)
var_sens = ones(nop,1)
ms_err = ones(nop,1)

μ = 0.96
NLy = 50
Ly = linspace(0.5,5.0,NLy)
opt_t = zeros(NLy,1)
err_tol = 1.e-6
for k = 1:NLy
	
	
	M = 1:5:floor(Int64,600/Ly[k])
	nop = length(M)
	θ = zeros(nop,1)
	var_sens = zeros(nop,1)

	for i = 1:nop

	AD_OUT = Lorenz.get_dzbardr_EA(floor(Int64,M[i]),Ly[k])
	θ[i] = AD_OUT[1]

	var_sens[i] = AD_OUT[2]
end
	opt_t[k] = findmin(abs(θ-μ))[2]
end
μ = 0.96
d = floor(Int64,N[2]-N[1])
M = floor(Int64,N[end])
nos = floor(Int64,N[end]/M)
M_N = ones(nop,1)
while(nos < 1000)
	#M= 90
	AD_OUT = Lorenz.get_dzbardr_EA_variableCost(floor(Int64,N[end]),M)									
	#istart = ceil(Int64, (M- 500)/d + 1) 	
	for i=1:nop			
		if(N[i]>M)	
			θ1 = mean(AD_OUT[1:floor(Int64,N[i]/M)])
			var_sens1 = var(AD_OUT[1:floor(Int64,N[i]/M)])
			ms_err1 = (θ1 -μ)^2.0 + var_sens1 
			#	print(ms_err1)
			if(ms_err1 < ms_err[i])
			#if((abs(θ1-μ) < abs(θ[i]-μ))&&(var_sens1 < var_sens[i]))
				θ[i] = θ1
				var_sens[i] = var_sens1
				M_N[i] = M
				ms_err[i] = ms_err1

			end
		end
	end

	nos += 1
	M = floor(Int64,N[end]/nos)
	print(M," ")
end

