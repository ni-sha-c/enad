include("Lorenz.jl")

μ = 0.96
NLy = 50
Ly = linspace(0.5,5.0,NLy)
opt_t = zeros(NLy,1)
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

