include("Lorenz.jl")
using PyPlot
nop = 200
M = linspace(3,2000,nop)
θ = zeros(nop,1)
var_sens = zeros(nop,1)
μ = 0.96
for i = 1:nop
	print("i= ",i)	
	AD_OUT = Lorenz.get_dzbardr_EA_fixedCost(floor(Int64,M[i]))
	θ[i] = AD_OUT[1]
	var_sens[i] = AD_OUT[2]
end
err = log(abs(θ - μ))
#plot(M,err)
