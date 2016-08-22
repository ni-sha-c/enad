include("Lorenz.jl")

nop = 50
N = linspace(500,40000,nop)
M = 90
θ = zeros(nop,1)
var_sens = zeros(nop,1)
μ = 0.96
AD_OUT = Lorenz.get_dzbardr_EA_variableCost(floor(Int64,N[end]),M)
for i = 1:nop

	θ[i] = mean(AD_OUT[1:floor(Int64,N[i]/M)])
	var_sens[i] = var(AD_OUT[1:floor(Int64,N[i]/M)])
	
end

