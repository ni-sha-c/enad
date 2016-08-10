include("Lorenz.jl")

nop = 200
M = linspace(10,5000,nop)
θ = zeros(nop,1)
μ = 0.96
for i = 1:nop

	θ[i] = Lorenz.get_dzbardr_EA(floor(Int64,M[i]))

end
err = log(abs(θ - μ))
plot(M,err)
