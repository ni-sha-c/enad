@everywhere include("test_ensemble_sensitivity.jl")
using DelimitedFiles

u0, factor_exp_thetaN = test_rare_event(1000,"adjoint")
open("data/scale_sensitivities.bin","w") do io
    writedlm(io, factor_exp_thetaN)
end

open("data/initial_conditions.bin","w") do io
    writedlm(io, u0)
end

