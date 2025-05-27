using RadiationPyrometers,PlanckFunctions
using Test

@testset "RadiationPyrometers.jl" begin
    p_vector = RadiationPyrometers.produce_pyrometers() # creating all default pyrometers vector
    p_new = RadiationPyrometers.Pyrometer(type="TEST2",λ=[1.1; 2.1],ϵ=0.99) # creating narrow-band pyrometer
    p_new2 = RadiationPyrometers.Pyrometer(type="TEST1",λ=18.0,ϵ=0.99) # creating narrow-band pyrometer
    push!(p_vector,p_new)
    push!(p_vector,p_new2)
    @show p_vector
    N = length(p_vector)
    Treal = 1235.0 # this is the real temperature of the surface
    Tmeasured = fill(1375.0,N) #input data, all pyrometers "measure" the same temperature
    RadiationPyrometers.fit_ϵ!(p_vector,1235.0,Tmeasured) # fitting the emissivity of each pyrometer
    Tmeasured_calc = Vector{Float64}(undef,N) #
    for (i,p) in enumerate(p_vector)
        if RadiationPyrometers.is_narrow_band(p)
            Imeas = PlanckFunctions.band_power(Treal,λₗ=p.λ[1],λᵣ=p.λ[2])*p.ϵ[]
        else
            Imeas = PlanckFunctions.ibb(p.λ[1],Treal)*p.ϵ[]
        end
        Tmeasured_calc[i] = RadiationPyrometers.measure(p,Imeas)
        println("Tmeasured_calc[$(i)] = $(Tmeasured_calc[i]); ϵ=$(p.ϵ[])")
        @test Tmeasured_calc[i] ≈ Treal rtol=1e-3 
    end
end

