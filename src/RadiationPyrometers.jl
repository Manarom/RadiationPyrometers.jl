
module RadiationPyrometers
    using   Optimization,
            OptimizationOptimJL,
            LinearAlgebra,
            StaticArrays,
            OrderedCollections 
    import  PlanckFunctions as Planck
    export Pyrometer,
        DefaultPyrometersTypes,
        fit_ϵ!,
        fit_ϵ_wavelength!
    """
    Default pyrometers types 

"""
const DefaultPyrometersTypes = OrderedDict(
                    "P"=> SVector{2}([2.0; 2.6]),
                    "M"=> SVector{1}([3.4]), 
                    "D"=> SVector{1}([3.9]),
                    "L"=> SVector{1}([4.6]),
                    "E"=> SVector{2}([4.8, 5.2]),
                    "F"=> SVector{1}([7.9]),
                    "K"=> SVector{2}([8.0, 9.0]),
                    "B"=> SVector{2}([9.1,14.0])
    )
    
    struct Pyrometer{N} # this type supports methods for radiative pyrometers
        type::String
        λ::SVector{N,Float64}
        ϵ::Base.RefValue{Float64}
        """
    Pyrometer(type::String)

Pyrometer object Constructor, the type of pyrometer can chosen from the DefaultPyrometersTypes dictionary
Input:
type - pyrometer type, must be member of DefaultPyrometersTypes 
"""
        Pyrometer(type::String) = begin
            if haskey(DefaultPyrometersTypes,type) 
                N = length(DefaultPyrometersTypes[type])
                return new{N}(type,
                           DefaultPyrometersTypes[type],
                           Ref(1.0)) 
            else
                 error("Unknown pyrometer type")
            end
        end
        Pyrometer(;type::String,λ::Union{Vector{Float64},Float64},ϵ::Float64) = begin
            @assert 0<ϵ<=1.0 "Emissivity should be within the (0..1] interval"
            if λ isa Vector
                N = length(λ)
                @assert N==1 || N==2 "λ should be a vector of two  Floats or a single Float number"
            else
                N=1
            end
            new{N}(type,SVector{N}(λ),Ref(ϵ))
        end
    end
    """
    wlength(::Pyrometer{N}) where N

Returns the number of wavelengths
"""
wlength(::Pyrometer{N}) where N = N
    """
    is_narrow_band(p::Pyrometer)

True if pyrometer `p` is a narrow-band pyrometer (worsk on a fixed wavelengh region)
"""
is_narrow_band(p::Pyrometer)  = wlength(p) == 2
"""
    is_fixed_wavelength(::Pyrometer{N}) where N

True if Pyrometer is single wavelength
"""
is_fixed_wavelength(::Pyrometer{N}) where N = N==1
    """
    measure(p::Pyrometer,i::Float64)

Calculates the "measured" temperature from "mesaured" intensity by fitting the Planck function.
The intensity units should be consistent with PlanckFunctions.ibb(λ,T) function for single wavelength pyrometer 
and with PlanckFunctions.band_power(λ,T),
it should be in [W/m²⋅sr⋅μm]
Input:
p - pyrometer object
i - measured intensity in [W/(m²⋅sr⋅μm)] or integral (over wavelength spectral intensity) in [W/(m²⋅sr)] 
returns the temperature "measured" by this pyrometer  
(optional)
T_starting  - starting temperature value
"""
    function measure(p::Pyrometer,i::Float64;T_starting::Float64=600.0)
        tup = (p,i)
        if !is_narrow_band(p) # single wavelength pyrometer
            fun = OptimizationFunction((t,tup)-> norm(tup[1].ϵ[]*Planck.ibb(tup[1].λ[1],t[]) - tup[2]))
        else# narrow band pyrometer
            fun = OptimizationFunction((t,tup)-> norm(tup[1].ϵ[]*Planck.band_power(t[],λₗ=tup[1].λ[1],λᵣ=tup[1].λ[2]) - tup[2]))
        end
        prob = OptimizationProblem(
                fun, #fun
                [T_starting],#starting temperature
                tup)
        sl = solve(prob,NelderMead())
        return sl.u[]   
    end
    """
    Base.isless(p1::Pyrometer,p2::Pyrometer)

Vector of Pyrometer objects can be sorted using isless
"""
function Base.isless(p1::Pyrometer,p2::Pyrometer) # is used to sort the vector of pyrometers
        return all(p1.λ.<p2.λ)
    end
    """
    wavelength_number()

Returns the length of wavelengths vector covered by the [`DefaultPyrometersTypes`](@ref) (all default pyrometers wavelengh region)
"""
function wavelengths_number()
        return mapreduce(x->length(x),+,DefaultPyrometersTypes)
    end
"""
    wavelengths_number(p::Vector{Pyrometer})

Returns the total number of wavelength for the vector of pyrometers
"""
function wavelengths_number(p::Vector{Pyrometer})
    return sum(wlength,p)
end
    """
    full_wavelength_range()

Creates the wavelengths vector covered by default pyrometers see [`DefaultPyrometersTypes`](@ref) 
"""
function full_wavelength_range()
        #sz = mapreduce(x->length(x),+,DefaultPyrometersTypes)
        λ = Vector{Float64}()
        pyr_names = Vector{String}()
        for l in DefaultPyrometersTypes
            append!(λ,l[2])
            push!(pyr_names,l[1])
            if length(l[2])>1
                push!(pyr_names,l[1])
            end
        end
        inds = sortperm(λ)
        pyr_names.=pyr_names[inds]
        λ.=λ[inds]
        return λ,pyr_names
    end
    """
    full_wavelength_range()

Creates the wavelengths vector covered by all pyrometers in vector `p`
"""
function full_wavelength_range(p::Vector{Pyrometer})
        #sz = mapreduce(x->length(x),+,DefaultPyrometersTypes)
        λ = Vector{Float64}(undef,wavelengths_number(p))
        counter = 0
        for pj in p
            if is_narrow_band(pj)
                counter += 2
                λ[counter-1] = pj.λ[1]
            else
                counter+=1
            end    
            λ[counter] = pj.λ[end]
        end
        return λ
    end   
    """
    produce_pyrometers()

Creates the vector of all default pyrometers 
"""
function produce_pyrometers()
        pyr_vec = Vector{Pyrometer}()
        for l in DefaultPyrometersTypes
            push!(pyr_vec,Pyrometer(l[1]))
        end
        sort!(pyr_vec) # sorting according to the wavelength increase
        return pyr_vec
    end
    """
    set_emissivity(p::Pyrometer,em_value::Float64)

Setter for spectral emissivity
"""
function set_emissivity(p::Pyrometer,em_value::Float64)
        if !(0.0<em_value<=1.0)
            em_value = 1.0
        end
        p.ϵ[]=em_value
    end

    """
    fit_ϵ!(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})

Fits the emissivity of pyrometers to make measured temperature `Tmeasured` fit
fit the real temperature `Treal`

Input:
p - pyrometer objects vector , [Nx0]
Treal - real temperature of the surface, Kelvins
Tmeasured - temperatures measured by the pyrometers, in Kelvins, [Nx0]

"""
function fit_ϵ!(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})
        @assert length(p)==length(Tmeasured)  "Vectors must be of the same size"
        N = length(p)
        e_out = Vector{Float64}(undef,N)
        Threads.@threads for i in 1:N
            @inbounds e_out[i] = fit_ϵ!(p[i],Tmeasured[i],Treal)
        end
        return e_out
    end

"""
    fit_ϵ_wavelength!(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})

The same as [`fit_ϵ!`](@ref) except that it returns the vector of fitted emissivities 
of the same length to the total number of wavelength in all pyrometers in vaector `p`,
e.g. if p[i] is the narrow-band pyrometer 
"""
function fit_ϵ_wavelength!(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})  
        total_wavelength_number =  sum(wlength,p)
        e_out= Vector{Float64}(undef,total_wavelength_number)
        counter = 0
        for (i,e) in enumerate(fit_ϵ!(p,Treal,Tmeasured))
            if is_narrow_band(p[i]) 
                counter +=2 
                e_out[counter-1] = e
            else
                counter +=1 
            end
            e_out[counter] = e
        end
        return e_out
    end
    """
    fit_ϵ!(p::Pyrometer,Tmeasured::Float64,Treal::Float64)

    Optimizes the emissivity of the pyrometer to make measured by the pyrometer temperature fit
    fit the real temperature

    Input:
    p - pyrometer object
    Treal - real temperature of the surface, Kelvins
    Tmeasured - temperature measured by the pyrometer, Kelvins
"""
function fit_ϵ!(p::Pyrometer,Tmeasured::Float64,Treal::Float64;optimizer = NelderMead())
        # fits the emssivity of the pyrometer
        # Tmeasured - is the temperature measured by the pyrometer 
        # Treal  - is the real temperature of the surface
        # the real temperature of the surface should be higher due to the emissivity
        if Treal<Tmeasured
            (Treal,Tmeasured)=(Tmeasured,Treal)
        end
        # tup = (p,Tmeasured,Treal)
        if length(p.λ)==1 # for pyrometers with fixed wavelength
            fun = OptimizationFunction((ϵ,tup)->norm(Planck.ibb(tup[1].λ[1],tup[2]) - ϵ[1]*Planck.ibb(tup[1].λ[1],tup[3])),AutoForwardDiff())
        else # pyrometer with limited band
            fun = OptimizationFunction((ϵ,tup)->norm(Planck.band_power(tup[2],λₗ=tup[1].λ[1],λᵣ=tup[1].λ[2]) .- ϵ[1]*Planck.band_power(tup[3],λₗ=tup[1].λ[1],λᵣ=tup[1].λ[2])),AutoForwardDiff())
        end
        prob = OptimizationProblem(
            fun, #fun
           [0.1],#starting vectors
           (p,Tmeasured,Treal),
           lb=[0.00],
           ub=[1.0])
        res = solve(prob,optimizer)
        set_emissivity(p,res.u[1])
        return p.ϵ[]
    end
    """
    fit_ϵ_wavelength!(p::Pyrometer,Tmeasured::Float64,Treal::Float64)

    Fits emissivity and returns it as a vector of the same size as pyrometer's wavelength region
    Some pyrometers has 2-wavelength, other work on a single wavelength, for two-wavelength pyrometers
    e_out will be a two-element vector
    Input:
        p - pyrometer object
        Tmeasured - temperature measured by the pyrometer, Kelvins
        Treal - real temeprature of the surface, Kelvins
     
"""
function fit_ϵ_wavelength!(p::Pyrometer,Tmeasured::Float64,Treal::Float64) # this is the same as fit_ϵ! with the exception that 
        # this function returns a vector of values, if pyrometer is single wavelength it returns one -element array
        e_out = similar(p.λ)
        e_out .=fit_ϵ!(p,Tmeasured,Treal)
        return e_out
    end
    """
    switch_the_type(λ::Float64)

    Returns the type (which can be used as a key of DefaultPyrometersTypes dict) depending on the wavelengh
    Input:
        λ - wavelength in μm
"""
function switch_the_type(λ::Float64)
        for (k,λp) in DefaultPyrometersTypes
            if length(λp)==1
                !isapprox(λ,λp[1],atol=0.05) ? continue : return k
            else
               !(λp[1]<=λ<=λp[2]) ? continue : return k
            end
        end 
        return ""
    end

    function Base.show(io::IO, p::Pyrometer{N}) where N 
        if N==2
            print(io, "Narrow-band pyrometer:λ ∈ $(p.λ[1]) ... $(p.λ[1]) μm,ϵ = $(p.ϵ[])")
        elseif N==1
            print(io, "Fixed-wavelength pyrometer:λ = $(p.λ[1]) μm,ϵ = $(p.ϵ[])")
        else
            print(io,p)
        end
    end
end