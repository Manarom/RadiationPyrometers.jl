### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 5e712312-0fc7-4205-84cc-834d57b814a3
begin # need this block because the PlanckFuntions is unregistered
	notebook_dir = @__DIR__()
	src_dir = joinpath(abspath(joinpath(notebook_dir,"..")),"src")
	import Pkg
	Pkg.activate(notebook_dir)
	# as far as PlanckFunctions.jl is unregistered package it should be loaded directly from the github repository using julia package manager 
	Pkg.add(url=raw"https://github.com/Manarom/PlanckFunctions.jl.git")
	Pkg.add("Revise")
	Pkg.add(["StaticArrays";"OrderedCollections";"OptimizationOptimJL";"Optimization";"Interpolations";"PrettyTables";"Plots";"LaTeXStrings";"NumericalIntegration"])
	
end

# ╔═╡ 7247c8a8-672e-4ea8-9fbb-9e38e7dea4cc
using PlutoUI,Plots,PrettyTables,Interpolations,DelimitedFiles,Revise,NumericalIntegration,LaTeXStrings

# ╔═╡ 15a5265e-61bc-440d-9a7d-ff10773b78d8
using Main.RadiationPyrometers #this line returns not defined error on the first Pluto run (probably, because of the Pluto running all "using"'s before the cells) just re-run this cell manually

# ╔═╡ 30743a02-c643-4bdc-837e-b97299f9520a
md"""
#  Pluto notebook for the demontration of `RadiationPyrometry.jl` package usage  

##### Package **`RadiationPyrometry.jl`** contains methods for virtual radiation pyrometers of different types. 
____________________
### Installation

To run this notebook, you need:
1) Install `julia` language itself from its official [download page](https://julialang.org/downloads) 
2) Install [Pluto](https://plutojl.org/) notebook from `julia` REPL by entering the following commands line-by-line:
```julia
import Pkg
Pkg.add("Pluto")
using Pluto
Pluto.run()
```
The last line will launch the Pluto starting page in your default browser 
3) Copy the entire GitHub [repository](https://github.com/Manarom/RadiationPyrometry.jl.git) to your local folder
4) Open this notebook file located at  `project_folder\notebooks\radiation_pyrometry.jl` in `Pluto` by providing the full path to the *"Open a notebook"* text field on `Pluto`'s starting page.

"""

# ╔═╡ 89a11dcd-b3b5-4349-930d-a66ad74e8fa2
import PlanckFunctions as Planck

# ╔═╡ 9cd8fe6d-dcf9-472e-a019-19b4c1a182ed
includet(joinpath(src_dir,"RadiationPyrometers.jl"))

# ╔═╡ 171409eb-22b5-4bc5-a8e2-eac0932a24f3
PlutoUI.TableOfContents(indent=true, depth=4, aside=true)

# ╔═╡ d5ee3913-66be-47d7-a755-699ba64b4f98
md"""
### Introduction

This notebook demonstrates three main examples of using the RadiationPyrometers.jl package. The purpose of this small package is to create a virtual pyrometer that can be used to calculate the emissivity of a real-life pyrometer, enabling the measured temperature to be adjusted to match the actual temperature of the heated object.  

"""

# ╔═╡ d442014a-20e6-4be4-ac7f-f13de329dec5
md"""
### I. `Blackbody` vs `real surface` thermal emission 
_______________________

All heated bodies emit thermal radiation. According to Planck's law, the spectrum of ideal emitter (called the *blackbody*) is governed solaly by its temperature. The blackbody spectral intensity can be calculated as follows \


``I_{blackbody}(\lambda , T) =  \frac{C_1}{\lambda ^5} \cdot \frac{1}{e^{\frac{C_2}{\lambda T } } - 1}``, \


where ``C_1`` = $(Planck.C₁), ``W \cdot μm/m² \cdot sr`` and ``C_2`` = $(Planck.C₂), ``μm \cdot K`` , ``\lambda`` - wavelength in ``\mu m``, ``T`` - temperature in Kelvins

A real surface thermal emission intensity is lower than the one of the blackbody. The fraction of blackbody thermal radiation intensity emitted by a real surface is characterized by directional spectral emissivity ``\epsilon (\lambda, T,\vec{\Omega})``:

``I_{real\ \ surface}(\lambda , T, \vec{\Omega} ) = \epsilon (\lambda, T,\vec{\Omega}) \cdot  I_{blackbody}(\lambda , T)``, \


here  ``\vec{\Omega}`` stays for direction. \


It is interesting that, unlike the blackbody, the real surface thermal emission (in general) depends  on the direction of radiation. Therefore, the most general characteristic for thermal radiation of a real surface is the `directional spectral emissivity`.

In [PlanckFunctions.jl](https://manarom.github.io/PlanckFunctions.jl) module there are several functions to calculate the blackbody thermal emission spectra (and various derivatives, integrals etc.).
The following figure show the impact of spectral emissivity on the real surface thermal emission intensity.
"""

# ╔═╡ 27b3c586-9eb0-4a51-b9ca-a9c0379fccdf
@bind  λ_BB PlutoUI.combine() do Child
	md"""
	Blackbody spectral range, ``\mu m`` : \
	``\lambda_{left}`` = $(
		Child(Slider(0.1:0.1:50,default=0.1,show_value = true))
	)   -- 
	 $(
		Child(Slider(0.1:0.1:50,default=15.0,show_value = true))
	)  ``\lambda_{right}`` 
	"""
end

# ╔═╡ f22d22b6-5d98-4cc4-998f-a53e92809618
@bind  T_BB PlutoUI.combine() do Child
	md"""
	Blackbody temperatures: \
	T₁ = $(
		Child(Slider(-273.0:1:2500,default=800,show_value = true))
	) ``^o C`` \
	T₂ = $(
		Child(Slider(-273.0:1:2500,default=900,show_value = true))
	)  ``^o C`` 
	"""
end

# ╔═╡ 7cc110e9-7655-4dfc-b1e0-ab3905866425
@bind  scales_BB PlutoUI.combine() do Child
	md"""
	Scales : \
	xscale = $(
		Child(Select([:identity,:ln,:log10]))
	)   yscale  
	 $(
		Child(Select([:identity,:ln,:log10]))
	)   
	"""
end

# ╔═╡ 2ad3ec82-54a2-49ac-94ef-579f808dfb1a
begin 
	rt_emissivity_data = readdlm("real_surface_emissivity.txt") # loading file, actually this file is two-column data with no headers, but this is ok
	rt_emissivity_interpolation = linear_interpolation(rt_emissivity_data[:,1],rt_emissivity_data[:,2],extrapolation_bc = Interpolations.Flat())
end;

# ╔═╡ 4d6337aa-cfc7-4154-a395-5aa53e23d01a
begin 
	T₁ = T_BB[1] + Planck.Tₖ # converting to Kelvins
	T₂ = T_BB[2] + Planck.Tₖ # converting to Kelvins
	λ_bb = collect(range(λ_BB...,length=500));
	ibb1 = Planck.ibb.(λ_bb,T₁ );ibb2 =  Planck.ibb.(λ_bb,T₂)
	
	p_bb = plot(λ_bb,ibb1,label="blackbody T=$(T₁),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.3)
	
	plot!(λ_bb,ibb2,label="blackbody T=$(T₂),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.3)
	
	plot!(λ_bb,ibb1.*rt_emissivity_interpolation(λ_bb),label="real surface T=$(T₁),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.2)	
	
	plot!(λ_bb,ibb2.*rt_emissivity_interpolation(λ_bb),label="real surface T=$(T₂),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.2)
	xlabel!("Wavelength, μm")
	ylabel!("Spectral intensity, W/m²⋅sr⋅μm")
end

# ╔═╡ 03d76e64-ebf4-432b-b9be-d4cb26275f55
begin 
	e_real_BB = rt_emissivity_interpolation(λ_bb) # this spectral emissivity was measured up to 18 μm, thus for higher wavelengths it uses flat extrapolation
	plot(λ_bb, e_real_BB, label=nothing, linewidth=3.0)
	title!("Real (measured) surface spectral emissivity")
	xlabel!("Wavelength, μm");ylabel!("Spectral emissivity (ϵ)")
end

# ╔═╡ 8a066ee5-80e9-462f-9a61-15851468aa63
md"""
### II. Partial radiation pyrometry
_______________________

As far as the `blackbody` thermal radiation energy strongly depends on temperature, this quantity can be used to measure the temperature of a real surface. This is the general idea of partial radiation pyrometry: **measure intensity to get the temperature**. As far as the intensity is a directional quantity, a pyrometer needs collimating optics (a telescope).The real surfaces emissivity often varies sufficiently with the wavelength, at the same time, partial radiation pyrometers assume constant emissivity (so-called `grey`-band approximation). Thus, for industrial purposes, it is useful to have several pyrometers, each working within a relatively narrow spectral band. In the spectral range of a partial radiation pyrometer, emissivity should not vary significantly to make the assumption of constant emissivity relevant.

The  [RadiationPyrometers.jl](https://manarom.github.io/RadiationPyrometers.jl) package provides several function to work with `virtual` partial radiation pyrometers.
"""

# ╔═╡ e2a9aa39-2490-4681-89d3-a01f058f6feb
pretty_table(HTML,RadiationPyrometers.DefaultPyrometersTypes,standalone=false,top_left_str="Table of default pyrometers types provied by `RadiationPyrometers.jl` package and corresponding wavelength regions",wrap_table_in_div=true,row_labels=["λₗ","λᵣ" ])


# ╔═╡ 9b08b767-7e8f-4483-9f2f-226022ce10e4
md"""
	It is interesting to look how various pyrometers (with their emissivity set to one) `measure` the temperature of a real surface. The following figure shows the real surface thermal emission intensity and several common pyrometers types working regions. In the legend their `measured` temperature is shown. The `mesured` temperature for each pyrometer type is obtained by fitting the blackbody power to the real surface power both integrated over pyrometer's working spectral range.  
	"""

# ╔═╡ 712828a7-fb54-42e6-95fc-233243190f59
md"Real surface temperature $(@bind T_pyr Slider(range(10,3000,1000),default=1500,show_value=true) ) "

# ╔═╡ c69acbf6-94fb-4ac3-8d56-d1f9dda11440
begin
	λ_pyr = collect(range(0.1,18,1000))
	pyrometers_vector = sort(RadiationPyrometers.produce_pyrometers())# returns a vector of all default pyrometers 
	N = length(pyrometers_vector) + 1
	# calculationg the real surface thremal radiation spectrum
	real_i = Planck.ibb.(λ_pyr,T_pyr).*rt_emissivity_interpolation(λ_pyr)
	
	data_legend = Matrix{String}(undef,N,1)
	data_legend[1] = "$(L"T_{real}") =$(round(T_pyr))"
	# poltting surface thremal emission spectrum
	plot_pyrometers = plot(λ_pyr,real_i,xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.3,dpi=600,label = data_legend[1],legend_background_color=:white,legend_foreground_color = :black,legend_position=:right)

	real_i_interp = linear_interpolation(λ_pyr,real_i)
	xlabel!("Wavelength , μm")
	ylabel!("Thermal radiation intensity")
	#ylabel!()
	max_val = maximum(real_i)
	t_em_unity = Vector{Float64}(undef,length(pyrometers_vector))
	t_em_acttual = Vector{Float64}(undef,length(pyrometers_vector))
	for (j,ppp) in enumerate(pyrometers_vector)# iterating over virtual pyrometers vector
		# the following function checks if current pyrometer is narrow band
		is_two_wavelength_pyrometer = RadiationPyrometers.is_narrow_band(ppp)
		# 
		l_cur = is_two_wavelength_pyrometer ? ppp.λ : [ppp.λ[]-0.2,ppp.λ[]+0.2 ]
		λ_pyr_interp = collect(range(l_cur...,length=30))
		# calculating the measured by the pyrometer value 
		measure_intensity =  is_two_wavelength_pyrometer ? NumericalIntegration.integrate(λ_pyr_interp,real_i_interp(λ_pyr_interp)) : real_i_interp(ppp.λ[1])
		# temperature measured by the current pyrometer
		measured_temp = round(RadiationPyrometers.measure(ppp,measure_intensity,T_starting= T_pyr))
		
		data_legend[j+1] = ppp.type*": T="*string(measured_temp)

		# plotting current pyrometer spectral range
		region_flag = 
		plot!(l_cur,[max_val,max_val],fillrange=0, fillalpha=0.5,label=data_legend[j+1])

		# remember the value of temperature with unit emissivity
		t_em_unity[j] = measured_temp 
		# calculating the averaged gray-band emissivity
		ppp.ϵ[] =measure_intensity/(is_two_wavelength_pyrometer ? Planck.band_power(T_pyr,λₗ=l_cur[1],λᵣ=l_cur[2]) : Planck.ibb(ppp.λ[],T_pyr)   )
		# calcaulting the averaged emissivity wthin the pyrometers spectral band (or at fixed wavelength)
		 t_em_acttual[j] =  round(RadiationPyrometers.measure(ppp,measure_intensity,T_starting= T_pyr))
		 
	end
	plot!(twinx(),λ_pyr,rt_emissivity_interpolation(λ_pyr),linewidth=4,linecolor=:red,label=nothing,alpha=0.3,ylabel ="Real surface spectral emissivity" )
	plot_pyrometers
end

# ╔═╡ f763d449-2a7a-4008-a183-823a774bc25e
savefig(plot_pyrometers,joinpath(notebook_dir,"Pyrometers.png"));

# ╔═╡ 667f7c30-56e0-461f-b35b-c924007eb9f2
md"""
For this particular material the type -`F` pyrometer readings are closer to the real temperature, because of the real surface emissivity being closer to unity for this pyrometer's spectral region. All results are summarized in the following table.
"""

# ╔═╡ 2d12b1a7-9474-44dc-8a39-e13bf451928e
begin # creating output table
	data = Matrix{Any}(undef,length(pyrometers_vector),5)
	e_grey = [p.ϵ[] for p in pyrometers_vector]
 	data[:,3:end] .= hcat(t_em_unity,e_grey, t_em_acttual)
	data[:,1] .= [p.type for p in pyrometers_vector]
	data[:,2] .= [string(p.λ) for p in pyrometers_vector]
	pretty_table(HTML,data,header = ["type","λ region,μm","T₀ (ϵ=1),K","grey-ϵ", "T₁ (grey-ϵ),K"],top_left_str="Table of temperatures `measured` by different pyrometers with the spectral emissivity settled to one (T₀) and to the calculated grey-ϵ and temperature measured after setting gray band emissivity to the right value (T₁) the real temperature is Tᵣ=$(T_pyr)",)
end

# ╔═╡ d08ec8f2-7689-4043-9f28-da06ab0124b9
md"""
### III. Blackbody reference source emissivity
_______________________

A real-life pyrometer does not directly measure radiation intensity due to the spectral dependence of its sensitivity and various other factors; therefore, it requires calibration. To calibrate the pyrometer, it is placed in front of a reference source that closely approximates blackbody thermal emission. Typically, this blackbody reference is a specially designed furnace with a heated cavity and a precise temperature controller.

During the calibration process, the pyrometer operator measures the temperature of the reference source. Ideally, the measured temperature should match the temperature set on the reference controller; however, in practice, these values rarely coincide exactly. This discrepancy arises from the non-ideal nature of the reference source. The spectral emissivity of a real reference is not exactly equal to one and varies with both wavelength and temperature. To account for this deviation, the reference source is supplied along with a calibration table, which looks like the following:

"""

# ╔═╡ c5ac80ee-8143-4c28-bbff-2ac761c71fac
begin
	bb_calibration_table_data = readdlm("BBethalon")
	table_header =["Tref";]
	all_types = [getfield(p,:type) for p in pyrometers_vector]
	table_header = vcat(table_header,all_types)
	pr_tbl = pretty_table(HTML,bb_calibration_table_data,header = table_header,max_num_of_rows=10,title = "Example of table data for the blackbody reference source (all tempeatures are in Celsius)")
	bb_calibration_table_data .+= Planck.Tₖ # converting table data to 
	ref_T = @view bb_calibration_table_data[:,1] # reference source temperature
	pr_tbl
end

# ╔═╡ 0c9fe7b1-374c-4fb8-9cfe-9337389713bf
md"""
The first row in the table represents the reference temperatures, while the subsequent columns show temperatures measured by different pyrometers. Each column is labeled according to the pyrometer type. It should be noted that the temperatures listed in the table vary both with the pyrometer type (across each row) and the temperature values themselves (down each column). This indicates that the spectral emissivity of the reference source depends on both the wavelength (corresponding to the pyrometer type) and the temperature. Since the various pyrometer types collectively cover a broad spectral range from 2 to 14 μm, the measured temperature data can be utilized to determine the spectral emissivity of the reference and its temperature dependence.

The following figure shows the spectral emissivity of the blackbody reference, calculated from the temperature calibration table shown above.
"""


# ╔═╡ bc2d93ae-6c30-462c-96e0-30fdb84d7c63
md"""
Adjust blackbody reference temperature $(@bind T_reference Slider(bb_calibration_table_data[:,1],show_value=true)) K

"""

# ╔═╡ d3199b6e-9779-4def-b701-fe85d1035045
begin 
	#wavlength = RadiationPyrometers.wlength
	full_wavelengths_range = RadiationPyrometers.full_wavelength_range(pyrometers_vector)
	jj = indexin(T_reference,ref_T)[]
	ej = RadiationPyrometers.fit_ϵ_wavelength!(pyrometers_vector,ref_T[jj],bb_calibration_table_data[jj,2:end])
	plot(full_wavelengths_range, ej,label=nothing, title="Spectral emissivity for T = $(ref_T[jj])")
end

# ╔═╡ Cell order:
# ╟─30743a02-c643-4bdc-837e-b97299f9520a
# ╟─5e712312-0fc7-4205-84cc-834d57b814a3
# ╟─89a11dcd-b3b5-4349-930d-a66ad74e8fa2
# ╟─7247c8a8-672e-4ea8-9fbb-9e38e7dea4cc
# ╟─9cd8fe6d-dcf9-472e-a019-19b4c1a182ed
# ╟─15a5265e-61bc-440d-9a7d-ff10773b78d8
# ╟─171409eb-22b5-4bc5-a8e2-eac0932a24f3
# ╟─d5ee3913-66be-47d7-a755-699ba64b4f98
# ╟─d442014a-20e6-4be4-ac7f-f13de329dec5
# ╟─27b3c586-9eb0-4a51-b9ca-a9c0379fccdf
# ╟─f22d22b6-5d98-4cc4-998f-a53e92809618
# ╟─7cc110e9-7655-4dfc-b1e0-ab3905866425
# ╟─4d6337aa-cfc7-4154-a395-5aa53e23d01a
# ╟─03d76e64-ebf4-432b-b9be-d4cb26275f55
# ╟─2ad3ec82-54a2-49ac-94ef-579f808dfb1a
# ╟─8a066ee5-80e9-462f-9a61-15851468aa63
# ╟─e2a9aa39-2490-4681-89d3-a01f058f6feb
# ╟─9b08b767-7e8f-4483-9f2f-226022ce10e4
# ╟─712828a7-fb54-42e6-95fc-233243190f59
# ╟─c69acbf6-94fb-4ac3-8d56-d1f9dda11440
# ╟─f763d449-2a7a-4008-a183-823a774bc25e
# ╟─667f7c30-56e0-461f-b35b-c924007eb9f2
# ╟─2d12b1a7-9474-44dc-8a39-e13bf451928e
# ╟─d08ec8f2-7689-4043-9f28-da06ab0124b9
# ╟─c5ac80ee-8143-4c28-bbff-2ac761c71fac
# ╟─0c9fe7b1-374c-4fb8-9cfe-9337389713bf
# ╟─bc2d93ae-6c30-462c-96e0-30fdb84d7c63
# ╟─d3199b6e-9779-4def-b701-fe85d1035045
