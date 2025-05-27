push!(LOAD_PATH,"../src/")
#include("../src/RadiationPyrometers.jl")
using Documenter, RadiationPyrometers
mathengine = Documenter.MathJax3()
makedocs(
        sitename = "RadiationPyrometers.jl",
        highlightsig = false,
        checkdocs = :none,
        format=Documenter.HTML(size_threshold = 2000 * 2^10),
        pages=[
                "Home" => "index.md"
                "RadiationPyrometers" => "pyrometers.md"
               ]#
			   )
#deploydocs(;
#         repo="github.com/Manarom/BandPyrometry"
#)