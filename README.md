# RadiationPyrometers

[![Build Status](https://github.com/Manarom/RadiationPyrometers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/RadiationPyrometers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/ThermovisorData.jl)

The purpose of this small package is to create a virtual pyrometer that can be used to calculate the emissivity of a real-life pyrometer, enabling the measured temperature to be adjusted to match the actual temperature of the heated object. There are several 
default pyrometers with specified spectral ranges, see the following figure:
<p float="left">
  <img src="./notebooks/Pyrometers.png" width="400"/>
</p>
It is also possible to create self-configured pyrometer with user-specified spectral range
Notebook file with `ThermovisorData-test.jl` usage examples are available at [Pluto notebooks](https://github.com/Manarom/RadiationPyrometers.jl/blob/main/notebooks).