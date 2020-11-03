# MoistAir - Thermodynamic properties of moist air

[![Build Status](https://github.com/longemen3000/WaterIF97.jl/workflows/CI/badge.svg)](https://github.com/longemen3000/WaterIF97.jl/actions)
[![Build Status](https://travis-ci.com/longemen3000/WaterIF97.jl.svg?branch=master)](https://travis-ci.com/longemen3000/WaterIF97.jl)

This package provides a model, `IndustrialWater` to compute thermodynamic properties of water, using the [ThermoState](https://github.com/longemen3000/ThermoState.jl) interface. The model uses the correlations described by IAPWS IF-97



This is a fork of [SteamTables.jl](https://github.com/braamvandyk/SteamTables.jl).

## instalation

```julia-repl
julia> ]
(v1.5) pkg> add https://github.com/longemen3000/WaterIF97.jl
```
## User interface

The package exports a single model: `IndustrialWater`, that accepts the following `ThermoState` specifications:

 - Pressure-Temperature (`state(p=p0,t=t0)`)
 - Pressure-Enthalpy (`state(p=p0,h=h0)`)
 - Pressure-Entropy (`state(p=p0,s=s0)`)
 - Saturation-Temperature-Phase (`state(sat=true,t=t0,phase=:liquid)`)
 - Saturation-Pressure-Phase (`state(sat=true,p=p0,phase=:liquid)`)


## Usage

```julia
using ThermoState,Unitful,WaterIF97
water = IndustrialWater() #Water model
λ = VariableSpec() #to create a variable model
st_t = state(p=1u"atm",t=λ,moles=2.5)
mol_enthalpy(water,st_t(65u"°C"))
mass_enthalpy(water,st_t(65u"°C"))
total_enthalpy(water,st_t(65u"°C"))

satv_t = state(sat=true,t=λ,moles=2.5,phase=:gas)

mol_volume(water,satv_t(100u"°C"))
```
