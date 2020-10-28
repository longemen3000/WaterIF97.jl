module WaterIF97

    using Unitful
    using ThermoState
    using Roots

    using ThermoState.Types
    using ThermoState.QuickStates
    using ThermoState.StatePoints

    import ThermoState: pressure,temperature,mass,moles,molar_mass
    import ThermoState: mass_volume, mol_volume, total_volume
    import ThermoState: mass_density, mol_density
    import ThermoState: mass_enthalpy, mol_enthalpy,total_enthalpy
    import ThermoState: mass_gibbs, mol_gibbs, total_gibbs
    import ThermoState: mass_helmholtz, mol_helmholtz, total_helmholtz
    import ThermoState: mass_internal_energy, mol_internal_energy, total_internal_energy
    import ThermoState: mass_entropy, mol_entropy, total_entropy
    import ThermoState: mass_fraction, mol_fraction
    import ThermoState: mass_number, mol_number
    import ThermoState: ThermoModel
    import ThermoState: mol_cp,mass_cp,mol_cv,mol_cp,sound_speed


    struct IF97Region{T} <: ThermoModel end
    
    const IF97_R =  0.461526
    const IF97_Tc = 647.096      #K          Critical temperature of water
    const IF97_Pc = 22.064       #Pa         Critical pressure of water
    const IF97_ρc = 322.         #kg/m3      Critical density of water
    const IF97_T3 = 273.16       #K          Triple point temperature of water
    const IF97_P3 = 611.657E-6   #MPa        Triple point pressure of water
    const IF97_Mr = 18.01528     #kg/kmol    Molecular weight of water
    const IF97_Tc_inv = 1.0/647.096
    const IF97_ρc_inv = 1.0/322.     
    struct IndustrialWater <: ThermoModel end
    
    temperature(::IndustrialWater,st::CriticalPoint,unit=u"K")  = convert_unit(u"K",unit,647.096)
    pressure(::IndustrialWater,st::CriticalPoint,unit=u"Pa")  = convert_unit(u"Pa",unit,2.2064e7)
    mol_volume(::IndustrialWater,st::CriticalPoint,unit=u"m^3/mol")  = convert_unit(u"m^3/mol",unit,5.594803726708074e-5)
    mol_density(::IndustrialWater,st::CriticalPoint,unit=u"mol/m^3")  = convert_unit(u"mol/m^3",unit,17873.72799560906)
    mass_density(::IndustrialWater,st::CriticalPoint,unit=u"kg/m^3")  = convert_unit(u"kg/m^3",unit,322.0)
    mass_volume(::IndustrialWater,st::CriticalPoint,unit=u"m^3/kg")  = convert_unit(u"m^3/kg",unit,0.003105590062111801)
    temperature(::IndustrialWater,st::TriplePoint,unit=u"K")  = convert_unit(u"K",unit,273.16)
    temperature(::IndustrialWater,st::NormalBoilingPoint,unit=u"K")  = convert_unit(u"K",unit,373.15)
    pressure(::IndustrialWater,st::TriplePoint,unit=u"Pa")  = convert_unit(u"Pa",unit,611.657)
    molecular_weight(::IndustrialWater) = 18.01528
    
    include("region1.jl")
    include("region2.jl")
    include("region3.jl")
    include("region4.jl")
    include("region5.jl")
    include("region_id.jl")
    include("industrialwater.jl")
    include("sat.jl")

    export IndustrialWater
    export mol_cp,mass_cp,mol_cv,mol_cp,sound_speed
end