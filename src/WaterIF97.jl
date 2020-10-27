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
    #import ThermoState: mass_cp,mass_cv,mol_cp,mol_cv
    #import ThermoState: sound_speed

    struct IF97Region{T} <: ThermoModel end

    include("region1.jl")
    include("region2.jl")
    include("region3.jl")
    include("region4.jl")
    include("region5.jl")
    include("industrialwater.jl")
    include("region_id.jl")

    export IndustrialWater
    export mol_cp,mass_cp,mol_cv,mol_cp,sound_speed
end