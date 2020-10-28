function tsat_check(t)
    if !(IF97_T3 ≤ T ≤ IF97_Tc)
        return throw(DomainError(T, "Temperature not between tripple and critical points"))
    else
        return nothing
    end
end

function is_liquid(x::Symbol)::Bool
    x in (:liquid,:l,:L,:LIQUID)
end

function is_gas(x::Symbol)::Bool
    x in (:gas,:g,:G,:GAS,:vapor,:VAPOR,:V,:v,:steam,:STEAM)
end

const satvl_b = [   1.992_740_64,
1.099_653_42,
-0.510_839_303,
-1.754_934_79,
-45.517_035_2,
-6.746_944_50E5]

satvv_c = [  -2.031_502_40,
-2.683_029_40,
-5.386_264_92,
-17.299_160_5,
-44.758_658_1,
-63.920_106_3
]
function massrho_l_impl(::SingleSatT,model::IndustrialWater,t::T) where T
    t =normalize_units(T)
    tsat_check(t)
    θ = T*IF97_Tc_inv
    τ = 1 - θ
    b = satvl_b
    return ρc*(T(1) + b[1]*τ^T(1/3) + b[2]*τ^T(2/3) + b[3]*τ^T(5/3) + b[4]*τ^T(16/3)
           + b[5]*τ^T(43/3) + b[6]*τ^T(110/3))
end

function massrho_v_impl(::SingleSatT,model::IndustrialWater,t::T) where T
    t =normalize_units(T)
    tsat_check(t)
    c = satvv_c
    θ = T*IF97_Tc_inv
    τ = 1 - θ

    return ρc*exp(c[1]*τ^T(2/6) + c[2]*τ^T(4/6) + c[3]*τ^T(8/6) + c[4]*τ^T(18/6)
    + c[5]*τ^T(37/6) + c[6]*τ^T(71/6))
end

function mass_density(mt::SingleSatT,model::IndustrialWater,st::ThermodynamicState,unit="kg/m3")
    sat = phase(FromState(),st)
    t = temperature(FromState(),t)
    if is_liquid(sat)
        res = massrho_l_impl(mt,model,t)
        return convert_unit("kg/m3",unit,res)
    elseif is_gas(sat)
        res =  massrho_v_impl(mt,model,t)
        return convert_unit("kg/m3",unit,res)
    else
        throw(error("unspecified required phase. Stablish phase = :gas or phase = :liquid"))
    end
end


function mass_volume(mt::SingleSatT,model::IndustrialWater,st::ThermodynamicState,unit)
    rho = mass_density(mt,model,st)
    res = one(res)/res
    return convert_unit(u"m^3/kg",unit,res)
end

#mol_volume does the following dispatch, in the case of SingleSatT

#mol_volume(model,st) - macro
#mol_volume(mt,model,st) - macro, multiplies mass volume with molar mass
#mass_volume(mt,model,st) - normally to mass_volume_impl, here to mass_density
#mass_density(mt,model,st)

#to add total and mol units, it suffices just to add the mass unit.


"""
    Psat2

    Do not call this function, but rather use Psat().
    Returns the saturated vapour pressure according to IAPWS SR1-86(1992)
    with T in K and P in MPa.
    This replaces the Region 4 equations for Psat and is used to calcualate
    the saturated enthalpy and entropy to match the values in IAPWS SR1-86(1992)
    It is not used in Tsat(), as the difference in vapour pressure is in the
    8th decimal at the boiling point.
"""

const dpsat_a = [ -7.85951783000,
1.84408259000,
-11.78664970000,
22.68074110000,
-15.96187190000,
1.80122502000
]

function ∂Psat∂T(T)
    a = dpsat_a
    θ = T*IF97_Tc_inv 
    τ = 1 - θ
    ∂τ = - IF97_Tc_inv
    ∂T = one(T)
    τhalf = sqrt(τ)
    τ15 = τ*τhalf
    τ2 = τ*τ
    τ25 = τ2*τhalf
    τ3 = τ2*τ
    τ35=τ3*τhalf
    τ4=τ2*τ2
    τ65 = τ35*τ3
    τ75 = τ35*τ4

    x = IF97_Tc*(a[1]*τ + a[2]*τ15 + a[3]*τ3 + a[4]*τ35 + a[5]*τ4 + a[6]*τ75)
    ∂x = IF97_Tc*∂τ*(a[1] + 1.5*a[2]*τhalf + 3.0*a[3]*τ2 + 3.5*a[4]*τ25 + 4.0*a[5]*τ3 + 7.5*a[6]*τ65)
    xT = x/T
    #res =  Pc*exp(xT)
    ∂res = IF97_Pc*exp(xT)*(∂x - xT*∂T)/T
end



