const IF97_R =  0.461526
const IF97_Tc = 647.096      #K          Critical temperature of water
const IF97_Pc = 22.064       #Pa         Critical pressure of water
const IF97_ρc = 322.         #kg/m3      Critical density of water
const IF97_T3 = 273.16       #K          Triple point temperature of water
const IF97_P3 = 611.657E-6   #MPa        Triple point pressure of water
const IF97_Mr = 18.01528     #kg/kmol    Molecular weight of water

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



#SinglePT
    for op in [:helmholtz, :gibbs, :internal_energy, :enthalpy,:cp,:cv,:volume,:entropy]
        mol_op_impl = Symbol(:mol_,op,:_impl)
        mass_op_impl = Symbol(:mass_,op,:_impl)
        total_op_impl = Symbol(:total_,op,:_impl)
        mol_op = Symbol(:mol_,op)
        mass_op = Symbol(:mass_,op)
        total_op = Symbol(:total_,op)
        if op == :volume
            _unit = u"m^3/kg"
            mol_unit = u"m^3/mol"
            total_unit = u"m^3"

        elseif op in (:cv,:cp,:entropy)
            _unit = u"J/(kg*K)"
            mol_unit = u"J/(mol*K)"
            total_unit = u"J/(K)"
        else
            _unit = u"J/(kg)"
            mol_unit = u"J/(mol)"
            total_unit = u"J"
        end
     
        @eval begin
            function $mass_op_impl(mt::SinglePT,model::IndustrialWater,p,t)
                id = region_id(mt,model,p,t)
                return $mass_op_impl(mt,IF97Region{id}(),p,t)
            end


            function $mass_op(model::IndustrialWater,st::ThermodynamicState,unit=$_unit)
                return $mass_op(state_type(st),model,st,unit)
            end

            function $mass_op(mt::SinglePT,model::IndustrialWater,st::ThermodynamicState,unit)
                p = pressure(FromState(),st)
                t = temperature(FromState(),st)
                res = $mass_op_impl(mt,model,p,t)
                return convert_unit($_unit,unit,res)
            end

            function $mol_op(model::IndustrialWater,st::ThermodynamicState,unit=$mol_unit)
                prod = molar_mass(FromState(),st,u"kg/mol",molecular_weight(model))
                res =  $mass_op(model,st)*prod
                return convert_unit($mol_unit,unit,res)
            end 
        end

        if !(op in (:cv,:cp))
            
            @eval begin
                function $total_op(model::IndustrialWater,st::ThermodynamicState,unit=$total_unit)
                    prod = mass(FromState(),st,u"kg",molecular_weight(model))
                    res =  $mass_op(mt,model,st,unit)*prod
                    return convert_unit($total_unit,unit,res)
                end
            end
        end

    if op != :enthalpy
        @eval begin
            function $mass_op_impl(mt::SinglePH,model::IndustrialWater,p,h)
                id = region_id(mt,model,p,h)
                return $mass_op_impl(mt,IF97Region{id}(),p,h)
            end

            function $mass_op(mt::SinglePH,model::IndustrialWater,st::ThermodynamicState,unit)
                p = pressure(FromState(),st)
                h = mass_enthalpy(FromState(),st)
                res = $mass_op_impl(mt,model,p,h)
                return convert_unit($_unit,unit,res)
            end
        end   
    end
end

function sound_speed(model::IndustrialWater,st::ThermodynamicState,unit=u"m/s")
    return sound_speed(state_type(st),model,st,unit)
end

function sound_speed(mt::SinglePT,model::IndustrialWater,st::ThermodynamicState,unit)
    p = pressure(FromState(),st)
    t = temperature(FromState(),st)
    id = region_id(mt,model,p,t)
    res = sound_speed_impl(mt,IF97Region{id}(),p,t)
    return convert_unit(u"m/s",unit,res)
end

function sound_speed(mt::SinglePH,model::IndustrialWater,st::ThermodynamicState,unit)
    p = pressure(FromState(),st)
    h = mass_enthalpy(FromState(),st)
    id = region_id(mt,model,p,t)
    res = sound_speed_impl(mt,IF97Region{id}(),p,h)
    return convert_unit(u"m/s",unit,res)
end

    #==
for op in (:helmholtz, :gibbs, :internal_energy, :enthalpy)
    mol_op_impl = Symbol(:mol_,op,:_impl)
    mass_op_impl = Symbol(:mass_,op,:_impl)
    total_op_impl = Symbol(:total_,op,:_impl)
    mol_op = Symbol(:mol_,op)
    mass_op = Symbol(:mass_,op)
    total_op = Symbol(:total_,op)
    for spec in [:SinglePT,:SinglePH,:SinglePS]
    @eval begin 
        function $mass_op_impl(mt,model::HelmholtzModel,p,t)
            id = region_id(mt,model,p,t)
            return $mass_op_impl()
        end 

        function ThermoState.$mol_op(model::HelmholtzModel,st::ThermodynamicState,unit=u"J/mol")
            return $mol_op(state_type(st),model,st,unit)
        end    

        function $mol_op(mt::SingleVT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/mol")
            mw = molecular_weight(model)
            v = mol_volume(FromState(),st,u"m^3/mol",mw)
            t = temperature(FromState(),st)
            val = $mol_op_impl(mt,model,v,t)
            return convert_unit(u"J/mol",unit,val)
        end

    end
end

==#
