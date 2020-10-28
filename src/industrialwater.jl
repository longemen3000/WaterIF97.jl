



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

            #dispatch basic mass impl to the correct region, P,T
            function $mass_op_impl(mt::SinglePT,model::IndustrialWater,p,t)
                id = region_id(mt,model,p,t)
                return $mass_op_impl(mt,IF97Region{id}(),p,t)
            end


            function $mass_op(model::IndustrialWater,st::ThermodynamicState,unit=$_unit)
                return $mass_op(state_type(st),model,st,unit)
            end
            # P T impl
            function $mass_op(mt::SinglePT,model::IndustrialWater,st::ThermodynamicState,unit)
                p = pressure(FromState(),st)
                t = temperature(FromState(),st)
                res = $mass_op_impl(mt,model,p,t)
                return convert_unit($_unit,unit,res)
            end

            #mol op
            function $mol_op(model::IndustrialWater,st::ThermodynamicState,unit=$mol_unit)
                prod = molar_mass(FromState(),st,u"kg/mol",molecular_weight(model))
                res =  $mass_op(state_type(st),model,st,$_unit)*prod
                return convert_unit($mol_unit,unit,res)
            end 
        end

        if !(op in (:cv,:cp))
            #total ops
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
            #if not enthalpy, define PH impl
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

    if op != :entropy
        @eval begin
            #if not entropy, define PS impl
            function $mass_op_impl(mt::SinglePS,model::IndustrialWater,p,s)
                id = region_id(mt,model,p,h)
                t = temperature_impl(mt,IF97Region{id}(),p,s) #transform to SinglePT
                _mt = QuickStates.pt()
                return $mass_op_impl(_mt,model,p,t)
            end

            function $mass_op(mt::SinglePS,model::IndustrialWater,st::ThermodynamicState,unit)
                p = pressure(FromState(),st)
                s = mass_entropy(FromState(),st)
                res = $mass_op_impl(mt,model,p,s)
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
function mass_density(model::IndustrialWater,st::ThermodynamicState,unit=u"kg/(m^3)")
    return mass_density(state_type(st),model,st,unit)
end

function mol_density(model::IndustrialWater,st::ThermodynamicState,unit=u"mol/(m^3)")
    return mol_density(state_type(st),model,st,unit)
end

function mass_density(mt,model::IndustrialWater,st::ThermodynamicState,unit)
    res = mass_volume(model,st)
    rho = one(res)/res
    return convert_unit(u"kg/(m^3)",unit,rho)
end

function mol_density(mt,model::IndustrialWater,st::ThermodynamicState,unit)
    res = mol_volume(model,st)
    rho = one(res)/res
    return convert_unit(u"mol/(m^3)",unit,rho)
end

