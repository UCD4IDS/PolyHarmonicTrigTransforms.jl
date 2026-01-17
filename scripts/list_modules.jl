using PolyHarmonicTrigTransforms

function walk(mod, prefix=string(mod))
    for nm in names(mod, all=true, imported=false)
        if isdefined(mod, nm)
            obj = getfield(mod, nm)
            if isa(obj, Module) && obj !== Main && obj !== Base
                println(prefix * "." * string(nm))
                walk(obj, prefix * "." * string(nm))
            end
        end
    end
end

walk(PolyHarmonicTrigTransforms)
