module Configs

export SpaceCfg, DynamicCfg, IntCfg

@kwdef mutable struct SpaceCfg
    length::Float64
    height::Float64
end

@kwdef struct DynamicCfg 
    ko::Float64
    ro::Float64
    ra::Float64
end

@kwdef struct IntCfg 
    dt::Float64
end

end