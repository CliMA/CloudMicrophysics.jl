export MiddleEasternDust

"""
    MiddleEasternDust{FT}

Parameters for Middle Eastern Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
CP.@params struct MiddleEasternDust{FT} <: AerosolType
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end
