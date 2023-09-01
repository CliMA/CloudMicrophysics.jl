import Test as TT
import CloudMicrophysics.CommonTypes as CMT

@info "Test broadcasting over types"

function test_common_types_broadcasts()
    for some_type in (
        CMT.AbstractCloudType,
        CMT.AbstractPrecipType,
        CMT.LiquidType,
        CMT.IceType,
        CMT.RainType,
        CMT.SnowType,
        CMT.AbstractAerosolDistribution,
        CMT.Abstract2MPrecipType,
        CMT.KK2000Type,
        CMT.B1994Type,
        CMT.TC1980Type,
        CMT.LD2004Type,
        CMT.SB2006Type,
        CMT.AbstractTerminalVelocityType,
        CMT.Blk1MVelType,
        CMT.SB2006VelType,
        CMT.Chen2022Type,
        CMT.AbstractAerosolType,
        CMT.ArizonaTestDustType,
        CMT.DesertDustType,
        CMT.KaoliniteType,
        CMT.IlliteType,
    )
        Base.materialize(Base.Broadcast.broadcasted(identity, some_type))
    end
end

test_common_types_broadcasts()
