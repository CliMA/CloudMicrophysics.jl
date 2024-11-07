import Test as TT
import CloudMicrophysics.Parameters as CMP

@info "Test broadcasting over types"

function test_common_types_broadcasts(FT)
    for some_type in (
        CMP.ParametersType{FT},
        CMP.AerosolType{FT},
        CMP.AerosolDistributionType,
        CMP.CloudCondensateType{FT},
        CMP.PrecipitationType{FT},
        CMP.TerminalVelocityType{FT},
        CMP.Precipitation2MType{FT},
        CMP.ArizonaTestDust{FT},
        CMP.DesertDust{FT},
        CMP.Kaolinite{FT},
        CMP.Illite{FT},
        CMP.Seasalt{FT},
        CMP.Sulfate{FT},
        CMP.H2SO4SolutionParameters{FT},
        CMP.H2S04NucleationParameters{FT},
        CMP.OrganicNucleationParameters{FT},
        CMP.MixedNucleationParameters{FT},
        CMP.AerosolActivationParameters{FT},
        CMP.Koop2000{FT},
        CMP.Mohler2006{FT},
        CMP.AirProperties{FT},
        CMP.WaterProperties{FT},
        CMP.Blk1MVelTypeRain{FT},
        CMP.Blk1MVelTypeSnow{FT},
        CMP.SB2006VelType{FT},
        CMP.Chen2022VelTypeSmallIce{FT},
        CMP.Chen2022VelTypeLargeIce{FT},
        CMP.Chen2022VelTypeRain{FT},
        CMP.Parameters0M{FT},
        CMP.ParticlePDFSnow{FT},
        CMP.ParticlePDFIceRain{FT},
        CMP.ParticleMass{FT},
        CMP.ParticleArea{FT},
        CMP.Ventilation{FT},
        CMP.Acnv1M{FT},
        CMP.CloudLiquid{FT},
        CMP.CollisionEff{FT},
        CMP.AcnvKK2000{FT},
        CMP.AccrKK2000{FT},
        CMP.AcnvB1994{FT},
        CMP.AccrB1994{FT},
        CMP.AcnvTC1980{FT},
        CMP.AccrTC1980{FT},
        CMP.LD2004{FT},
        CMP.VarTimescaleAcnv{FT},
        CMP.RainParticlePDF_SB2006{FT},
        CMP.CloudParticlePDF_SB2006{FT},
        CMP.AcnvSB2006{FT},
        CMP.AccrSB2006{FT},
        CMP.SelfColSB2006{FT},
        CMP.BreakupSB2006{FT},
        CMP.EvaporationSB2006{FT},
    )
        Base.materialize(Base.Broadcast.broadcasted(identity, some_type))
    end
end

test_common_types_broadcasts(Float32)
test_common_types_broadcasts(Float64)
