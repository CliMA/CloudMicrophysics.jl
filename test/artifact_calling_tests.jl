import Test as TT

import CloudMicrophysics as CM
import CloudMicrophysics.ArtifactCalling as AFC


@info "Artifact Calling Tests"

function test_artifact_calling()

    TT.@testset "AIDA ice nucleating experiments" begin

        data_file_names = ["in05_17_aida.edf"]

        # File exists and is callable
        for data_file_name in data_file_names
            file_path = AFC.AIDA_ice_nucleation(data_file_name)
            TT.@test isfile(file_path)
        end

    end
end

test_artifact_calling()
