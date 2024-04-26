"""
Call artifacts from Artifacts.toml
"""
module ArtifactCalling

using LazyArtifacts

export AIDA_ice_nucleation

"""
    AIDA_ice_nucleation(data_file_name)

 - `data_file_name` - name of the data file on Caltech box.

Returns the filepath of the data file in Caltech box.
"""
function AIDA_ice_nucleation(data_file_name::String)
    return joinpath(artifact"AIDA_ice_nucleation", data_file_name)
end

end
