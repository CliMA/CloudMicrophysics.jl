# Artifact Calling

## Using Artifacts
Calling the artifact needs to be done in the `CloudMicrophysics` source code so that the 
    `Artifact.toml` file can be seen. The file containing the list of artifacts is
    `ArtifactCalling.jl`. To call an artifact in a certain file outside of the `src` folder,
    simply call the calling function just as you would other functions in the source code.
    For example:
    
```
import CloudMicrophysics as CM

data_file_name = "in05_17_aida.edf"
CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name)
```
This will automatically download the `in05_17_aida.edf` file from the Caltech box,
    if not already available, and output the filepath.

## More information
For more information about artifacts, visit the `ClimaArtifacts` [repo on GitHub](https://github.com/CliMA/ClimaArtifacts).
