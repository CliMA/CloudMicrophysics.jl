# The following comments are instructions to use ClimaArtifacts within CloudMicrophysics.jl
# Directions will differ if you expect to use the artifact in MPI runs (i.e. as part
# of a function in CloudMicrophysics that is called by ClimaAtmos). For more info, please
# visit the ClimaArtifacts.jl repo on GitHub!

# To use in VS Code, first activate the AIDA_data_artifact environment.
# Then, copy and paste the following script to the Julia REPL.
# To use in terminal, activate Julia at --project=papers/ice_nucleation_2024/AIDA_data_artifact.
# Then, copy and paste the following script in your now opened Julia REPL.

# If you do not already have ClimaArtifactsHelper, go into the Julia REPL and type:
# Pkg.develop(url="https://github.com/CliMA/ClimaArtifacts.git", subdir="ClimaArtifactsHelper.jl")
using ClimaArtifactsHelper

# path of the folder in the artifact in which the actual data files are stored
data_folder = "AIDA_ice_nuc_data"

# naming this new set of data (AKA the artifact) for CloudMicrophysics.jl's Artifact.toml
data_set_name = "AIDA_ice_nucleation"

create_artifact_guided(data_folder; artifact_name = data_set_name)

# You will be guided in the terminal/REPL to upload the files on the Caltech box.
# Make sure to copy the downloadable link, *not* the normal share link.

# After getting an updated "OutputArtifacts.toml" file in this folder, 
# Add it to the list of artifacts in CloudMicrophysics.jl's Artifact.toml file.
# Add `lazy = true` under the `git-tree-sha1` line in Artifact.toml i fyou would like
# the artifact to be downloaded lazily.

# To use data, make sure your script has "using LazyArtifacts" at the top.
# The following command will return the artifact folder: artifact"data_set_name"
# To choose a specific data file, use joinpath to append the name of the file as
# uploaded to the artifacts folder (in this case, "AIDA_ice_nuc_data").
