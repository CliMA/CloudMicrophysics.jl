# Run every file `runtests.jl` includes, each in its own testset, and report per file.
#
# `runtests.jl` wraps its includes in ONE testset, so a file that throws at load propagates out
# and every include below it is skipped. The suite still prints a plausible tree and total with
# no sign that a third of it never ran, which is how two separate defects hid for a whole
# rebuild. This driver isolates each file so one bad file costs one file, and enumerates every
# failure in a single run instead of one per iteration.
using Test
import CloudMicrophysics
const TESTDIR = @__DIR__

files = String[]
for line in eachline(joinpath(TESTDIR, "runtests.jl"))
    m = match(r"include\(\"([^\"]+)\"\)", line)
    m === nothing || push!(files, m.captures[1])
end
println("=== ", length(files), " files listed in runtests.jl")
flush(stdout)

results = Tuple{String, String}[]
for f in files
    t0 = time()
    try
        @testset "$f" begin
            include(joinpath(TESTDIR, f))
        end
        push!(results, (f, "OK"))
    catch e
        msg = first(sprint(showerror, e), 220)
        push!(results, (f, "THREW: " * replace(msg, "\n" => " | ")))
    end
    println("---- ", rpad(f, 44), " ", round(time() - t0, digits = 1), " s  ", results[end][2][1:min(end, 90)])
    flush(stdout)
end

println("\n=== SUMMARY")
for (f, r) in results
    println("  ", rpad(f, 44), " ", r === "OK" ? "OK" : r)
end
println("=== ", count(r -> r[2] == "OK", results), " of ", length(results), " files completed without throwing")
