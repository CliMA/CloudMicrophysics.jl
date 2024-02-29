# Developer's Guide

This doc is meant to guide any beginner developers with organization
  and common surprises associated with contributing to `CloudMicrophysics.jl`.

## GitHub workflow

We work on branches of the main repository, rather than personal forks.
Each branch name should ideally start with your initials, followed by short
  name relevant to what will be added in that branch.
It is considered a good practice to open an issue that describes the planned work
  before starting the implementation.
This allows you to get feedback from other developers and advertise the planned work to the group.
When the implemented code is ready for review, create a pull request (PR) on GitHub
  and tag the issue it is addressing.
In the PR you can describe which parts of the relevant issue are solved and what needs to
  be done to reach that goal.
After creating a PR on GitHub, you will find that every following commit & push you make
  will go through continuous integration (CI) checks:

- The, `ci` runs the tests on Ubuntu, Windows and OSX machines using GitHub Actions.
  The tests include some unit tests and very simple performance tests.
  It may happen that the tests pass on your local machine,
  but fail on one of the machines provided in the cloud for the CI.
  This is especially true for the performance tests,
  as the execution times may vary a lot depending on the machine.
  The CI is the source of truth over local tests, and all tests must pass in the CI before code can be merged.
  See the [Tests](https://clima.github.io/CloudMicrophysics.jl/dev/DevelopersGuide/#Tests)
  section below for debugging tips.
- `CloudMicrophysics.jl` has a small set of tests that are run on the GPUs using `buildkite`.
  They are triggered automatically when trying to merge a PR.
  If the GPU tests fail, check first for use of any "out of place" functions
  (for example, `@warn` will not work on GPU).
- The `Documentation / docbuild` builds the documentation.
  Most common errors in the docbuild are related to missing docstrings,
  see the [Documentation](https://clima.github.io/CloudMicrophysics.jl/dev/DevelopersGuide/#Documentation) section for hints.
  This check may fail if the source code itself is not compiling,
  so please handle the compilation errors first.
  If the documentation build was successful, the `documenter/deploy`
  will display the documentation page based on the PR (click on the details).
- The `JuliaFormatter / format` ensures consistent formatting throughout the repository.
  If this check fails, you can click on details
  to see which files are not following the formatting rules.
  You can apply the formatter by running
  `julia --project=.dev .dev/climaformat.jl file_name` in the terminal.
  The formatter might break if the code does not compile,
  so make sure you handle the compilation errors first.
- The `codecov` shows what percentage of source code lines are exercised inside the tests.
  We strive to keep the test coverage high, but this test is not strictly required to pass before merging a PR.

Once the PR is opened, you can request reviewers to look over and approve your work.
To keep things tidy, we want PRs to have few (ideally only one) commit
  before merging to the main branch.
You can squash and rebase multiple commits into one
  in your favorite editor (VS Code, Vim etc).
In case of doubt, see Git tutorials on how to squash and rebase
  or reach out to other developers in our team.
The first couple of times it pays off to create a "backup branch" before starting your rebase,
  and then comparing afterwards if the rebased and backup branches are identical.
If the main branch has been updated after your branch was created,
you will need to rebase onto the main branch. To do this,
run git rebase origin/main and solve any conflicts manually.
This is done more easily if you only have one commit.
You can merge the PR if all the tests are passing
  and the PR branch is up to date with the main branch.

The easiest way to use the new additions in the `main` branch of `CloudMicrophysics.jl`
  from another package is to do a package release.
To do that, you need to change the package version in the `Project.toml`.
We do a patch release if the API did not change (for example `0.11.1 -> 0.11.2`)
  and a minor version release if it did (for example `0.11.1 -> 0.12.0`).
We use the `JuliaRegistrator` bot to register new versions in the `Julia` package ecosystem,
  by commenting `@JuliaRegistrator register` under the merged PR that changes the package version.
You can also use specific branches from the repository,
  if you are working with commits that were not yet merged into `main`.
See [Pkg.jl docs](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-registered-packages)
  for more details.

## New contributions

`CloudMicrophysics.jl` is a collection of point-wise functions grouped
  in different modules depending on which aspect of aerosol and cloud microphysics
  they address.
When adding a new function, decide if a new module is required, or add to the existing one.
If the added function will be used in other modules, export the function.
Avoid exporting functions that are only used within the module it is defined in
  or functions that have the same name as a function already in the API.
Avoid shortening words or phrases when naming new functions.
The name of the function should be self-explanatory yet brief.
All functions should have docstrings describing the API, as well as
  documentation focusing on the scientific aspects of what they do.
All functions should have their own unit, performance and GPU tests.
All free parameters should be stored in [ClimaParams](https://github.com/CliMA/ClimaParams.jl).
It is usually faster to prototype defining the free parameters locally,
  and move them to `ClimaParams` at a last step.

Other files that may require editing after you make a new function are:
 - `CloudMicrophysics.jl` (found in `src` folder) if you need to include a new source file,
 - `index.md` (found in `docs/src` folder) if you want to mention the new additions in the main documentation page,
 - `make.jl` (found in `docs` folder) if you are adding a new file to the documentation,
 - `API.md` (found in `docs/src` folder) if you are adding functions (see ``Documentation`` section for more details),
 - `runtests.jl` (found in `test` folder) if you are adding new test file for unit tests.

## Documentation

Each new addition to the library should be accompanied by a documentation page
  summarizing the derivation/assumptions and it's potential uses.
If possible, it's really appreciated to also add short code snippets that reproduce
  results from the literature.
Those code snippets are executed every time documentation is build and provide
  great examples on how different available parameterizations work.

Additionally, each function in the source code is required to have a docstring
  and should be added to the `API` documentation page.
The docstring formatting is pretty strict:
  (i) no empty lines between the docstring and the function/module,
  (ii) docstring starts and ends with a line consisting of three quotation marks,
  (iii) the second line has the function's name preceded by exactly 4 spaces.
Examples can be found in the source files.
Be sure to add the function to `API.md` in the docs folder.
Missing docstrings and functions in the API cause the documentation build to break.

To build and work the documentation locally, you can use [LiveServer.jl](https://github.com/tlienart/LiveServer.jl#serve-docs).
This will compile the documentation to a local server that is updated whenever you make changes.
Alternatively, you can save the documentation to a static webpage: `julia --project=docs docs/make.jl`.
The index page will be saved in `docs/build` folder.

## Tests

### Unit Tests

Unit tests aim to ensure that parameterizations make physical sense.
They can be found in the `test` folder under files named after corresponding source files.
If you create a new function, please also create a new test that checks it.
If creating a new file for unit tests, make sure you import `Test` and any other necessary libraries.
There is some boilerplate code needed to create the sets with free parameters
  based on the default `toml` file from [ClimaParams](https://github.com/CliMA/ClimaParams.jl).

Some possible tests include checking if the returned values agree with values
  in the literature, if something is smaller/greater at warmer/cooler
  temperatures, if assertion errors are returned when a function is used outside its
  valid range of parameters, or if a function is zero at certain input values.
In general, writing good tests is difficult and we are always on the lookout for new good candidates.
We strive to exercise all functions in some way in tests,
  so that at minimum we can catch changes in the API.

You can run the tests locally: `julia --project=test test/runtests.jl`.

### Performance Tests

Performance tests check the memory allocations (there should be none) and execution times
  of some of the functions.
They are found in the `test` folder under a single file named `performance_tests.jl`.

You can add performance tests for your new functions in the `benchmark_test(FT)` function.
Parameters are listed at the start of the function.
Next, add your performance tests under the comment with the file name which your
  new function is in.
This is done by calling `bench_press(function_name, (Parameter1, Parameter2), min_time)`.
The last argument is an estimate of the minimum run-time of the function.
It may take some trial and error to find a number
  that will satisfy the `ci` tests that are run on different operating systems
  and for different floating-point precision.
You can identify which specific performance test is failing on GitHub
  by clicking on details next to a failed check from one of the `ci / ci 1.8.1` checks.
Adjust the estimated minimum run-time appropriately.

### GPU Tests

Graphic Processing Unit (GPU) tests check that `CloudMicrophysics.jl` functions are able to run on GPU.
They are as simple as checking that a certain input returns a known value.
Right now, we do not test the whole library on the GPUs,
  so the support is limited.
GPU tests can be found in the `test` folder under a single file named `gpu_tests.jl`.
There are two things to add: a kernel for your function and the actual test.

Kernels are added at the top half of the file using `@kernel`.
The naming convention follows `test_<function name>_kernel!`.
Within the kernel, use `@inbounds` to place any outputs into an output array.

The test itself should be added in the `test_gpu(FT)` function.
The test starts with defining `data_length` and ends after an `@test` macro.
`data_length` corresponds to the number of outputs your function has.
Add an array for each input required by your function.
For example, if you want to test at temperature of 230K,
  you can add `T = ArrayType([FT(230)])`.
Add a comment of what you are testing and use `@test` to create your test.
The GPU tests are ran twice: for `Float64` and `Float32`.
Similar as with performance tests, some trial and error is needed
  to find good tolerances for both options.
