function prepare_results(base::AbstractString; overwrite::Bool = false)
    target = abspath(base)
    if ispath(target) && !overwrite
        throw(
            ArgumentError(
                "output path already exists: $target; choose a new path or pass overwrite=true explicitly",
            ),
        )
    end
    mkpath(target)
    for directory in ("tables", "figures", "structures", "logs")
        mkpath(joinpath(target, directory))
    end
    target
end

function _file_checksum(path)
    isfile(path) || return nothing
    bytes2hex(open(SHA.sha256, path))
end

const _DIRECT_DEPENDENCY_MODULES = Dict(
    "BioStructures" => BioStructures,
    "CSV" => CSV,
    "DataFrames" => DataFrames,
    "JSON3" => JSON3,
    "NearestNeighbors" => NearestNeighbors,
    "StaticArrays" => StaticArrays,
    "StatsBase" => StatsBase,
)

function _dependency_versions()
    Dict(
        name => string(Base.pkgversion(package)) for
        (name, package) in _DIRECT_DEPENDENCY_MODULES
    )
end

function write_metadata(
    output_dir::AbstractString;
    inputs = String[],
    parameters = Dict(),
    random_seed = nothing,
    command = join(Base.ARGS, " "),
)
    input_records =
        [Dict("path" => abspath(path), "sha256" => _file_checksum(path)) for path in inputs]
    metadata = Dict{String,Any}(
        "created_at" => string(Dates.now(Dates.UTC)),
        "julia_version" => string(VERSION),
        "package_version" => string(Base.pkgversion(@__MODULE__)),
        "dependency_versions" => _dependency_versions(),
        "inputs" => input_records,
        "parameters" => parameters,
        "random_seed" => random_seed,
        "command" => command,
        "coordinate_unit" => "angstrom",
    )
    open(joinpath(output_dir, "metadata.json"), "w") do io
        JSON3.pretty(io, metadata)
    end
    metadata
end

function write_summary(output_dir::AbstractString, summary)
    open(joinpath(output_dir, "summary.json"), "w") do io
        JSON3.pretty(io, summary)
    end
    joinpath(output_dir, "summary.json")
end

function write_table(output_dir::AbstractString, filename::AbstractString, table)
    endswith(lowercase(filename), ".csv") ||
        throw(ArgumentError("table filename must end in .csv"))
    path = joinpath(output_dir, "tables", filename)
    CSV.write(path, table)
    path
end
