using Downloads
using SHA
using TOML

function file_sha256(path)
    bytes2hex(open(SHA.sha256, path))
end

function main(arguments = ARGS)
    output_directory =
        isempty(arguments) ? joinpath(@__DIR__, "data") : abspath(first(arguments))
    mkpath(output_directory)
    manifest = TOML.parsefile(joinpath(@__DIR__, "public_structures.toml"))
    for record in values(manifest["structures"])
        target = joinpath(output_directory, record["filename"])
        if !isfile(target)
            Downloads.download(record["url"], target)
        end
        observed = file_sha256(target)
        observed == record["sha256"] || error(
            "checksum mismatch for $(record["filename"]): expected $(record["sha256"]), observed $observed",
        )
        println("Verified $(record["accession"]): $target")
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
