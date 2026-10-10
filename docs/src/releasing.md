# Releasing and registration

Releases follow semantic versioning and are validated on Julia 1.10 and the
current stable Julia release.

## Release checklist

1. Update `Project.toml` and `CHANGELOG.md`.
2. Run the offline test suite and controlled examples.
3. Run the checksum-pinned public-structure validation.
4. Confirm that documentation builds without warnings.
5. Register the release with Julia Registrator.
6. Confirm the tag, GitHub release, documentation deployment, and archival DOI.

New Julia General registrations are triggered from a commit in this repository
using `@JuliaRegistrator register`. TagBot is configured for future registered
versions. The v0.1.0 tag and release were created manually before initial
registry submission.

The scheduled public-data workflow periodically re-runs the RCSB 1J7T and 4V7U
checks without committing downloaded structures.
