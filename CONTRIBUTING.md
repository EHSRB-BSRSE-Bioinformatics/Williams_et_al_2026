# Contributing

This repository is a publication companion and prioritizes traceability to the
reported analysis over feature development. Corrections that improve
reproducibility, documentation, or error handling are welcome.

Before proposing a change:

1. Create a focused branch from the current default branch.
2. Run `pixi run check`.
3. If plotting code or compact summaries changed, run `pixi run figures` and
   inspect all PDF and PNG outputs.
4. Update `docs/data_manifest.tsv` when an included data file changes.
5. Explain whether the change affects numerical results or only packaging.

Do not submit credentials, private URLs, confidential data, personally
identifying information, or unsanitized logs. Large trial-level objects belong
in the archival deposit, not Git.

