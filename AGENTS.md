# scHeatmap Agent Guidelines

This file defines how coding agents should work in this repository. It applies
to the entire project.

## Project goal

`scHeatmap` is an R package that provides simple, reproducible single-cell
heatmaps from Seurat objects while retaining the flexibility of
ComplexHeatmap. Keep the common workflow easy, but preserve support for the
specialized CR2020 literature-reproduction workflow.

## Communication

- Communicate with the maintainer in concise Chinese unless asked otherwise.
- Lead with the result or proposed behavior, especially when discussing an API
  change.
- For changes that materially alter plotting semantics, explain the proposed
  design and trade-offs before implementation unless the user already approved
  a concrete design.
- Agents may commit completed implementation work automatically when the
  conditions under [Automatic commit policy](#automatic-commit-policy) are
  satisfied. Never push without explicit user authorization. Push directly to
  `main` unless the user asks for another branch or pull request.

## Automatic commit policy

An agent may create a local Git commit without a separate request only when all
of the following are true:

- The user requested an implementation, modification, or bug fix rather than
  analysis, diagnosis, review, or a proposed plan.
- The requested task is complete and no material design choice is waiting for
  user confirmation.
- Relevant tests and checks pass, or a documentation-only change has passed the
  appropriate formatting and example validation.
- The diff has a clear scope and contains no unrelated user changes.
- The commit contains no temporary files, credentials, sensitive information,
  unexpected large files, or accidental generated output.
- The work is not in the middle of iterative debugging or waiting for the user
  to test and refine the result.
- The user has not said to avoid committing, to show a plan first, or to wait
  for approval.

Do not commit automatically when any of these conditions is uncertain. Leave
the changes in the working tree and report what remains to be decided or
validated.

Automatic commit permission does not include `git push`. Pushing changes to a
remote always requires an explicit user request. After authorization, push
directly to `main` unless another branch or pull-request workflow was requested.

## Repository layout

- `R/`: exported functions and internal package helpers.
- `tests/testthat/`: automated tests for public behavior and regressions.
- `man/`: roxygen-generated documentation; do not edit generated `.Rd` files
  independently of their roxygen source.
- `development/demo.R`: complete executable development examples, including
  PBMC3K and the CR2020 reproduction.
- `development/legacy_utils.R`: reference implementation and legacy palettes;
  migrate useful behavior into `R/` rather than sourcing this file at runtime.
- `inst/examples/scHeatmap-gallery.md`: compact gallery of supported workflows.
- `data-raw/`: development data excluded from the installed package.
- `man/figures/`: selected, versioned images that README files may reference.

## Development environment

The maintainer uses the conda environment `scptools`. Prefer commands such as:

```bash
conda run -n scptools Rscript -e 'devtools::document()'
conda run -n scptools Rscript -e 'devtools::test()'
conda run -n scptools Rscript -e 'devtools::check(document = FALSE, manual = FALSE, cran = FALSE)'
```

During interactive development, `devtools::load_all(".")` is preferred over
reinstalling after every change. Do not claim the package is installed merely
because it has been loaded with `load_all()`.

## Implementation rules

- Keep `sc_heatmap()` returning a standard `ComplexHeatmap::Heatmap` object.
- Preserve compatibility with documented legacy arguments unless removal is
  explicitly approved.
- Put reusable plotting, palette, sizing, and validation helpers inside the
  package, not only in `development/` scripts.
- Prefer explicit argument names that describe what is displayed, for example
  `show.feature.names`, `show.column.names`, and `show.group.names`.
- Preserve user-supplied ComplexHeatmap arguments passed through `...`; an
  explicit native argument should override a package default when documented.
- Validate metadata fields, features, group mappings, and custom colors early
  and return actionable error messages.
- Keep random operations reproducible through a user-visible seed without
  unnecessarily changing the caller's random state.

## Expression and aggregation semantics

- Default expression comes from the active assay's `data` layer unless the
  user selects another assay or layer.
- In `mode = "average"`, calculate the arithmetic mean for each feature across
  cells in every `aggregate.by` group, then apply the requested scaling and
  clipping to the aggregated matrix.
- The default average is cell-weighted mean log-normalized expression; do not
  describe it as pseudobulk.
- Do not silently change averaging to raw-count summation, sample-aware
  averaging, or pseudobulk. Add an explicit mode or argument, documentation,
  and tests if such behavior is requested.
- Keep scaling order and clipping order documented and covered by tests.

## Plotting behavior to preserve

- `label.features` labels only selected genes but keeps all requested features
  in the matrix. Linked labels must use `ComplexHeatmap::anno_mark()` so lines
  point to their heatmap rows and survive `save_sc_heatmap()`.
- Linked gene labels are italic by default; users can select plain text or
  provide `label.gp` for advanced styling.
- `cell_type` annotations use the migrated `divergentcolor` palette by default,
  while explicit `annotation.colors` always takes precedence.
- Heatmap slices have black borders by default. Users can disable or override
  the border color, and an explicit ComplexHeatmap `border` argument takes
  precedence.
- Display switches for feature labels, column labels, split titles, and top
  annotations must operate independently.

## Documentation and examples

- Public API changes require synchronized roxygen comments, regenerated `man/`
  documentation, and tests.
- Keep `README.md` and `README.zh-CN.md` synchronized. Both must retain the
  language-switch badges.
- README examples should be runnable. Use the concise PBMC3K example for normal
  usage and retain the CR2020 section because the README hero image is produced
  by that workflow.
- Put lengthy marker discovery and export workflows in `development/demo.R`;
  keep README examples focused and link to the complete script.
- Update `inst/examples/scHeatmap-gallery.md` when a new major visualization
  pattern is introduced.
- Only reference images that are intentionally stored under `man/figures/`.

## Tests and validation

- Add or update `testthat` tests for every public argument, bug fix, default
  change, and important interaction between display options.
- At minimum, run `devtools::test()` after code changes.
- Run `devtools::document()` whenever roxygen or exports change, then verify the
  resulting `NAMESPACE` and `.Rd` changes.
- Before publishing a substantial change, run `devtools::check()` and report
  errors, warnings, and notes accurately.
- Documentation-only changes require `git diff --check`; executable examples
  should also be validated locally when the required data is available.
- Empty `Rplots.pdf` files are test artifacts, not results. Never commit them.

## Data and generated files

- The CR2020 RDS is intentionally versioned with Git LFS. Preserve its LFS
  tracking and do not move or re-add it as a normal Git blob.
- Reorganizing local directories does not require uploading unchanged LFS
  content again unless the tracked path or file content changes.
- Files generated under `development/figures/` are local outputs and are
  ignored. Promote a selected image to `man/figures/` only when it is intended
  for package documentation.
- Do not delete user data, RDS files, or generated results without explicit
  confirmation.

## Git hygiene

- Inspect `git status`, the diff, and `git diff --check` before staging.
- Preserve unrelated user changes and stage only files belonging to the current
  task.
- Use concise commits that describe the outcome.
- Never use destructive commands such as `git reset --hard` or discard user
  work.
- After pushing, verify that the local branch matches `origin/main` and report
  the commit hash and validation performed.
