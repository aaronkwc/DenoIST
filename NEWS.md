## Version 1.0.1 (2026-05-24)

* Added `on_disk` mode to `denoist()` for large datasets that do not fit in memory.

## Version 0.99.4 (2026-03-17)

* Added `local_offset_distance_with_background_fast` function which uses `dbscan` for faster neighour finding.

* `denoist()` now defaults to fast neighbour finding via the `neighbour_mode` option.

* Fixed minor bug where background offset cannot be calculated because an entire gene gets filtered out because of low qv. This should not change existing usage as the issue only arises in extremely small toy datasets.

* `n_inits` can now be tuned in the `denoist()` function for speed.

## Version 0.99.3 (2026-02-24)

* Removed `print` from README.

* Fixed minor bug where background offset cannot be calculated because an entire gene gets filtered out because of low qv.

* `n_inits` can now be tuned in the `denoist()` function for speed.

## Version 0.99.2 (2026-02-05)

* Removed commented out chunks during development.

* Added `verbose` option to make progress messages off by default.

* Replaced `dplyr`, `tidyr`, and `tibble` dependency with base R `xtabs`.

* Specified functions imported from `flexmix` instead of importing the whole package.

* Updated vignette to use `BiocStyle`.

* Trimmed README content for a more pleasant reading experience.

* Vignette now uses an even smaller dataset.

## Version 0.99.1 (2026-01-19)

* `local_offset_distance_with_background` can now handle non-Xenium transcript
data with no QV column without data modification from user.

## Version 0.99.0 (2025-12-01)

* Submitted to Bioconductor
