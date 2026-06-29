
# Changes

## [1.2.0] - 2026-06-29

### Features

- New option 'rename' to relabel variables for display ('orig -> Display') and/or rename factor levels in the data ('Var: old -> new, ...') in a single string.
- Replaced the chocolate/reaction-time demos with a systematic demo set (demo_01 ... demo_12) ported from kbstatpy, each on a classic R dataset and illustrating one concept (t-tests, LMM, random slopes, transforms, gamma/binomial GLMM, partial interaction, outliers, multiple dependent variables). Added run_demos.m to run them all.

### Changes

- Pre- and post-fit outliers are now set to NaN in place instead of removing their rows. This keeps each data point's position in the data array fixed (which encodes its identity), and is statistically equivalent because the model fit ignores NaN responses.
- Demo output is no longer tracked in the repository (gitignored under demo/results/).

### Bugs

- connectPoints in the violin plots now works with unbalanced data (unequal group sizes, e.g. after outlier removal), pairing points by their data slot instead of over-indexing.

## [1.1.6] - 2026-02-02

### Bugs

- When options.xOrder or options.xOrder<n> were provided, kbstat would exit with an error.

## [1.1.5] - 2025-12-16
### Features

- New option 'connectPoints' to draw a line between associated data points in the violin plots. This only works when the group sizes are equal, and it should only be applied when the data points in the group actually correspond, i.e. for paired data.
- Added possibility to not plot the data points by setting 'markerSize' to zero.

### Bugs

- When posthocLevel > 1, data groups were wrongly labeled, and an error would be thrown when the number of levels would not be equal.

## [1.1.4] - 2025-06-12

### Features

- Enhanced options.xOrder to also use group names rather than numerical
- Added options.plotGroupSize to plot group size (N=...) in DataPlots
- Added options.closeFigures to flag if the figures created during the analysis should be closed or left open
- Added this CHANGELOG file to the repository

### Bugs

- Fixed postHocTable: same row/col/group was printed for all rows/cols/groups

### Other

- Allow options.levelOrder to contain uppercase letters

## [1.1.3] - 2025-06-06

### Features

- Using 'transform' now induces back-transformation of data in plots and tables
- Added boxcharts as plot option

### Bugs

- Fixed bug when levels are numerical, and include trial random variable in renaming scheme
- Fixed a weird bug that emmeans calculates wrong values if some levels are substrings of other levels
- Fix flawed results when data table contains numeric data for fixed effect levels
- Fixed wrong plot sorting

### Other

- Simplify options
- Updated demo files

## [1.1.2] - 2024-12-13

-  kbstat_demo_live.mlx: Set new option 'closeFigures' to false so that figures show in live script

## [1.1.1] - 2024-12-13

- kbstat: FIX covariate should not be in factor list
- excel2csv: Change delimiter of CSV to semicolon to support Excel opening the file directly
- tableLong2Wide, tableWide2Long: Bug fixes and improvements
- Updated README and demo output files

## [1.1.0] - 2024-11-12

- Improve option showVarNames
- Improve demo: level display
- Output argument is now a structure with the relevant items
- If no outDir is specified or empty, do not save anything
- New option 'closeFigures' to close figures after creating them

## [1.0.0] - 2024-11-06

- Project published