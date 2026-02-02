
# Changes

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