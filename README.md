# Atoll seabird analysis

Text about main findings 

(Tried adding the summary figure, but this currently contains a lot of whitespace at the bottom. Can add it back later after removing whitespace)

## Workflow

This would be nice, but we did not go with a figure for this.
Drop or describe in text format?

## Installation

### System-independent prerequisites

Install the R programming language ([Windows](https://cran.r-project.org/bin/windows/), [MacOS](https://cran.r-project.org/bin/macosx/), [Linux](https://cran.r-project.org/bin/linux/)), and the [Julia programming language](https://julialang.org/downloads/).

### Windows-specific prerequisites

Install an environment that allows you to run bash scripts, like [Git Bash](https://gitforwindows.org/)

### Unix-based systems

Unix-based systems support executing bash scripts natively.
However, since the default shell used by MacOS is now the Z shell, we recommend running the commands under [Reproducing results](#reproducing-results) in the bash shell.

## Reproducing results

You can reproduce the analysis and corresponding visualisations by executing the `reproduce.sh` script.
To do so, start by making this file executable on your machine:

```console
$ chmod +x /path/to/repository/reproduce.sh
```

Replace the segment `/path/to/repository` below with the path to the project folder on your machine.

Next, run all scripts by executing `reproduce.sh`:

```console
$ ./reproduce.sh /path/to/repository true true
```

Replace the segment `path/to/repository` with the path to the project folder on your machine.

**Note on Z shell:**
To run the above command in Z shell, prefix with `bash`, i.e., `% bash ./reproduce.sh ...`.

**Additional arguments:**
The extra arguments to `reproduce.sh` (here, `true true`) are forwarded to the Julia scripts.
The first argument determines if the analysis scripts sample from the posterior (`true`) or attempt to load previously saved chains (`false`). Loading saved chains requires having sampled from the posterior at least once on your machine.
The second argument determines if [cross validation](https://mc-stan.org/loo/articles/online-only/faq.html) is performed (`true`) or skipped (`false`).
Both options are intended to allow the user to quickly rerun the analysis after having sampled from the posterior once.

**Run-time:**
The analyses were performed on a MacBook Pro (M1) and took *TODO measure time!*.

(If we / reviewers think that this runtime is too long, we can make running sensitivity analyses optional with another argument to `reproduce.sh`)

## Project structure

A sentence or two about the project structure.

```
├── R                            # R scripts
│   ├── wrangle                      # R scripts used for data wrangling
│   └── create_figures.R             # Creates raw figures for article
├── data                         # CSV files which are *inputs to* the model
├── figures                      # Final figures used in the article
├── julia                        # Julia scripts
│   ├── scripts                      # Analysis scripts
│   ├── src                          # Modules defining functions and variables
│   └── reproduce.jl                 # Runs all julia scripts
├── manuscript                   # Directory with manuscript
├── renv                         # Renv for storing R package versions
├── results                      # Outputs of modeling
│   ├── chains                       # Chains will be saved here
│   ├── data                         # CSV files which are *outputs of* the model
│   ├── png                          # PNG files
│   └── svg                          # SVG files
.
.
.
└── reproduce.sh                 # Execute entire model pipeline, see "Installation" above for instructions
```
