# TIEGCM v3.0

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.20076374-blue)](https://doi.org/10.5281/zenodo.20076374)
[![Docs](https://img.shields.io/badge/Docs-readthedocs-blue)](https://tiegcm-docs.readthedocs.io/en/latest/)
[![Paper](https://img.shields.io/badge/JGR%3A%20Space%20Physics-10.1029%2F2025JA034219-blue)](https://doi.org/10.1029/2025JA034219)

The **NCAR Thermosphere–Ionosphere–Electrodynamics General Circulation Model (TIEGCM)** is a comprehensive, first-principles, three-dimensional, non-linear representation of the coupled thermosphere and ionosphere system that includes a self-consistent solution of the middle and low-latitude dynamo field.

It is developed and maintained by the [High Altitude Observatory](https://www2.hao.ucar.edu/) at NCAR.
## Contents

- [Quick Start](#quick-start)
- [Repository Layout](#repository-layout)
- [What's New in v3.0](#whats-new-in-v30)
  - [New Features and Functional Changes](#new-features-and-functional-changes)
  - [Changes in Physics](#changes-in-physics)
- [Data Files](#data-files)
- [Utility Tools](#utility-tools)
- [Documentation and Support](#documentation-and-support)
- [Citation](#citation)
- [License](#license)

## Quick Start

For a brief set of instructions to build the model and make a short default run, see the [TIEGCM documentation](https://tiegcm-docs.readthedocs.io/en/latest/). The recommended path is via [TIEGCMrun](#tiegcmrun), a Python tool that automates compilation and execution.

> **Note:** Download the [input and supporting data files](#data-files) before building or running the model.

## Repository Layout

| Subdirectory | Description           | Summary of Contents                                                |
|--------------|-----------------------|--------------------------------------------------------------------|
| `scripts/`   | Support scripts       | Job scripts, Make files, utilities                                 |
| `src/`       | Source code           | Source files (`*.F`, `*.F90`, `*.h`)                               |
| `tiegcmrun/` | TIEGCM user interface | Python tool for pre-processing, build, and job execution           |

## What's New in v3.0

This is a summary of modifications made to TIEGCM since the release of v2.0.

### New Features and Functional Changes

- **Flexible resolutions** — the job script supports arbitrary combinations of horizontal and vertical resolutions; magnetic-grid changes are also supported.
- **Extended upper boundary** — the job script, `defs.h`, and several altitude-dependent variables (`xfac`, `ar_glbm`, `aureff`, `bdriz`) have been rewritten to support raising the upper boundary.
- **High-cadence model time** — input/output timestamps are now 4 digits (day/hour/minute/second) instead of the old 3-digit format (day/hour/minute).
- **Unified N2 / MBAR / SCHT calculations** — the N2 mixing ratio, mean molecular mass, and scale height are computed once and shared. Artificial caps on the N2 mixing ratio have been removed.
- **Rewritten Helium module** — included at all resolutions (default on); Helium effects on heating rates and elsewhere are now accounted for consistently throughout the code.
- **Ring filter** — replaces the old Fourier filter.
- **O+ sub-cycling** — controlled by a new input parameter `NSTEP_SUB`.
- **NetCDF4 parallel I/O** — reduces memory usage on the root task.
- **Bit-for-bit reproducibility** — ESMF calls have been modified to ensure reproducibility.
- **Updated IGRF** — the geomagnetic field is updated to the latest International Geomagnetic Reference Field.
- **Rewritten magnetospheric coupling module** — supports in-memory MPI data transfer.
- **MPI subroutine optimizations** — several MPI subroutines have been rewritten for a speed boost.
- **Consistent `dipmin` calculation** — set to `sin(dlat*2*dtr)` instead of being manually set per resolution.
- **Code simplifications** — unused parameters, arguments, and variables removed from some functions.
- **Miscellaneous bug fixes**.

### Changes in Physics

- **Solar heating coefficients** — modified (*Astrid Maute*).
- **Height variation of equatorward electric field** — a scaling factor is added to account for the height variation of `elam` (*Astrid Maute*).
- **Field-aligned ion drag** — now included in the momentum equation (*Jiuhou Lei*).
- **Collision frequency** — `lamdas.F` now includes all ion species (O+, O2+, N+, N2+, NO+) instead of just O+, O2+, and NO+ (*Haonan Wu*).
- **N(2D) transport** — the minor-species solver now includes N(2D), previously assumed in (photo)chemical equilibrium; affects N chemistry at very high altitudes (z > 7) (*Haonan Wu*).
- **Electron heat flux parameterization** — the parameterization of `fed` near the equator is changed in `settei` (*Tong Dang, Wenbin Wang, Kevin Pham*).
- **O+ number flux parameterization** — the parameterization of `opflux` near the equator is changed in `oplus` (*Haonan Wu, Wenbin Wang*).
- **Thermal electron heating efficiency** — a sixth-order polynomial is now used (*Yihui Cai*).
- **Electrojet turbulent heating** — included, default off (*Jing Liu*).
- **Empirical SAPS** — included, default off (*Wenbin Wang*).
- **Eclipse solar EUV masking** — support added (*Tong Dang, Jiuhou Lei*).
- **Lower-boundary forcing by external data** — SD nudging support added (*Haonan Wu, Xian Lu*).

## Data Files

The input and supporting data files required for TIEGCM v3.0 are available here:

- [TIEGCM Data Files](https://bit.ly/4lddXZC)

These include the datasets needed to run the model and support various configuration options.

> Download and place these files in the appropriate location before building or running the model.

Additional data may be available on the [HAO public FTP site](http://download.hao.ucar.edu/pub/tgcm).

## Utility Tools

### TIEGCMrun

[TIEGCMrun](./tiegcmrun) is a Python tool used to compile and execute TIEGCM in an automated fashion. It can be run interactively on the command line. See [QuickStart](https://tiegcm-docs.readthedocs.io/en/latest/tiegcm/quickstart.html) for example usage.

The questions can start from an earlier run's answers. If the current directory holds a `tiegcmrun_parameters.json` (or a `tiegcmrun_input.json`), TIEGCMrun asks `Found <file>: start the wizard from its values? [Y/n]` and, on Enter, offers each saved answer as the default. `--prefill PATH` names the file (or a directory holding one) and does not ask; `--no-prefill` never looks. The question is asked only at a terminal; with piped input the file is not used (a NOTE says so) unless `--prefill` is given. Only answers are taken: values derived from them (vertres, mres, STEP, SOURCE, the history cadences and file names, the run directories, the model, data and make paths, the PE layout) are recomputed, and the run directory is not taken from the file. The saved executables (`modelexe`, `coupled_modelexe`) are offered again unless they are the saved run's own link or build; a job block kept per machine (`job.<machine>`, `job._common`) is read as with `-o`; a GPI or IMF file that was generated into the saved run is offered as `gen` again (a file for this run's dates), and a saved dated GPI/IMF file that does not cover the new dates is not offered; a saved walltime longer than the chosen queue allows gives way to the queue's default (`config/machines.yaml`). A bare data file name in a `tiegcmrun_input.json` is looked up under `$TIEGCMDATA`, as with `-o`. A prefilled KP, POWER, CTPOTEN, F107 or F107A that the chosen forcing replaces (a GPI file, for example) is not used; a NOTE names it. Under `--engage` the current directory is not searched: `engage.py` passes `--prefill` or `--no-prefill`. Unknown keys are skipped with a NOTE. A value that fails its check is replaced by the default with a NOTE at a terminal and refused with one ERROR when stdin is not a terminal. Under `--engage` a refused executable is asked again at most three times (as the forcing prompts); then, or at end of input, the run stops with the refusal. When the input ends at any prompt, tiegcmrun stops with one ERROR naming the option. `--prefill` cannot be combined with `-o`, which replays a file as it is.

`-o` and engage `-to` replays are checked like the prompts: a dated GPI/IMF file must cover the run (set it to 'gen' or a covering file; a file generated by this run only gets the WARNING when its data stop short of the run), KP with both POWER and CTPOTEN, an empty SECFLDS list and a key given twice are refused, as are a queue of another machine with `-o` (under engage the queue is engage's) and a null inp/job section under engage, a bare data-file name is looked up under `$TIEGCMDATA`, and with `--check` a pending 'gen' is refused before anything is written when gcmprocpy cannot be imported. Job fields of another machine (for example `job_priority` on aitken) are dropped with a NOTE, with `-o` and under engage. When simulation.hpc_system is another cluster than the host, the data files that were looked up under this host's $TIEGCMDATA (and those not found under that machine's tgcmdata) are named in a WARNING; set them to that machine's paths or generate the run there.

Under `--engage`, when the coupled grid differs from the spin-up grid, the last spin-up job runs the hand-off script in the Python environment named by engage's `coupling.conda_env`. A path to a venv (a directory holding `bin/activate`) is activated with the machine's `venv_activate` lines in `config/machines.yaml` (`source <path>/bin/activate.csh` in the tcsh jobs of derecho and aitken, after setting the `prompt` variable that a batch tcsh lacks, `source <path>/bin/activate` in bash jobs); any other value is a conda environment, activated with the machine's `conda_activate` lines. A machine without the lines needed stops the run with one ERROR naming the missing key. `machine_env.py conda <path>` prints the bash form of the same venv lines (`venv_activate_sh`, else `venv_activate`).

### GCMProcpy

[GCMProcpy](https://github.com/NCAR/gcmprocpy) is a Python tool for post-processing and visualization of TIEGCM output. It can be used interactively on the command line or as an API in a Python script. See the [GCMProcpy docs](https://gcmprocpy.readthedocs.io/en/latest/) for examples.

## Documentation and Support

- **User's Guide, Model Description, and Release notes:** [TIEGCM ReadtheDocs](https://tiegcm-docs.readthedocs.io/en/latest/)
- **TGCM website:** http://www.hao.ucar.edu/modeling/tgcm
- **Discussion / questions:** tgcmgroup@ucar.edu

## Citation

If you use TIEGCM in your work, please cite **both** the software and the accompanying paper. Citation metadata is in [`CITATION.cff`](./CITATION.cff).

- **Software (this release):** [10.5281/zenodo.20076374](https://doi.org/10.5281/zenodo.20076374)
- **Paper:** Wu, H., Wang, W., Pham, K. H., et al. (2025). *The NCAR-TIEGCM Version 3.0.* Journal of Geophysical Research: Space Physics. [10.1029/2025JA034219](https://doi.org/10.1029/2025JA034219)

## License

TIEGCM is distributed under the **NCAR TIE-GCM Open Source Academic Research License Agreement**. The full text is in [`LICENSE`](./LICENSE) and applies to v2, v3, and subsequent versions.
  
