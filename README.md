# Single-Pipe Pipe1 Case

This repository contains the single-pipe heating-network case used for local
dynamic calculation comparisons.

## What Is Included

- `single_pipe_case/`
  - MATLAB scripts for the single-pipe calculation workflows.
  - `single pipe.xlsx`, the local boundary input table, including inlet flow,
    inlet temperature, outlet/initial temperature, and the inlet `NCI_source`
    boundary condition.
  - `P_out_*.mat` and `T_out_300.mat`, pressure/temperature boundary inputs used by the scripts.
- `supplementary_data/`
  - `BC_case_DHN-1.xlsx`, copied from the paper supplementary data.
  - `Topological data_DHN-1.xlsx`, copied from the paper supplementary data and modified so the heat coefficient matches the local Pipe1 setting.

## Heat-Loss Setting

The paper supplementary data gives the DHN-1 heat coefficient as:

```text
0.67 W/(m*K)
```

The local Pipe1 case uses:

```matlab
data_lamda_w = 0.67*5;
```

Therefore this repository uses:

```text
3.35 W/(m*K)
```

in `supplementary_data/Topological data_DHN-1.xlsx`, and the script copies in
`single_pipe_case/` have been aligned to the same setting.

## Suggested Run Order

Run the MATLAB scripts from inside `single_pipe_case/`. Generated result files
are ignored by Git.

## NCI Boundary

`single_pipe_case/single pipe.xlsx` includes an explicit inlet NCI boundary
column:

```text
NCI_source
```

The value series follows the original waveform used inside the scripts, scaled
to the range `0.2` to `0.3`. The MATLAB helper
`single_pipe_case/load_nci_source_boundary.m` reads this column first. If the
column is missing, it falls back to the historical default waveform.

## Source Note

The DHN-1 supplementary data was originally published at:

```text
https://github.com/lyq-0106/Case-data
```

This repository keeps a modified local case version for the Pipe1 heat-loss
setting.
