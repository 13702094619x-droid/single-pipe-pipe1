# Single-Pipe Pipe1 Case

This repository contains the single-pipe heating-network case used for the
SAS/FDM comparison in the local `SAS_Pipe1.m` workflow.

## What Is Included

- `single_pipe_case/`
  - MATLAB scripts for the single-pipe SAS and FDM calculations.
  - `single pipe.xlsx`, the local boundary input table.
  - `P_out_*.mat` and `T_out_300.mat`, pressure/temperature boundary inputs used by the scripts.
- `supplementary_data/`
  - `BC_case_DHN-1.xlsx`, copied from the paper supplementary data.
  - `Topological data_DHN-1.xlsx`, copied from the paper supplementary data and modified so the heat coefficient matches `SAS_Pipe1.m`.

## Heat-Loss Setting

The paper supplementary data gives the DHN-1 heat coefficient as:

```text
0.67 W/(m*K)
```

The local `SAS_Pipe1.m` case uses:

```matlab
data_lamda_w = 0.67*5;
```

Therefore this repository uses:

```text
3.35 W/(m*K)
```

in `supplementary_data/Topological data_DHN-1.xlsx`, and the FDM script copies in
`single_pipe_case/` have been aligned to the same setting.

## Suggested Run Order

Run the MATLAB scripts from inside `single_pipe_case/`:

```matlab
FDM_pipe1_liangtiao
FDM_pipe1_AE_liangtiao
SAS_Pipe1
```

Generated result files such as `SAS_Pipe1_NCI_results.mat` are ignored by Git.

## Source Note

The DHN-1 supplementary data was originally published at:

```text
https://github.com/lyq-0106/Case-data
```

This repository keeps a modified local case version for the SAS Pipe1 heat-loss
setting.
