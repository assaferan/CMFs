# Regenerating the Eisenstein CMF data tables

This regenerates all `mf_*_eis` data tables for `Nk^2 <= B` with the corrected
Eisenstein constant term (`a0_num`), fixing [LMFDB/lmfdb#7034](https://github.com/LMFDB/lmfdb/issues/7034):
Eisenstein forms with a non-rational Hecke field previously stored `a0_num = 0`
(so the newform page showed no constant term). The generation fix is CMFs
commit `a4b8403`; this run reproduces the data with it.

## Files

- `data_to_tables.m` — `DoEverythingNk2UpTo` now takes `Jobs`/`JobId`/`DecomposeOnly`/`FormatOnly`
  so the expensive `DecomposeSpace` step can be split across cores. `Jobs=1` (the default) is
  the original single-process behavior, unchanged.
- `run_eisenstein.m` — command-line entry point (attaches `../../magma/mf.spec`, then runs
  `DoEverythingNk2UpTo`). Run it **from `eisenstein/magma/`**.
- `run_parallel.sh` — orchestrator: N parallel `DecomposeSpace` workers → concat+sort → format → manifest.

## Running on a many-core server

Works on any server with Magma and enough cores/disk (pick `NJOBS` to match the
available cores). Run it under `tmux`/`screen` so it survives disconnects.

```sh
cd CMFs/eisenstein/magma
tmux new -s eisregen                     # so it survives disconnects
./run_parallel.sh 4000 64 ../data        # B=4000, 64 workers, output to ../data
```

Phases (all logged under `eisenstein/magma/logs/`):
1. 64 Magma workers, each `DecomposeSpace` on its slice → `../data/Nk2_4000.m.txt_<i>`
2. concatenate + sort by `(N,k,char_orbit)` → `../data/Nk2_4000.m.txt` (matches a serial run)
3. one Magma process: dimension table + `Format*` → the `Nk2_mf_*_4000.m.txt` tables
4. writes `../data/MANIFEST_Nk2_4000.txt` (file list with sizes + line counts)

Notes:
- The conrey-labels file `../data/Nk2_conrey_4000.m.txt` already exists and is reused
  (`DO_CONREY` auto-detects; set `DO_CONREY=true ./run_parallel.sh ...` to force regeneration).
- `mf_hecke_cc_4000.m.txt` is large (~17 GB); ensure enough disk in `../data`.
- Env overrides: `MAGMA=/path/to/magma`, `KEEP_SLICES=true`, `SORT_TMP=/scratch`.

## Validation already done (B=100, local)

- Regenerated `a0_num` is now nonzero **and consistent** with `trace_a0_num`
  (e.g. `7.3.E.d.a`: `a0_num={-1,1}`, `Tr = -1 = trace_a0_num`); `qexp_display` shows the constant term.
- Parallel (`NJOBS=4`) vs serial output: **all tables byte-identical** (infile identical modulo the per-space timing field).

## Loading the results

The output `Nk2_mf_*_<B>.m.txt` files are colon-separated with a 2-line header
(column names, then Postgres types), matching the existing `mf_*_eis` tables:

| file | table |
|------|-------|
| `Nk2_mf_newspaces_<B>.m.txt`      | `mf_newspaces_eis` |
| `Nk2_mf_gamma1_<B>.m.txt`         | `mf_gamma1_eis` |
| `Nk2_mf_newforms_<B>.m.txt`       | `mf_newforms_eis` |
| `Nk2_mf_hecke_nf_<B>.m.txt`       | `mf_hecke_nf_eis`  ← carries the fixed `a0_num` |
| `Nk2_mf_hecke_traces_<B>.m.txt`   | `mf_hecke_traces_eis` |
| `Nk2_mf_hecke_lpolys_<B>.m.txt`   | `mf_hecke_lpolys_eis` |
| `Nk2_mf_hecke_charpolys_<B>.m.txt`| `mf_hecke_charpolys` |
| `Nk2_mf_hecke_cc_<B>.m.txt`       | `mf_hecke_cc_eis` |

After loading, spot-check on a previously-broken form:
`.../ModularForm/GL2/Q/holomorphic/7/3/E/d/a/` should show a nonzero constant term,
and `mf_hecke_nf_eis.a0_num` for `7.3.E.d.a` should be `[-1,1]` (not `[0,0]`).
