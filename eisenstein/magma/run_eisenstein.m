/*
  Parameterized entry point for regenerating the Eisenstein CMF data tables.
  Wraps DoEverythingNk2UpTo in data_to_tables.m so it can be driven from the
  command line by run_parallel.sh (job-split across many cores).

  Must be run from the eisenstein/magma/ directory (the imports in
  data_to_tables.m use paths relative to it, and folder defaults to ../data/).

  Command-line variables (Magma `name:=value` syntax), all optional except B:
      B             Nk^2 bound (required),                    e.g. B:=4000
      Jobs          number of parallel workers,               default 1
      JobId         this worker's id in [0..Jobs-1],          default 0
      DecomposeOnly true  -> run only DecomposeSpace (writes  <infile>_<JobId>)
      FormatOnly    true  -> run only dim table + Format*      (reads assembled <infile>)
      folder        output folder (with trailing slash),      default "../data/"
      do_conrey     true  -> (re)generate the conrey labels file first, default false

  Examples:
      magma B:=4000 Jobs:=64 JobId:=7 DecomposeOnly:=true folder:="../data/" run_eisenstein.m
      magma B:=4000 FormatOnly:=true folder:="../data/" run_eisenstein.m
      magma B:=100  run_eisenstein.m                 // plain single-process full run
*/

// chars.m (via the spec) supplies ConreyCharacterOrbitReps etc. used by
// data_to_tables.m and its imported mf.m / mfpgdata.m. Paths are relative to
// eisenstein/magma/, which must be the working directory.
AttachSpec("../../magma/mf.spec");
load "data_to_tables.m";

error if not assigned B, "Missing required command-line variable B (Nk^2 bound), e.g. B:=4000.";

// Magma passes command-line `name:=value` as strings; coerce to the intended types.
// (Type checks keep this correct even if a value is supplied as a literal integer/bool.)
if Type(B) eq MonStgElt then B := StringToInteger(B); end if;
if not assigned folder then folder := "../data/"; end if;
if not assigned Jobs then Jobs := 1; elif Type(Jobs) eq MonStgElt then Jobs := StringToInteger(Jobs); end if;
if not assigned JobId then JobId := 0; elif Type(JobId) eq MonStgElt then JobId := StringToInteger(JobId); end if;
if not assigned do_conrey then do_conrey := false; elif Type(do_conrey) eq MonStgElt then do_conrey := (do_conrey eq "true"); end if;
if not assigned DecomposeOnly then DecomposeOnly := false; elif Type(DecomposeOnly) eq MonStgElt then DecomposeOnly := (DecomposeOnly eq "true"); end if;
if not assigned FormatOnly then FormatOnly := false; elif Type(FormatOnly) eq MonStgElt then FormatOnly := (FormatOnly eq "true"); end if;

printf "run_eisenstein: B=%o folder=%o Jobs=%o JobId=%o DecomposeOnly=%o FormatOnly=%o do_conrey=%o\n",
       B, folder, Jobs, JobId, DecomposeOnly, FormatOnly, do_conrey;

DoEverythingNk2UpTo(B : folder := folder, do_conrey := do_conrey,
                    Jobs := Jobs, JobId := JobId,
                    DecomposeOnly := DecomposeOnly, FormatOnly := FormatOnly);

exit;
