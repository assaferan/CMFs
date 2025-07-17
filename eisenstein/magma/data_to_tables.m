import "../../magma/mfpgdata.m" : newspaces_columns, gamma1_columns, hecke_traces_columns, 
                                  newforms_columns, hecke_nf_columns, hecke_lpolys_columns, 
                                  twists_nf_columns, twists_cc_columns, hecke_cc_columns, char_dir_orbits_columns;

import "../../magma/mfpgdata.m" : FormatNewspaceData, CreateDimensionTableTwo, 
                                  FormatNewformData, FormatHeckeCCData;

import "../../magma/mf.m" : DecomposeSpace;

input_columns := ["level", "weight", "char_orbit", "time", "dims", "traces", "AL_signs", "hecke_fields",
                  "hecke_cutters", "hecke_ev_data", "cm", "twists", "is_polredabs", "zero_rate", "moments", 
                  "trace_hashes", "char_gens_values", "self_dual", "embedded_ap"];

space_input_columns := ["level", "weight", "char_orbit", "time", "dims", "traces", "hecke_fields", "AL_signs", "hecke_cutters"];

dim_columns := ["level", "weight", "char_orbit", "dim_cusp", "dim_eis", "dim_cusp_new", "dim_eis_new"];

/*
* 
* Need to generate file with conrey_labels !?
* 
* Writing for myself the flow of the Magma code (hypothesized):
* 0) Set names for infile, infile_dim, dimfile, newspace_outfile, gamma1_outfile, trace_outfile, 
      outfile_prefix, outfile_suffix, outfile_hecke_cc, 
* 1) Run DecomposeSpace(infile_dim, B : DimensionsOnly) where B is the bound on Nk^2.
* 2) Run DecomposeSpace(infile, B)
* 3) Run CreateDimensionTable(infile_dim, dimfile) // Also has CreateDimensionTableTwo which gets N:k:i as input?
* 4) Run FormatNewspaceData(infile, newspace_outfile, gamma1_outfile, trace_outfile, dimfile)
* 5) Run FormatNewformData(infile, outfile_prefix, outfile_suffix)
* 6) Run FormatHeckeCCData(infile, outfile_hecke_cc) ?
* 
* Missing pieces - 
* 1) input data for FormatTwistDataNF (character orbit labels and multiplicities) ?
* 2) input data for FormatTwistDataCC (embedded newform labels and conrey labels) ?
* 3) conrey_file (conrey_XXX.txt) ? 
* 4) input for FormatHeckeCCData ?
* 5) "lmfdb_nf_labels.txt" - coeffs:label for number fields from the LMFDB ?
* 6) "mf_ranks.txt" - label:rank:proved ?
* 7) "mf_twists_inner.txt" - source:target:chars:mults ?
* 8) "mf_twists_minimal.txt" - source:target ?
*    "mf_twists_cc_minimal.txt" - source:target ?
* 9) "mf_related_objects.txt" - label:URL ?
* 10) "mfsplit.txt" ? 
*/

// !! TODO - some functions need more data than that
// format is N:o:[n1,n2,...]:M:pmo:order:deg:parity:is_real (list of conrey chars chi_N(n,*) in orbit o, M=cond, pmo=primitive_orbit_index
procedure WriteConreyLabelsFile(outfile, maxN : minN := 1, Quiet := false)
    fp := Open(outfile, "w");
    A := AssociativeArray();
    for N in [minN..maxN] do
        if not Quiet then printf "Constructing character orbit table for modulus %o...", N; t:=Cputime(); end if;
        G := FullDirichletGroup(N);
        orbit_table := AssociativeArray();
        conductors := AssociativeArray();
        pmo := AssociativeArray();
        T := AssociativeArray();
        orders := AssociativeArray();
        degs := AssociativeArray();
        pars := AssociativeArray();
        is_real := AssociativeArray();
        for chi in Elements(G) do
            o := CharacterOrbit(chi);
            if not IsDefined(orbit_table, o) then orbit_table[o] := []; end if;
            if not IsDefined(conductors, o) then conductors[o] := Conductor(chi); end if;
            orbit_table[o] := Append(orbit_table[o], ConreyIndex(chi));
            if not IsDefined(pmo, o) then 
                T[chi] := o;
                M := Conductor(chi);
                pmo[o] := (M eq N) select o else A[M][AssociatedPrimitiveCharacter(chi)];
            end if;
            if not IsDefined(orders, o) then orders[o] := Order(chi); end if;
            if not IsDefined(degs, o) then degs[o] := Degree(chi); end if;
            if not IsDefined(pars, o) then pars[o] := Parity(chi); end if;
            if not IsDefined(is_real, o) then is_real[o] := (IsReal(chi) select 1 else 0); end if;
        end for;
        A[N] := T;
        for o in Sort([k : k in Keys(orbit_table)]) do
            str := Sprintf("%o:%o:%o:%o:%o:%o:%o:%o:%o", N, o, orbit_table[o], conductors[o], pmo[o],orders[o],degs[o],pars[o],is_real[o]);
            Puts(fp, str);
            Flush(fp);
        end for;
        if not Quiet then printf "took %o secs\n",Cputime()-t; end if;
    end for;
    delete fp;
    return;
end procedure;

procedure DoEverythingNk2UpTo(B : folder := "../data/")
    infile := folder cat "Nk2_";
    infile_dim := infile cat "dim_only";
    dimfile := infile cat "dim";
    suffix := Sprintf("%o.m.txt", B);
    newspace_outfile := infile cat "mf_newspaces_" cat suffix;
    gamma1_outfile := infile cat "mf_gamma1_" cat suffix;
    trace_outfile := infile cat "mf_hecke_traces_" cat suffix;
    outfile_prefix := folder cat "Nk2_";
    outfile_hecke_cc := outfile_prefix cat "mf_hecke_cc_" cat suffix;
    conrey_labels := infile cat "conrey_" cat suffix;
    WriteConreyLabelsFile(conrey_labels, B);
    // DecomposeSpace(infile_dim, B : DimensionsOnly, Timings := false);
    DecomposeSpace(infile, B : Eisenstein);
    // CreateDimensionTable(infile_dim, dimfile);
    CreateDimensionTableTwo(dimfile, B);
    FormatNewspaceData(infile, newspace_outfile, gamma1_outfile, trace_outfile, dimfile : conrey_labels := conrey_labels, Eisenstein);
    FormatNewformData(infile, outfile_prefix, suffix : conrey_labels := conrey_labels, Eisenstein);
    FormatHeckeCCData(infile, outfile_hecke_cc : Eisenstein, conrey_labels := conrey_labels);
    return;
end procedure;


// This procedure takes as input the output of DecomposeSpace from mf_eis.m
// and generates a file whose lines are formatted as N:k:o:time:dims:traces:hecke_fields:AL_signs:cutter_data 
// for input to FormatNewSpaceData
procedure MakeSpaceInputFile(infile, space_infile)
    f := Open(infile);
    lines := Split(Read(f),"\n");
    delete f;
    fp := Open(space_infile, "w");
    for line in lines do
        input_values := Split(line, ":");
        space_dict := AssociativeArray();
        for j->fld in input_columns do
            if fld notin space_input_columns then continue; end if;
            space_dict[fld] := input_values[j];
        end for;
        outline := Join([space_dict[fld] : fld in space_input_columns], ":");
        Puts(fp,outline);
        Flush(fp);
    end for;
    delete fp;
end procedure;

// Now producing the dimension file together with the input file
/*
procedure MakeDimFile(infile, dim_file)
    f := Open(infile);
    lines := Readlines(f);
    delete f;
    fp := Open(dim_file, "w");
    for line in lines do
        input_values := Split(line, ":");
        dim_dict := AssociativeArray();
        for j->fld in input_columns do
            if fld not in dim_columns then continue; end if;
            dim_dict[fld] := input_values[j];
        end for;
        dim_dict["dim_cusp"] := QDimensionCuspForms();
        outline := Join([dim_dict[fld] : fld in dim_input_columns], ":");
        Puts(fp,outline);
        Flush(fp);
    end for;
    delete fp;
end procedure;
*/

procedure DataToTables(infile, dimfile, outfile_prefix)
    suffixes := ["newspaces", "gamma1", "hecke_traces", "newforms", "hecke_nf", "hecke_lpolys", 
                 "twists_nf", "twists_cc", "hecke_cc", "char_dir_orbits"];
    outfiles := AssociativeArray();
    columns := AssociativeArray();
    for suffix in suffixes do
        outfiles[suffix] := outfile_prefix cat suffix;
        columns[suffix] := eval(suffix cat "_columns");
    end for;
   
    // !! TODO - need to create a dimension file with lines N:k:o:dS:dE:dNS:dNE
    space_infile := infile cat "_space";
    MakeSpaceInputFile(infile, space_infile);
    FormatNewspaceData(space_infile, outfiles["newspaces"], outfiles["gamma1"], outfiles["hecke_traces"], dimfile);
    FormatNewformData (infile, infile, "m.txt");
end procedure;