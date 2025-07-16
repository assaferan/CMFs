import "../../magma/mf.m" : ComputeHeckeCutters, InnerTwistData, SortEmbeddings, IsSelfDual;


/*
Format of data is N:k:i:t:D:T:A:F:C:E:cm:tw:pra:zr:mm:h:X:sd:eap
The data depends on a degree bound (determines forms with exact eigenvalue data, a coefficient count (number of a_n to compute), and a complex precision (for forms without exact eigevnalue data)

N = level, a positive integer
k = weight, a positive integer (for .m.txt files, k > 1)
i = character orbit of chi (Galois orbits of Dirichlet characters chi of modulus N are lex-sorted by order and then by trace vectors [tr(chi(n)) for n in [1..N]], taking traces from Q(chi) to Q; the first orbit index is 1, corresponding to the trivial character, the second orbit will correspond to a quadratic character).
t = time in secs to compute this line of data
D = sorted list of dimensions [d1,d2,...] of Galois stable subspaces of S_k^{new}(N,chi), ordered by dimension
T = lex-sorted list of trace vectors [[tr(a_1),...tr(a_n)],...] for Galois conjugacy classes of eigenforms f corresponding to the subspaces listed in D, traces are from the coefficient field of the form down to Q (note that lex-sorting trace vectors sorts subspaces by dimension because tr(a_1)=tr(1) is the degree of the coefficient field)
A = Atkin-Lehner signs (empty list if chi is not the trivial character (i.e. i=1)) [[<p,sign> for p in Divisors(N)],...], one list of Atkin-Lehner signs for each subspace listed in D.
F = Hecke field polys [[f0,f1,...,1],...] list of coeffs (constant coeff first), one list for each subspace listed in D of dimension up to the degree bound (currently 20); note that F[n] corresponds to the space D[n] but F may be shorter than D
C = Hecke cutters [[<p,[g0,g1,...,1]>,...],...] list of minimal lists of coefficients of charpolys g(x) of T_p sufficient to distinguish all the subspaces listed in D up to the degree bound.
E = Hecke Eigenvalue data [<g,b,n,m,e>,...] list of tuples <g,b,n,m,e> of Hecke eigenvalue data for each subspace listed in D of dimension greater than 1 up to the degree bound where:
    g is a polredbestified field poly for the coefficient field (should be the same as the corresponding poly in F),
    b is a basis for the Hecke ring R:=Z[a_n] in terms of the power basis of K:=Q[x]/(f(x)) (a list of lists of rationals),
    n is an integer that divides the index [O_K:R] of the Hecke ring R in the ring of integers O_K
    d is a pair that is either <0,[]> or <D,[<p,e>]> giving the discriminant of the Hecke ring and its prime factorization (if known)
    e is a list of eigenvalues specified in terms of the basis b (list of deg(f) integers for each a_n)
    x is a pair <u,v> where u is a list of integers generating Z/NZ* and v is a list of values of chi on u in written in the basis b
    m is the list integer such that teh first m eigenvalues generate the Hecke ring (as a ring)
cm = list of cm discriminants, one for each subspace listed in D up to the degree bound, 0 indicates non-CM forms (rigorous)
tw = list of lists of quadruples <b,n,m,i> identifying char orbits m.i of non-trivial inner twists with multiplicity n, b=0,1 indicates proved or not
pra = list of boolean values (0 or 1) such that pra[i] is 1 if F[i] is the polredabs polynomial for the Hecke field
zr = list of proportions of zero a_p over primes p < 2^13 (decimal number), one for each subspace
mm = list of list of moments of normalized a_p over primes p < 2^13 (decimal numbers), one for each subspace
h = list of trace hashes (linear combination of tr(a_p) over p in [2^12,2^13] mod 2^61-1), one for subspace
X = list of pairs <u,v> one for each entry in F where u is a list of integers r generating Z/NZ* and v is a list of values of chi on r in power basis defined by corresponding entry in F
sd = list of booleans, one for each entry in D, indicating whether newform is self dual or not (i.e. a_n are real)
eap = list of lists of lists of real or complex valued a_p's for p up to the coefficient bound for each embedding of each form where exact eigenvalues have not been computed
      if character is trivial embedded a_p's will always be real (this is actually the only case currently used)

This format is also documented in https://github.com/JohnCremona/CMFs/blob/master/README.md.
*/


function NewspaceData(chi, k, o: CharTable:=AssociativeArray(), TraceHint:=[], DimensionsOnly:=false, ComputeEigenvalues:=false, ComputeTwists:=false, ComputeTraceStats:=false,
                       NumberOfCoefficients:=0, DegreeBound:=0, EmbeddingPrecision:= 0, Detail:=0, ReturnDecomposition:=false, ComputeCutters:=false, 
                       ComputeCharacterValues:=true, Timings := true)
    start := Cputime();
    if o eq 0 then o := CharacterOrbit(chi); end if;
    N := Modulus(chi);
    dNS := QDimensionNewEisensteinForms(chi,k);
    if dNS eq 0 then
        if Detail gt 0 then printf "The space %o:%o:%o is empty\n",N,k,o; end if;
        s := Sprintf("%o:%o:%o", N, k, o);
        if Timings then s cat:= Sprintf(":%o", Cputime()-start); end if;
        s cat:= ":[]";
        if not DimensionsOnly then s cat:= ":[]:[]:[]:[]:[]:[]:[]:[]:[]:[]:[]:[]:[]:[]"; end if;
        return strip(s);
    end if;
    if Detail gt 0 then printf "Constructing space %o:%o:%o...", N,k,o; t:=Cputime(); end if;
    NS := NewSubspace(EisensteinSubspace(ModularSymbols(chi,k,0)));
    if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
    if NumberOfCoefficients eq 0 then
        if N le 1000 then NumberOfCoefficients := 1000; end if;
        if N gt 1000 and N le 4000 then NumberOfCoefficients := 2000; end if;
        if N gt 4000 and N le 10000 then NumberOfCoefficients := 3000; end if;
    end if;
    n := Max([SturmBound(N,k)+1,Floor(30*Sqrt(N)),NumberOfCoefficients]);
    assert QDimension(NS) eq dNS;
    if #TraceHint gt 0 then
        assert &and[t[1] ge 1:t in TraceHint] and &+[t[1]:t in TraceHint] eq dNS;
        assert #Set([#t:t in TraceHint]) eq 1;
        if #TraceHint eq 1 and DimensionsOnly then
            if Timings then 
                return strip(Sprintf("%o:%o:%o:%o:%o:", N, k, o, Cputime()-start, [TraceHint[1][1]]));
            else
                return strip(Sprintf("%o:%o:%o:%o:", N, k, o, [TraceHint[1][1]]));
            end if;
        end if;
        TraceHint := Sort(TraceHint);
        if #TraceHint[1] lt n then
            printf "*** Warning: ignoring TraceHint because it contains only %o < %o traces! ***", #TraceHint[1], n;
            TraceHint := [];
        end if;
        if #TraceHint eq 1 and DegreeBound gt 0 and dNS gt DegreeBound and EmbeddingPrecision eq 0 then
            if Detail gt 0 then printf "TraceHint implies that the space %o:%o:%o consists of a single orbit of dimension %o\n",N,k,o,dNS; end if;
            if Detail gt 0 and Order(chi) eq 1 then printf "Computing Atkin-Lehner signs for space %o:%o:%o...", N,k,o; t:=Cputime(); end if;
            AL := Order(chi) eq 1 select [[<p,ExactQuotient(Trace(AtkinLehnerOperator(NS,p)),dNS)>:p in PrimeDivisors(N)]] else [];
            if Detail gt 0 and Order(chi) eq 1 then printf "took %o secs.\n", Cputime()-t; printf "Atkin-Lehner signs %o\n", sprint(AL); end if;
            if Detail gt 0 then printf "Checking whether single form in space %o:%o:%o has CM...", N,k,o; t:=Cputime(); end if;
            cm := [<1,SelfTwists([],NS:TraceHint:=TraceHint[1],pBound:=SturmBound(N,k))>];
            if Detail gt 0 then printf "took %o secs.\n", Cputime()-t; printf "CM: %o\n", cm; end if;
            if Timings then
                s := Sprintf("%o:%o:%o:%o:%o:%o:%o:[]:[[]]:[]:%o:[]:[]", N, k, o, Cputime()-start, [dNS], TraceHint, AL, cm);
            else
                s := Sprintf("%o:%o:%o:%o:%o:%o:[]:[[]]:[]:%o:[]:[]", N, k, o, [dNS], TraceHint, AL, cm);
            end if;
            s cat:= Sprintf(":[]:[]:[]:[]:[%o]:[]",IsSelfDual(chi,dNS,TraceHint,[],NS) select 1 else 0);
            return strip(s);
        end if;
    end if;
    if Detail gt 0 then printf "Decomposing space %o:%o:%o of dimension %o...", N,k,o,dNS; t:=Cputime(); end if;
    S := NewformDecomposition(NS);
    if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
    D := [QDimension(S[i]): i in [1..#S]];
    if Detail gt 0 then printf "dims = %o\n", sprint(D); end if;
    if DimensionsOnly then
        if Timings then
            return strip(Sprintf("%o:%o:%o:%o:%o:", N, k, o, Cputime()-start, Sort(D)));
        end if;
        return strip(Sprintf("%o:%o:%o:%o:", N, k, o, Sort(D)));
    end if;
    if DegreeBound eq 0 then DegreeBound := Max(D); end if;
    if #TraceHint gt 0 then
        if Detail gt 0 then printf "Computing labels for forms in space %o:%o:%o using TraceHint...",N,k,o; t:=Cputime(); end if;
        assert Multiset([t[1]:t in TraceHint]) eq Multiset(D);
        // For forms with dimension in (1,DegreeBound], Compute trace bound sufficient to distinguish forms (if dimension is unique, this will be 1)
        M :=[1:d in D];
        for i := 1 to #M do
            if ComputeEigenvalues and D[i] gt 1 and D[i] le DegreeBound then
                M[i] := n;
            else
                m := n; for j:=1 to n do if #{t[1..j]:t in TraceHint|t[1] eq D[i]} eq #[t:t in TraceHint|t[1] eq D[i]] then m:=j; break; end if; end for;
                M[i] := m;
            end if;
        end for;
        F := [* M[i] gt 1 select Eigenform(S[i],M[i]+1) else 0 : i in [1..#S] *];
        T := []; DT := [t[1] : t in TraceHint];
        for i := 1 to #S do
            if M[i] eq 1 then
                ii := Index(DT,D[i]);
                T[i] := <[Integers()|TraceHint[ii][j] : j in [1..n]],i>;
            else
                // This is as fast (or faster) than using qExpnasion, SystemOfEigenvalues, or taking traces of HeckeOperator
                T[i] := <[Integers()|Parent(a) eq Rationals() select a else AbsoluteTrace(a) where a:=Coefficient(F[i],j) : j in [1..n]],i>;
            end if;
        end for;
        T := Sort(T);
        if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
    else
        if Detail gt 0 then printf "Computing %o traces for space %o:%o:%o...", n, N,k,o; t:=Cputime(); end if;
        // This does not compute the constant term for the Eisenstein series with one of the character being trivial
        // Is there an easy way to figure out which spaces are these?
        F := [*Eigenform(S[i],n+1):i in [1..#S]*];
        T := Sort([<[Integers()|Parent(a) eq Rationals() select a else AbsoluteTrace(a) where a:=Coefficient(F[i],j) : j in [1..n]],i> : i in [1..#F]]);
        if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
    end if;
    D := [D[T[i][2]] : i in [1..#T]];  S := [S[T[i][2]] : i in [1..#T]];  F := [*F[T[i][2]] : i in [1.. #T]*];
    T := [T[i][1] : i in [1..#T]];
    assert #Set(T) eq #T;
    if Detail gt 1 then printf "Lex sorted traces = %o\n", sprint(T); end if; 
    if Detail gt 0 and Order(chi) eq 1 then printf "Computing Atkin-Lehner signs for space %o:%o:%o...", N,k,o; t:=Cputime(); end if;
    AL := Order(chi) eq 1 select [[<p,ExactQuotient(Trace(AtkinLehnerOperator(S[i],p)),D[i])>:p in PrimeDivisors(N)]:i in [1..#S]] else [];
    if Detail gt 0 and Order(chi) eq 1 then printf "took %o secs.\n", Cputime()-t; end if;
    // In order to get rigorous results for twists we need to compute Hecke cutters (regardless of dimensions!)
    if ComputeTwists or (ComputeCutters and #[d:d in D|d le DegreeBound] gt 0) then
        HC:=[[]:d in D];   
        if #D gt 1 then
            if Detail gt 0 then printf "Computing Hecke cutters for space %o:%o:%o...",N,k,o; t:=Cputime(); end if;
            HC := ComputeHeckeCutters(S);
            if Detail gt 0 then printf "cutter length %o, took %o secs\n", #HC[1], Cputime()-t; end if;
        end if;
    else
        HC := [];
    end if;
    // Deal with 1-dimensional forms first
    HF := [[0,1]:d in D|d eq 1];  d1 := #HF;
    pra:=[1:i in [1..#HF]];
    u := UnitGenerators(chi);
    X := [<u,[Eltseq(v):v in ValuesOnUnitGenerators(chi)]>:i in [1..#HF]];
    cm := []; it := [];
    tb := SturmBound(N,k);
    if #HC gt 0 and #HC[1] gt 0 and HC[1][#HC[1]][1] gt tb then
        if ComputeTwists then printf "Increasing twist bound to %o past Sturm bound %o to hit all Hecke cutters\n", HC[#HC][1][1], tb; end if;
        tb := HC[#HC][1][1];
    end if;
    if ComputeTwists then
        for i:=1 to #HF do
            if Detail gt 0 then printf "Computing self twists for form %o:%o:%o:%o of dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
            cm[i] := <1,SelfTwists(T[i],S[i]:pBound:=tb)>;
            if Detail gt 0 then printf "found self twists %o in %o secs\n", cm[i], Cputime()-t; end if;
            if Detail gt 0 then printf "Computing inner twists for form %o:%o:%o:%o of dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
            it[i] := InnerTwistData(T[i],chi,k:CharTable:=CharTable);
            if Detail gt 0 then printf "found inner twists %o in %o secs\n", it[i], Cputime()-t; end if;
        end for;
    end if;
    // Now deal with forms of dimension 2 to DegreeBound
    E := [];
    if ComputeEigenvalues then
        R<x> := PolynomialRing(Rationals());
        for i:=d1+1 to #F do
            if D[i] gt DegreeBound then break; end if;
            if Detail gt 0 then printf "Computing %o exact Hecke eigenvalues form %o:%o:%o:%o of dimension %o...",n,N,k,o,i,D[i]; t:=Cputime(); end if;
            K := AbsoluteField(BaseRing(Parent(F[i])));
            f,b,a,c,d,pr,m := OptimizedOrderBasis(Eltseq(MinimalPolynomial(K.1)),[Eltseq(K!Coefficient(F[i],j)) : j in [1..n]]:Verbose:=Detail gt 0);
            if ComputeCharacterValues then
                if Detail gt 0 then printf "Computing character values in Hecke field for form %o:%o:%o:%o of dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
                aa := NFSeq(f,b,a);
                v := #u gt 0 select [Eltseq(z):z in EmbeddedCharacterValuesOnUnitGenerators(chi,k,aa:Detail:=Detail)] else [];
                w := #u gt 0 select [Eltseq(r):r in Rows(Matrix(Rationals(),v)*Matrix(Rationals(),b)^-1)] else [];
                if Detail gt 0 then printf "Computed character values in %o secs\n", Cputime()-t; end if;                
            else
                if Detail gt 0 then printf "Computed Hecke eigenvalues but not computing character values in Hecke field."; end if;  
                v := []; w := [];
            end if;
            Append(~E,<f,b,c,<d,d eq 0 select Factorization(1) else Factorization(d)>,a,<u,w>,m>);  Append(~pra,pr select 1 else 0);  Append(~HF,f);  Append(~X,<u,v>);
            if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
            if ComputeTwists then 
                if Detail gt 0 then printf "Computing self twists for form %o:%o:%o:%o of dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
                cm[i] := <1,SelfTwists(aa,S[i]:pBound:=tb)>;
                if Detail gt 0 then printf "found self twists %o in %o secs\n", cm[i], Cputime()-t; end if;
                if Detail gt 0 then printf "Computing inner twists for form %o:%o:%o:%o of dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
                assert #aa ge tb; // guarantee rigorous result
                it[i] := InnerTwistData(aa,chi,k:CharTable:=CharTable);
                if Detail gt 0 then printf "found inner twists %o in %o secs\n", it[i], Cputime()-t; end if;
            end if;
        end for;
    else
        if ComputeTwists then
            for i:=d1+1 to #F do
                a := [Coefficient(F[i],j):j in [1..n]];
                if Detail gt 0 then printf "Computing self twists for form %o:%o:%o:%o of dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
                cm[i] := <1,SelfTwists(a,S[i]:pBound:=tb)>;
                if Detail gt 0 then printf "found self twists %o in %o secs\n", cm[i], Cputime()-t; end if;
                if Detail gt 0 then printf "Computing inner twists for form %o:%o:%o:%o of dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
                assert #a ge tb; // guarantee rigorous result
                it[i] := InnerTwistData(a,chi,k:CharTable:=CharTable);
                if Detail gt 0 then printf "found inner twists %o in %o secs\n", it[i], Cputime()-t; end if;
            end for;
        end if;
    end if;
    if ComputeTwists then
        // Compute self twist data for remaining forms
        for i := #cm+1 to #F do
            if Detail gt 0 then printf "Computing self twists for form %o:%o:%o:%o of dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
            if F[i] ne 0 then
                a := [Coefficient(F[i],j):j in [1..n]];
                cm[i] := <1,SelfTwists(a,S[i]:pBound:=tb)>;
            else
                cm[i] := <1,SelfTwists([],S[i]:TraceHint:=TraceHint[i],pBound:=tb)>;
            end if;
            if Detail gt 0 then printf "found self twists %o in %o secs\n", cm[i], Cputime()-t; end if;
        end for;
    end if;
    if ComputeTraceStats then
        Z := []; M := []; H:=[];
        for i:=1 to #D do
            if D[i] gt DegreeBound then break; end if;
            if Detail gt 0 then printf "Computing trace stats for form %o:%o:%o:%o with dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
            z,m,h := TraceStats([Integers()!Trace(Trace(a)):a in SystemOfEigenvalues(S[i],8192)], k-1);
            z := Sprintf("%.3o",z lt 0.001 select 0 else z); m:=[Sprintf("%.3o",x lt 0.001 select 0 else x):x in m]; // trim precision (still ridiculous)
            Append(~Z,z); Append(~M,m); Append(~H,h);
            if Detail gt 0 then printf "took %o secs\n", Cputime()-t; printf "trace hash: %o\n", h; end if;
        end for;
    end if;
    eap := [];
    if EmbeddingPrecision gt 0 then
        P := PrimesInInterval(1,n);
        m := d1+#E;
        for i:=m+1 to #F do
            if Detail gt 0 then printf "Computing real/complex a_p for embeddings of form %o:%o:%o:%o of dimension %o...",N,k,o,i,D[i]; t:=Cputime(); end if;
            f := F[i];
            K := CoefficientRing(Parent(f)); dK := Degree(K); dB := Degree(BaseField(K));
            // list of lists of complex aps lex-ordered by [ii,jj] where ii indexes embeddings of Q(chi) and jj indexes relative embeddings of K
            if dB eq 1 then
                A := [[Conjugate(Coefficient(f,p),[jj]:Prec:=EmbeddingPrecision):p in P]:jj in [1..dK]];
            else
                A := [[Conjugate(Coefficient(f,p),[ii,jj]:Prec:=EmbeddingPrecision):p in P]:jj in [1..dK], ii in [1..dB]];
            end if;
            L := SortEmbeddings(f,chi,EmbeddingPrecision);
            // output embeddings in order (this implicitly labels them l.m where l is a Conrey index and m is an embedding index (lex order of a_n ordered by re,im)
            if o eq 1 then
                // for trivial character we just list real parts
                Append(~eap,[[Real(c):c in A[j]]:j in L]);
            else
                Append(~eap,[[[Real(c),Imaginary(c)]:c in A[j]]:j in L]);
            end if;
            // Sanity check that embedding data sums to traces for small p (we might not have enough precision for large p)
            CC := ComplexField(EmbeddingPrecision);
            assert [T[i][p]:p in P|p le 200] eq [Round(&+[CC!eap[#eap][j][jj]:j in [1..#eap[#eap]]]):jj in [1..#P]|P[jj] le 200];
            if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
        end for;
    end if;
    if Timings then 
        s := Sprintf("%o:%o:%o:%o:%o", N, k, o, Cputime()-start, D);
    else
        s := Sprintf("%o:%o:%o:%o", N, k, o, D);
    end if;
    s cat:= Sprintf(":%o:%o:%o:%o",T,AL,HF,HC);
    if ComputeEigenvalues then s cat:= Sprintf(":%o:%o:%o:%o",E,cm,it,pra); else s cat:= ":[]:[]:[]:[]"; end if;
    if ComputeTraceStats then s cat:= Sprintf(":%o:%o:%o", Z, M, H); else s cat:= ":[]:[]:[]"; end if;
    s cat:= Sprintf(":%o",X);
    if ComputeEigenvalues then s cat:= Sprintf(":%o",[IsSelfDual(chi,D[i],T[i],#HF ge i select HF[i] else [],S[i]) select 1 else 0:i in [1..#D]]); else s cat:= ":[]"; end if;
    s cat:= Sprintf(":%o",eap);
    if ReturnDecomposition then return strip(s),S; else return strip(s); end if;
end function;

function DimensionData(chi, k, o)
    N := Modulus(chi);
    dS := QDimensionCuspForms(chi,k);
    dNS := QDimensionNewCuspForms(chi,k);
    dE := QDimensionEisensteinForms(chi,k);
    dNE := QDimensionNewEisensteinForms(chi,k);
    s := Sprintf("%o:%o:%o:%o:%o:%o:%o", N, k, o, dS, dNS, dE, dNE);
    return s;
end function;

// Decompose spaces S_k(N,chi)^new into Galois stable subspaces for B0 < k^2*N <= B and k > 1.
procedure DecomposeSpace(outfile,B :TodoFile:="",B0:=0,Quiet:=false,DimensionsOnly:=false,Coeffs:=1000,DegBound:=20,Eigenvalues:=true,Cutters:=true,Twists:=true,
                          TrivialCharOnly:=false,TraceStats:=false,Precision:=0,Jobs:=1,JobId:=0,Timings:=true)
    if Jobs ne 1 then outfile cat:= Sprintf("_%o",JobId); end if;
    if Jobs ne 1 then dimfile cat:= Sprintf("_%o",JobId); end if;
    if TodoFile ne "" then
        TodoList := AssociativeArray();
        for s in Split(Read(TodoFile),"\n") do
            r := Split(s,":");
            TodoList[[StringToInteger(r[1]),StringToInteger(r[2]),StringToInteger(r[3])]] := #r ge 5 select eval(r[5]) else [];
        end for;
        L := {r[1]:r in Keys(TodoList)};
        printf "Loaded todo list with %o items from %o\n", #Keys(TodoList), TodoFile;
    else
        L := {};
    end if;
    st := Cputime();
    n := 0; cnt:=0;
    fp := Open(outfile,"w");
    // dim_fp := Open(dimfile, "w");
    A := AssociativeArray();
    for N:=1 to Floor(B/4) do
        if #L gt 0 and not N in L then continue; end if;
        if not TrivialCharOnly then
            if not Quiet then printf "Constructing character orbit table for modulus %o...", N; t:=Cputime(); end if;
            G, T := CharacterOrbitReps(N:RepTable); A[N] := <G,T>;
            if not Quiet then printf "took %o secs\n",Cputime()-t; end if;
        end if;
        for k := Max(2,Floor(Sqrt(B0/N))+1) to Floor(Sqrt(B/N)) do
            m := TrivialCharOnly select 1 else #A[N][1];
            for o in [1..m] do
                if #L gt 0 and not IsDefined(TodoList,[N,k,o]) then continue; end if;
                hint := #L gt 0 select TodoList[[N,k,o]] else [];
                chi := o gt 1 select A[N][1][o] else DirichletGroup(N)!1;
                n +:= 1;
                if ((n-JobId) mod Jobs) eq 0 then
                    if DimensionsOnly then
                        str := NewspaceData(chi,k,o:DimensionsOnly:=true,Detail:=Quiet select 0 else 1,Timings:=Timings);
                    else
                        // Note that we need character orbit tables even when TrivialCharOnly is set because we may have twists by non-trivial characters (e.g. for CM forms)
                        if not Quiet then printf "Constructing character orbit table for divisors of modulus %o...", N; t:=Cputime(); end if;
                        for M in Divisors(N) do if not IsDefined(A,M) then G, T := CharacterOrbitReps(M:RepTable); A[M] := <G,T>; end if; end for;
                        if not Quiet then printf "took %o secs\n",Cputime()-t; end if;
                        if not Quiet then printf "\nProcessing space %o:%o:%o with Coeffs=%o, DegBound=%o\n", N,k,o, Coeffs, DegBound; t:=Cputime(); end if;
                        str := NewspaceData(chi,k,o:CharTable:=A,TraceHint:=hint,NumberOfCoefficients:=Coeffs,ComputeEigenvalues:=Eigenvalues,ComputeCutters:=Cutters,EmbeddingPrecision:=Precision,ComputeTwists:=Twists,ComputeTraceStats:=TraceStats,DegreeBound:=DegBound,Detail:=Quiet select 0 else 1,Timings:=Timings);
                        if not Quiet then printf "Total time for space %o:%o:%o was %os\n\n", N,k,o,Cputime()-t; end if;
                    end if;
                    Puts(fp,str);
                    Flush(fp);
                    // Puts(dim_fp, DimensionData(chi, k, o));
                    // Flush(dim_fp);
                    cnt +:= 1;
                end if;
            end for;
        end for;
    end for;
    printf "Wrote %o records to %o using %os of CPU time.\n", cnt, outfile, Cputime()-st;
end procedure;