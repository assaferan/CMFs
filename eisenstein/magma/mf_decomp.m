// !! TODO - work out why Magma still has bugs in constructing Eisenstein series

import "../../magma/mf.m" : ComputeHeckeCutters, InnerTwistData, SortEmbeddings;

function GaussSum(chi)
    N := Conductor(chi);
    if N eq 1 then
	return 1;
    end if;
    Z_N_star := [i : i in [0..N-1] | GCD(i,N) eq 1];
    K<zeta> := CyclotomicField(N);
    return &+[chi(n)*zeta^n : n in Z_N_star];
end function;

function QLvalueNegative(chi, k)
    N := Modulus(chi);
    f := Conductor(chi);
    chi_0 := DirichletGroupFull(f)!chi;
    if IsTrivial(chi) then
	    prim_val := - BernoulliNumber(1-k)/(1-k);
    else
	    prim_val := - BernoulliNumber(1-k, chi_0) / (1-k);
    end if;
    fac_N := Factorization(N);
    primes := [p[1] : p in fac_N | f mod p[1] ne 0];
    if IsEmpty(primes) then
      prod := 1;
    else
      prod := &*[1 - (chi_0^(-1))(p) / p^k : p in primes];
    end if;
    return prim_val * prod * (N/f)^(1-k);
end function;

// Exact version giving the values divided by C_k = (-2pi i)^k / N^k (k-1)!
function QLvalue(chi, k)
    if (k lt 1) then 
        // use functional equation
        return QLvalueNegative(chi,k);
    end if; 
    N := Modulus(chi);
    f := Conductor(chi);
    chi_0 := DirichletGroupFull(f)!chi;
    if IsTrivial(chi) then
	    prim_val := - BernoulliNumber(k)/k;
    else
	    prim_val := - GaussSum(chi_0) * BernoulliNumber(k, chi_0^(-1)) / k;
    end if;
    fac_N := Factorization(N);
    primes := [p[1] : p in fac_N | f mod p[1] ne 0];
    if IsEmpty(primes) then
      prod := 1;
    else
      prod := &*[1 - chi_0(p) / p^k : p in primes];
    end if;
    return prim_val * prod * (N/f)^k;
end function;

function EisensteinSeries(psi, phi, k : Bound := 1000)
    assert k ge 3; // For now we focus on k >= 3, no special cases
    assert IsPrimitive(psi);
    assert IsPrimitive(phi);
    u := Conductor(psi);
    v := Conductor(phi);
    chi := psi*phi;
    assert chi(-1) eq (-1)^k;
    delta := IsTrivial(psi) select QLvalue(phi, 1-k) else 0;
    n_psi := Conductor(Codomain(MinimalBaseRingCharacter(psi)));
    n_phi := Conductor(Codomain(MinimalBaseRingCharacter(phi)));
    K := CyclotomicField(LCM(n_psi, n_phi));
    Kq<q> := PowerSeriesRing(K);
    // coeffs := [K | ]; 
    coeffs := [K | delta/2]; // this is the constant term, we write it separately 
    // E := Kq!(delta/2);
    for n in [1..Bound] do
        // This is half of the convention in [DS] as we normalize to have the coefficient of q be 1
        an := &+[psi(n div m) * phi(m) * m^(k-1) : m in Divisors(n)];
        Append(~coeffs, an);
        // E +:= an*q^n;
    end for;
    E := Kq!coeffs;
    return E;
end function;

function EisensteinAdmissibleCharacterPairs(chi, k : IdentifyConjugates := true, NewformsOnly := true)
  
  N := Modulus(chi);
  
  H := FullDirichletGroup(N);

  if (chi(-1) ne (-1)^k) then
      return [* *];
  end if;
  
  check_chi := func<eta, psi | (eta * psi eq chi)>;
  
  // By default we only produce pairs of characters corresponding to newforms,
  // i.e. such that Cond(eta) * Cond(psi) = N, but for backwards compatibility
  // and testing it is useful to be able to produce the pairs corresponding to
  // all the Eisenstein series.
  check_n := func<eta, psi, exact |\
    (exact) select (N eq Conductor(eta) * Conductor(psi)) else (N subset Conductor(eta) * Conductor(psi))>;

  pairs := &join{{<eta, psi> : psi in Elements(H) |\
    check_chi(eta, psi) and check_n(eta, psi, NewformsOnly)} : eta in Elements(H)};

  mth_power := func<pair, m | <pair[1]^m, pair[2]^m>>;
  n := Exponent(AbelianGroup(H));
  coprime_to_n := [m : m in [1 .. n] | IsCoprime(m, n)];

  if IdentifyConjugates then
    pairs_up_to_galois := {};
    while #pairs gt 0 do
      pair := Representative(pairs);
      triple := <pair[1], pair[2], 0>;
      for m in coprime_to_n do
        other := mth_power(pair, m);
        if other in pairs then 
            Exclude(~pairs, other);
            triple[3] +:= 1;
        end if;
      end for;
      Include(~pairs_up_to_galois, triple);
      if k eq 1 then
        Exclude(~pairs, <pair[2], pair[1]>);
      end if;
    end while;
    pairs := pairs_up_to_galois;
  end if;

  APC := func<pair | <AssociatedPrimitiveCharacter(pair[1]), AssociatedPrimitiveCharacter(pair[2]), pair[3]>>;
  pairs := [* APC(pair) : pair in pairs *];
  return pairs;
end function;

function QDimension(triples)
    return &+[triple[3] : triple in triples];
end function;

function decompose_character(chi, Q)
    dec := Decomposition(chi);
    triv := DirichletGroup(1)!1;
    chi_Q := triv;
    chi_R := triv;
    for psi in dec do
	if Q mod Modulus(psi) eq 0 then
	    chi_Q *:= psi;
	else
	    chi_R *:= psi;
	end if;
    end for;
    return chi_Q, chi_R;
end function;

function AL_action(chi1, chi2, Q)
    assert IsPrimitive(chi1) and IsPrimitive(chi2);
    chi := chi1 * chi2^(-1);
    N := Modulus(chi);
    N1 := Modulus(chi1);
    N2 := Modulus(chi2);
    Q1 := GCD(N1, Q);
    Q2 := GCD(N2, Q);
    assert N mod Q eq 0;
    R := N div Q;
    assert GCD(Q,R) eq 1;
    R1 := GCD(N1, R);
    R2 := GCD(N2, R);
    chi_Q, chi_R := decompose_character(chi, Q);
    chi_1_Q, chi_1_R := decompose_character(chi1, Q);
    chi_2_Q, chi_2_R := decompose_character(chi2, Q);
    chi_1_prime := chi_2_Q * chi_1_R;
    chi_2_prime := chi_1_Q * chi_2_R;
    scalar := Evaluate(chi_2_Q, -1) * Evaluate(chi_Q^(-1), R2) * Evaluate(chi_R^(-1), Q2);
    sqrt_arg := Q1 / Q2;
    return sqrt_arg, scalar, chi_1_prime, chi_2_prime;
end function;

function trace_AL_action(chi1, chi2, Q)
    _, scalar, chi1_prime, chi2_prime := AL_action(chi1, chi2, Q);
    if (chi1 eq chi1_prime) and (chi2 eq chi2_prime) then
        return scalar;
    end if;
    return 0;
end function;

/*
Format of data is N:k:i:t:D:T:A:F:C:E:cm:tw:pra:zr:mm:h:X:sd:eap
The data depends on a degree bound (determines forms with exact eigenvalue data, a coefficient count (number of a_n to compute), and a complex precision (for forms without exact eigenvalue data)

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

function NewspaceData (chi, k, o: CharTable:=AssociativeArray(), DimensionsOnly:=false, ComputeEigenvalues:=false, ComputeTwists:=false, ComputeTraceStats:=false,
                       NumberOfCoefficients:=0, DegreeBound:=0, EmbeddingPrecision:= 0, Detail:=0, ReturnDecomposition:=false, ComputeCutters:=false, ComputeCharacterValues:=true)
    start := Cputime();
    if o eq 0 then o := CharacterOrbit(chi); end if;
    N := Modulus(chi);
    dNS := QDimensionNewEisensteinForms(chi,k);
    if dNS eq 0 then
        if Detail gt 0 then printf "The space %o:%o:%o is empty\n",N,k,o; end if;
        s := Sprintf("%o:%o:%o:%o:[]", N, k, o, Cputime()-start);
        if not DimensionsOnly then s cat:= ":[]:[]:[]:[]:[]:[]:[]:[]:[]:[]:[]:[]:[]:[]"; end if;
        return strip(s);
    end if;
    if Detail gt 0 then printf "Constructing admissible pairs for space %o:%o:%o...", N,k,o; t:=Cputime(); end if;
    S := EisensteinAdmissibleCharacterPairs(chi, k);
    if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
    if NumberOfCoefficients eq 0 then
        if N le 1000 then NumberOfCoefficients := 1000; end if;
        if N gt 1000 and N le 4000 then NumberOfCoefficients := 2000; end if;
        if N gt 4000 and N le 10000 then NumberOfCoefficients := 3000; end if;
    end if;
    n := Max([SturmBound(N,k)+1,Floor(30*Sqrt(N)),NumberOfCoefficients]);
    assert QDimension(S) eq dNS;
    D := [t[3]: t in S];
    if Detail gt 0 then printf "dims = %o\n", sprint(D); end if;
    if DimensionsOnly then
        return strip(Sprintf("%o:%o:%o:%o:%o:", N, k, o, Cputime()-start, Sort(D)));
    end if;
    if DegreeBound eq 0 then DegreeBound := Max(D); end if;
    if Detail gt 0 then printf "Computing %o traces for space %o:%o:%o...", n, N,k,o; t:=Cputime(); end if;
    F := [* EisensteinSeries(t[1],t[2],k : Bound := n): t in S *];
    T := Sort([<[Integers()|Parent(a) eq Rationals() select a else AbsoluteTrace(a) where a:=Coefficient(F[i],j) : j in [1..n]],i> : i in [1..#F]]);
    if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
    
    D := [D[T[i][2]] : i in [1..#T]];  S := [* S[T[i][2]] : i in [1..#T] *];  F := [*F[T[i][2]] : i in [1.. #T]*];
    T := [T[i][1] : i in [1..#T]];
    assert #Set(T) eq #T;
    if Detail gt 1 then printf "Lex sorted traces = %o\n", sprint(T); end if; 
    if Detail gt 0 and Order(chi) eq 1 then printf "Computing Atkin-Lehner signs for space %o:%o:%o...", N,k,o; t:=Cputime(); end if;
    AL := Order(chi) eq 1 select [[<p,ExactQuotient(trace_AL_action(S[i][1],S[i][2],p^Valuation(N,p)),D[i])>:p in PrimeDivisors(N)]:i in [1..#S]] else [];
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
                cm[i] := <1,SelfTwists([],S[i]:TraceHint:=0,pBound:=tb)>;
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
    s := Sprintf("%o:%o:%o:%o:%o", N, k, o, Cputime()-start, D);
    s cat:= Sprintf(":%o:%o:%o:%o",T,AL,HF,HC);
    if ComputeEigenvalues then s cat:= Sprintf(":%o:%o:%o:%o",E,cm,it,pra); else s cat:= ":[]:[]:[]:[]"; end if;
    if ComputeTraceStats then s cat:= Sprintf(":%o:%o:%o", Z, M, H); else s cat:= ":[]:[]:[]"; end if;
    s cat:= Sprintf(":%o",X);
    if ComputeEigenvalues then s cat:= Sprintf(":%o",[IsSelfDual(chi,D[i],T[i],#HF ge i select HF[i] else [],S[i]) select 1 else 0:i in [1..#D]]); else s cat:= ":[]"; end if;
    s cat:= Sprintf(":%o",eap);
    if ReturnDecomposition then return strip(s),S; else return strip(s); end if;
end function;

// Still working on this one - 
// Decompose spaces S_k(N,chi)^new into Galois stable subspaces for B0 < k^2*N <= B and k > 1.
procedure DecomposeSpaces (outfile,B:TodoFile:="",B0:=0,Quiet:=false,DimensionsOnly:=false,Coeffs:=1000,Eigenvalues:=true,Cutters:=true,Twists:=true,
                           NBound:=0,kBound:=0,DegBound:=20,TrivialCharOnly:=false,TraceStats:=false,Precision:=0,Jobs:=1,JobId:=0)
    if Jobs ne 1 then outfile cat:= Sprintf("_%o",JobId); end if;
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
    A := AssociativeArray();
    for N:=1 to Floor(B/4) do
        if NBound gt 0 and N gt NBound then break; end if;
        if #L gt 0 and not N in L then continue; end if;
        if not TrivialCharOnly then
            if not Quiet then printf "Constructing character orbit table for modulus %o...", N; t:=Cputime(); end if;
            G, T := CharacterOrbitReps(N:RepTable); A[N] := <G,T>;
            if not Quiet then printf "took %o secs\n",Cputime()-t; end if;
        end if;
        for k := Max(2,Floor(Sqrt(B0/N))+1) to Floor(Sqrt(B/N)) do
            if kBound gt 0 and k gt kBound then break; end if;
            m := TrivialCharOnly select 1 else #A[N][1];
            for o in [1..m] do
                if #L gt 0 and not IsDefined(TodoList,[N,k,o]) then continue; end if;
                hint := #L gt 0 select TodoList[[N,k,o]] else [];
                chi := o gt 1 select A[N][1][o] else DirichletGroup(N)!1;
                n +:= 1;
                if ((n-JobId) mod Jobs) eq 0 then
                    if DimensionsOnly then
                        str := NewspaceData(chi,k,o:DimensionsOnly:=true,Detail:=Quiet select 0 else 1);
                    else
                        // Note that we need character orbit tables even when TrivialCharOnly is set because we may have twists by non-trivial characters (e.g. for CM forms)
                        if not Quiet then printf "Constructing character orbit table for divisors of modulus %o...", N; t:=Cputime(); end if;
                        for M in Divisors(N) do if not IsDefined(A,M) then G, T := CharacterOrbitReps(M:RepTable); A[M] := <G,T>; end if; end for;
                        if not Quiet then printf "took %o secs\n",Cputime()-t; end if;
                        if not Quiet then printf "\nProcessing space %o:%o:%o with Coeffs=%o, DegBound=%o\n", N,k,o, Coeffs, DegBound; t:=Cputime(); end if;
                        Nustr := NewspaceData(chi,k,o:CharTable:=A,TraceHint:=hint,NumberOfCoefficients:=Coeffs,ComputeEigenvalues:=Eigenvalues,ComputeCutters:=Cutters,EmbeddingPrecision:=Precision,ComputeTwists:=Twists,ComputeTraceStats:=TraceStats,DegreeBound:=DegBound,Detail:=Quiet select 0 else 1);
                        if not Quiet then printf "Total time for space %o:%o:%o was %os\n\n", N,k,o,Cputime()-t; end if;
                    end if;
                    Puts(fp,str);
                    Flush(fp);
                    cnt +:= 1;
                end if;
            end for;
        end for;
    end for;
    printf "Wrote %o records to %o using %os of CPU time.\n", cnt, outfile, Cputime()-st;
end procedure;

procedure WriteSpaceData(fp, chi, k, o: CharTable:=AssociativeArray(), ComputeEigenvalues:=true, 
                                        NumberOfCoefficients:=0, DegreeBound:=0, Detail:=0, Sep := ",", MinTrace := 0)
    start := Cputime();
    if o eq 0 then o := CharacterOrbit(chi); end if;
    N := Modulus(chi);
    dNS := QDimensionNewEisensteinForms(chi,k);
    if dNS eq 0 then
        if Detail gt 0 then printf "The space %o:%o:%o is empty\n",N,k,o; end if;
        return;
    end if;
    if Detail gt 0 then printf "Constructing admissible pairs for space %o:%o:%o...", N,k,o; t:=Cputime(); end if;
    S := EisensteinAdmissibleCharacterPairs(chi, k);
    if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
    /*
    if NumberOfCoefficients eq 0 then
        if N le 1000 then NumberOfCoefficients := 1000; end if;
        if N gt 1000 and N le 4000 then NumberOfCoefficients := 2000; end if;
        if N gt 4000 and N le 10000 then NumberOfCoefficients := 3000; end if;
    end if;
    n := Max([SturmBound(N,k)+1,Floor(30*Sqrt(N)),NumberOfCoefficients]);
    */
    n := NumberOfCoefficients;
    // assert QDimension(S) eq dNS;
    D := [t[3]: t in S];
    if Detail gt 0 then printf "dims = %o\n", sprint(D); end if;
    if DegreeBound eq 0 then DegreeBound := Max(D); end if;
    if Detail gt 0 then printf "Computing %o traces for space %o:%o:%o...", n, N,k,o; t:=Cputime(); end if;
    F := [* EisensteinSeries(t[1],t[2],k : Bound := n): t in S *];
    T := Sort([<[Integers()|Parent(a) eq Rationals() select a else AbsoluteTrace(a) where a:=Coefficient(F[i],j) : j in [1..n]],i> : i in [1..#F]]);
    if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
    
    D := [D[T[i][2]] : i in [1..#T]];  S := [* S[T[i][2]] : i in [1..#T] *];  F := [*F[T[i][2]] : i in [1.. #T]*];
    T := [T[i][1] : i in [1..#T]];
    assert #Set(T) eq #T;
    if Detail gt 1 then printf "Lex sorted traces = %o\n", sprint(T); end if; 
    // !! TODO - figure out where I can find first !!
    GalOrb := Sprintf("%o.%o", N, Base26Encode(o-1));
    prefix := Join([Sprintf("%o", fld) : fld in [* GalOrb, ConreyIndex(chi), k, NumberOfCoefficients *]], Sep);
    R<x> := PolynomialRing(Rationals());
    space_strs := [];
    for i:=1 to #F do
        if D[i] gt DegreeBound then break; end if;
        if Detail gt 0 then printf "Computing %o exact Hecke eigenvalues form %o:%o:%o:%o of dimension %o...",n,N,k,o,i,D[i]; t:=Cputime(); end if;
        K := AbsoluteField(BaseRing(Parent(F[i])));
        // f,b,a,c,d,pr,m := OptimizedOrderBasis(Eltseq(MinimalPolynomial(K.1)),[Eltseq(K!Coefficient(F[i],j)) : j in [1..n]]:Verbose:=Detail gt 0);
        assert IsCyclotomic(K);
        cyc := CyclotomicOrder(K);
        an := [K!Coefficient(F[i],j) : j in [0..n]];
        // an_coeffs := [Eltseq(a) : a in an];
        an_str := "[" cat Join(["[" cat Join([Sprintf("%o", x) : x in Eltseq(a)], ",") cat "]" : a in an], ",") cat "]";
        traces := [Trace(a) : a in an[MinTrace..#an]]; // !! TODO - note that this has already been computed for n > 0
        trace_str := "[" cat Join([Sprintf("%o", t) : t in traces],",") cat "]";
        label := NewformEisensteinLabel(N,k,o,i);
        suffix := Join([Sprintf("%o", fld) : fld in [* cyc,an_str,trace_str,label *]], Sep);
        str := Join([prefix, suffix], Sep);
        Puts(fp,str);
        Flush(fp);
        if Detail gt 0 then printf "took %o secs\n", Cputime()-t; end if;
    end for;
    return;
end procedure;

// Note - in terms of data generation, it will be easier to loop on k last
procedure WriteDataFile(outfile, k_min, k_max, N_min, N_max, lim : Quiet := false, Jobs := 1, JobId := 0, Eigenvalues := true, 
                                                                   DegBound := 0, MinTrace := 0)
    if Jobs ne 1 then outfile cat:= Sprintf("_%o",JobId); end if;
    st := Cputime();
    n := 0; cnt:=0;
    fp := Open(outfile,"w");  
    A := AssociativeArray();  
    for k in [k_min..k_max] do
        for N in [N_min..N_max] do
            if not Quiet then printf "Constructing character orbit table for modulus %o...", N; t:=Cputime(); end if;
            G, T := CharacterOrbitReps(N:RepTable); A[N] := <G,T>;
            if not Quiet then printf "took %o secs\n",Cputime()-t; end if;
            m := #A[N][1];
            for o in [1..m] do
                chi := o gt 1 select A[N][1][o] else DirichletGroup(N)!1;
                n +:= 1;
                if ((n-JobId) mod Jobs) eq 0 then
                    if not Quiet then printf "Constructing character orbit table for divisors of modulus %o...", N; t:=Cputime(); end if;
                    for M in Divisors(N) do if not IsDefined(A,M) then G, T := CharacterOrbitReps(M:RepTable); A[M] := <G,T>; end if; end for;
                    if not Quiet then printf "took %o secs\n",Cputime()-t; end if;
                    if not Quiet then printf "\nProcessing space %o:%o:%o with Coeffs=%o, DegBound=%o\n", N,k,o, lim, DegBound; t:=Cputime(); end if;
                    WriteSpaceData(fp,chi,k,o: CharTable:=A,NumberOfCoefficients:=lim,ComputeEigenvalues:=Eigenvalues,
                                               DegreeBound := DegBound, Detail:=Quiet select 0 else 1, MinTrace := MinTrace);
                    if not Quiet then printf "Total time for space %o:%o:%o was %os\n\n", N,k,o,Cputime()-t; end if;
                    cnt +:= 1;
                end if;
            end for;
        end for;
    end for;
    delete fp;
    printf "Wrote %o records to %o using %os of CPU time.\n", cnt, outfile, Cputime()-st;
end procedure;