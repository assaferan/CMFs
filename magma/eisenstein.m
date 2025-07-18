import "Lvalue.m" : QLvalue;

function EisensteinSeries(psi, phi, k : Bound := 1000)
    // assert k ge 3; // For now we focus on k >= 3, no special cases
    assert IsPrimitive(psi);
    assert IsPrimitive(phi);
    u := Conductor(psi);
    v := Conductor(phi);
    chi := psi*phi;
    assert chi(-1) eq (-1)^k;
    delta := IsTrivial(psi) select QLvalue(phi, 1-k) else 0;
    n_psi := Conductor(Codomain(MinimalBaseRingCharacter(psi)));
    n_phi := Conductor(Codomain(MinimalBaseRingCharacter(phi)));
    m := LCM(n_psi, n_phi);
    if IsEven(m) and IsOdd(m div 2) then m div:= 2; end if;
    K := (m eq 1) select Rationals() else CyclotomicField(m);
    // coeffs := [K | ]; 
    coeffs := [K | delta/2]; // this is the constant term, we write it separately 
    // E := Kq!(delta/2);
    for n in [1..Bound] do
        // This is half of the convention in [DS] as we normalize to have the coefficient of q be 1
        an := &+[psi(n div m) * phi(m) * m^(k-1) : m in Divisors(n)];
        Append(~coeffs, an);
        // E +:= an*q^n;
    end for;
    // This might be expensive. Can we know it ahead of time?
    if Type(K) ne FldRat then K := sub<K | coeffs>; end if;
    Kq<q> := PowerSeriesRing(K);
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
  /*
  pairs := &join{{<eta, psi> : psi in Elements(H) |\
    check_chi(eta, psi) and check_n(eta, psi, NewformsOnly)} : eta in Elements(H)};
  */

  // check_chi is taking very long - trying to cut down the time by not running it
  pairs := {<eta, psi> : psi in Elements(H) | check_n(eta, psi, NewformsOnly) where eta := H!(chi*psi^(-1))};

  // assert pairs eq pairs2;

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

function GaloisConjugacyNewEisensteinSeries(chi,k : Bound := 1000)
  pairs := EisensteinAdmissibleCharacterPairs(chi, k);
  eis := [* EisensteinSeries(p[1], p[2], k : Bound := Bound) : p in pairs *];
  if IsTrivial(chi) and (k eq 2) then
    N := Modulus(chi);
    if IsPrime(N) then
      chi0 := DirichletGroup(1)!1;
      E2<q> := EisensteinSeries(chi0, chi0, 2 : Bound := Bound);
      eis cat:= [* E2 - N*Evaluate(E2, q^N) *];
    end if;
  end if;
  return eis;
end function;

function NumOrbitsUpTo(B)
  cnt := 0;
  for N:=1 to Floor(B/4) do
    conrey_labels := ConreyCharacterOrbitReps(N);
    m := #conrey_labels;
    for o in [1..m] do
      chi := DirichletCharacter(conrey_labels[o]);
      max_k := Floor(Sqrt(B/N));
      num_k := IsEven(chi) select max_k div 2 else (max_k-1) div 2;
      k0 := IsEven(chi) select 4 else 3;
      cnt +:= num_k * #EisensteinAdmissibleCharacterPairs(chi, k0);
      if (o eq 1) and IsPrime(N) then cnt +:= 1; end if; // for weight 2
    end for;
  end for;
  return cnt;
end function;