
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
