
# include("inv_sqrt_lambda.jl")
include("isometryCorboz.jl")
include("renormalise.jl")

function evolution_ctm(T1,T2,T3,C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da,tensa,chi::Int,iso::String="Corboz")


    C1ab = C1cf*T1.E*T3.D*(tensa.A*tensa.B)
    C2cd = C2eb*T2.A*T1.F*(tensa.C*tensa.D)
    C3fe = C3da*T3.C*T2.B*(tensa.F*tensa.E)

    ind12 = commoninds(C1ab,C2cd)
    ind23 = commoninds(C2cd,C3fe)
    ind31 = commoninds(C3fe,C1ab)

    prime!(C1ab,ind12)
    rho1 = (C3fe*C1ab)*C2cd
    noprime!(C1ab)

    prime!(C2cd,ind23)
    rho2 = (C3fe*C1ab)*C2cd
    noprime!(C2cd)

    prime!(C3fe,ind31)
    rho3 = (C3fe*C2cd)*C1ab
    noprime!(C3fe)


    if iso == "Corboz"
        P,Pbar = isometryCorboz(C1ab,C2cd,chi)
    elseif iso == "nishino"
        u,s,v = svd(rho1, ind12', maxdim = chi) # rho1 = C1'C2C3
        noprime!(u)
        v = copy(u)
        P = dag(u)
        Pbar = u
    end

    C1ab = C1ab*P
    C2cd = C2cd*Pbar
    T1.C = T1.F*tensa.C*Pbar
    T1.B = T1.E*tensa.B*P

    if iso == "Corboz"
        P,Pbar = isometryCorboz(C2cd,C3fe,chi)
    elseif iso == "nishino"
        u,s,v = svd(rho2, ind23', maxdim = chi) # rho2 = C1C2'C3
        noprime!(u)
        P = dag(u)
        Pbar = u
    end

    C2cd = C2cd*P
    C3fe = C3fe*Pbar
    T2.E = T2.B*tensa.E*Pbar
    T2.D = T2.A*tensa.D*P

    if iso == "Corboz"
        P,Pbar = isometryCorboz(C3fe,C1ab,chi)
    elseif iso == "nishino"
        u,s,v = svd(rho3, ind31', maxdim = chi) # rho3 = C1C2C3'
        noprime!(u)
        P = dag(u)
        Pbar = u 
    end

    C3fe = C3fe*P
    C1ab = C1ab*Pbar
    T3.F = T3.C*tensa.F*P
    T3.A = T3.D*tensa.A*Pbar


    renormalise_tensor_to_norm!(C1ab,C2cd,C3fe,T1.C,T1.B,T2.E,T2.D,T3.F,T3.A)


    ################## SECOND STEP ##################

    C1ed = C1ab*T1.C*T3.F*(tensa.E*tensa.D)
    C2af = C2cd*T1.B*T2.E*(tensa.A*tensa.F)
    C3bc = C3fe*T3.A*T2.D*(tensa.B*tensa.C)


    ind12 = commoninds(C1ed, C2af)
    ind23 = commoninds(C2af, C3bc)
    ind31 = commoninds(C3bc, C1ed)

    prime!(C1ed, ind12)
    rho1 = C3bc*C2af*C1ed
    noprime!(C1ed)

    prime!(C2af, ind23)
    rho2 = C1ed*C3bc*C2af
    noprime!(C2af)

    prime!(C3bc, ind31)
    rho3 = C1ed*C2af*C3bc
    noprime!(C3bc)


    if iso == "Corboz"
        P,Pbar = isometryCorboz(C1ed,C2af,chi)
    elseif iso == "nishino"
        u,s,v = svd(rho1, (ind12'), maxdim = chi) # rho1 = C1'*C2*C3
        noprime!(u)
        P = dag(u)
        Pbar = u
    end

    C1ed = C1ed*P
    C2af = C2af*Pbar
    T1.A = T1.B*tensa.A*Pbar
    T1.D = T1.C*tensa.D*P

    if iso == "Corboz"
        P,Pbar = isometryCorboz(C2af, C3bc,chi)
    elseif iso == "nishino"
        u,s,v = svd(rho2, ind23', maxdim = chi) # rho2 = C1*C2'*C3
        noprime!(u)
        P = dag(u)
        Pbar = u
    end

    C2af = C2af*P
    C3bc = C3bc*Pbar
    T2.C = T2.D*tensa.C*Pbar
    T2.F = T2.E*tensa.F*P

    if iso == "Corboz"
        P,Pbar = isometryCorboz(C3bc, C1ed,chi)
    elseif iso == "nishino"
        u,s,v = svd(rho3, ind31', maxdim = chi) # rho3 = C1*C2*C3'
        noprime!(u)
        P = dag(u) 
        Pbar = u 
    end

    C3bc = C3bc*P
    C1ed = C1ed*Pbar
    T3.E = T3.F*tensa.E*Pbar
    T3.B = T3.A*tensa.B*P


    renormalise_tensor_to_norm!(C1ed,C2af,C3bc,T1.A,T1.D,T2.C,T2.F,T3.E,T3.B)

    
    ################## THIRD STEP ##################

    C1cf = C1ed*T1.A*T3.B*(tensa.C*tensa.F)
    C2eb = C2af*T1.D*T2.C*(tensa.E*tensa.B)
    C3da = C3bc*T3.E*T2.F*(tensa.D*tensa.A)
    
    ind12 = commoninds(C1cf, C2eb)
    ind23 = commoninds(C2eb, C3da)
    ind31 = commoninds(C3da, C1cf)


    prime!(C1cf,ind12)
    rho1 = C2eb*C3da*C1cf
    noprime!(C1cf)

    prime!(C2eb,ind23)
    rho2 = C3da*C1cf*C2eb
    noprime!(C2eb)

    prime!(C3da,ind31)
    rho3 = C1cf*C2eb*C3da
    noprime!(C3da)


    if iso == "Corboz"
        P,Pbar = isometryCorboz(C1cf,C2eb,chi)
    elseif iso == "nishino"
        u,s,v = svd(rho1, (ind12'), maxdim = chi) # rho1 = C1'C2C3
        noprime!(u)
        P = dag(u)
        Pbar = u
    end

    C1cf = C1cf*P
    C2eb = C2eb*Pbar
    T1.E = T1.D*tensa.E*Pbar
    T1.F = T1.A*tensa.F*P

    if iso == "Corboz"
        P,Pbar = isometryCorboz(C2eb, C3da,chi)
    elseif iso == "nishino"
        u,s,v = svd(rho2, ind23', maxdim = chi) # rho2 = C1C2'C3
        noprime!(u)
        P = dag(u)
        Pbar = u
    end

    C2eb = C2eb*P
    C3da = C3da*Pbar
    T2.A = T2.F*tensa.A*Pbar
    T2.B = T2.C*tensa.B*P


    if iso == "Corboz"
        P,Pbar = isometryCorboz(C3da, C1cf,chi)
    elseif iso == "nishino"
        u,s,v = svd(rho3, ind31', maxdim = chi) # rho3 = C1C2C3'
        noprime!(u)
        P = dag(u)
        Pbar = u
    end

    C3da = C3da*P
    C1cf = C1cf*Pbar
    T3.C = T3.B*tensa.C*Pbar
    T3.D = T3.E*tensa.D*P

    renormalise_tensor_to_norm!(C1cf,C2eb,C3da,T1.E,T1.F,T2.A,T2.B,T3.C,T3.D)




    return T1,T2,T3,C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da


end