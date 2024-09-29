
include("renormalise.jl")
include("get_obs_tensors.jl")

function energy(tensa, tensA, cx, bond_dimension, physical_legs, ancilla_legs, hamil,hex1,hex2)


    # CONVERGENCE CRITERIA
    sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
    obs_criteria = ITensor(physical_legs.A', dag(physical_legs.A))
    obs = kron(id,sx,sx) + kron(id,sy,sy) + kron(id,sz,sz)
    charges = [-3, -1, -1, 1, -1, 1, 1, 3];

    for i1 = 1:8, i2 = 1:8
        if norm(charges[i1]-charges[i2]) < 1e-9
            obs_criteria[physical_legs.A'=>i1,physical_legs.A=>i2] = obs[i1,i2]
        end
    end
    obs_criteraprime = prime(noprime(tensA.A*obs_criteria),noncommoninds(inds(tensA.A),physical_legs.A,ancilla_legs.A))
    obs_crit = obs_criteraprime*dag(tensA.A)
    obs_crit = obs_crit*cx.AF*cx.AB*cx.DA



    # MAGNETISATION

    sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];

    ssz = lattice("tensors")
    ssz_ = [lattice("tensors"),lattice("tensors"),lattice("tensors")]

    names = fieldnames(lattice)
    for i = 1:6
        setproperty!(ssz,names[i],ITensor(getproperty(physical_legs,names[i])', dag(getproperty(physical_legs,names[i]))))
    end

    szz = kron(id,id,sz) + kron(id,sz,id) + kron(sz,id,id)
    sz1 = kron(sx,sx,id) + kron(sy,sy,id) + kron(sz,sz,id);
    sz2 = kron(sx,id,sx) + kron(sy,id,sy) + kron(sz,id,sz);
    sz3 = kron(id,sx,sx) + kron(id,sy,sy) + kron(id,sz,sz);

    for i = 1:6
       
        obsz = ITensor(getproperty(physical_legs,names[i])',dag(getproperty(physical_legs,names[i])))

        all_magne = [ITensor(getproperty(physical_legs,names[i])',dag(getproperty(physical_legs,names[i]))),
        ITensor(getproperty(physical_legs,names[i])',dag(getproperty(physical_legs,names[i]))),
        ITensor(getproperty(physical_legs,names[i])',dag(getproperty(physical_legs,names[i])))]

        for i1 = 1:8, i2 = 1:8
            if abs(charges[i1] - charges[i2]) < 1e-6

                obsz[getproperty(physical_legs,names[i])'=>i1,getproperty(physical_legs,names[i])=>i2] = szz[i1,i2]
                all_magne[1][getproperty(physical_legs,names[i])'=>i1,getproperty(physical_legs,names[i])=>i2] = sz1[i1,i2]
                all_magne[2][getproperty(physical_legs,names[i])'=>i1,getproperty(physical_legs,names[i])=>i2] = sz2[i1,i2]
                all_magne[3][getproperty(physical_legs,names[i])'=>i1,getproperty(physical_legs,names[i])=>i2] = sz3[i1,i2]

            end

        end

        setproperty!(ssz, names[i], obsz)
        setproperty!(ssz_[1], names[i], all_magne[1])
        setproperty!(ssz_[2], names[i], all_magne[2])
        setproperty!(ssz_[3], names[i], all_magne[3])

    end

    mz = lattice("tensors")
    m_ = [lattice("tensors"),lattice("tensors"),lattice("tensors")]

    for i = 1:6

        tensor = getproperty(tensA, names[i])

        tmp_z_prime = prime(noprime(tensor*getproperty(ssz,names[i])),noncommoninds(inds(tensor),getproperty(physical_legs,names[i]),getproperty(ancilla_legs,names[i])))
        tmp_z = tmp_z_prime*dag(tensor)

        tmp_z_prime1 = prime(noprime(tensor*getproperty(ssz_[1],names[i])),noncommoninds(inds(tensor),getproperty(physical_legs,names[i]),getproperty(ancilla_legs,names[i])))
        tmp_z_prime2 = prime(noprime(tensor*getproperty(ssz_[2],names[i])),noncommoninds(inds(tensor),getproperty(physical_legs,names[i]),getproperty(ancilla_legs,names[i])))
        tmp_z_prime3 = prime(noprime(tensor*getproperty(ssz_[3],names[i])),noncommoninds(inds(tensor),getproperty(physical_legs,names[i]),getproperty(ancilla_legs,names[i])))
        tmp_z1 = tmp_z_prime1*dag(tensor)
        tmp_z2 = tmp_z_prime2*dag(tensor)
        tmp_z3 = tmp_z_prime3*dag(tensor)

        if i == 1 

            setproperty!(mz, names[i], tmp_z*cx.AF*cx.AB*cx.DA) 
            setproperty!(m_[1], names[i], tmp_z1*cx.AF*cx.AB*cx.DA) 
            setproperty!(m_[2], names[i], tmp_z2*cx.AF*cx.AB*cx.DA) 
            setproperty!(m_[3], names[i], tmp_z3*cx.AF*cx.AB*cx.DA) 

        elseif i == 2

            setproperty!(mz, names[i], tmp_z*dag(cx.AB)*dag(cx.BC)*dag(cx.EB)) 
            setproperty!(m_[1], names[i], tmp_z1*dag(cx.AB)*dag(cx.BC)*dag(cx.EB)) 
            setproperty!(m_[2], names[i], tmp_z2*dag(cx.AB)*dag(cx.BC)*dag(cx.EB)) 
            setproperty!(m_[3], names[i], tmp_z3*dag(cx.AB)*dag(cx.BC)*dag(cx.EB)) 

        elseif i == 3

            setproperty!(mz, names[i], tmp_z*cx.BC*cx.CD*cx.CF) 
            setproperty!(m_[1], names[i], tmp_z1*cx.BC*cx.CD*cx.CF) 
            setproperty!(m_[2], names[i], tmp_z2*cx.BC*cx.CD*cx.CF) 
            setproperty!(m_[3], names[i], tmp_z3*cx.BC*cx.CD*cx.CF) 

        elseif i == 4

            setproperty!(mz, names[i], tmp_z*dag(cx.DA)*dag(cx.CD)*dag(cx.ED)) 
            setproperty!(m_[1], names[i], tmp_z1*dag(cx.DA)*dag(cx.CD)*dag(cx.ED)) 
            setproperty!(m_[2], names[i], tmp_z2*dag(cx.DA)*dag(cx.CD)*dag(cx.ED)) 
            setproperty!(m_[3], names[i], tmp_z3*dag(cx.DA)*dag(cx.CD)*dag(cx.ED)) 

        elseif i == 5

            setproperty!(mz, names[i], tmp_z*cx.FE*cx.ED*cx.EB) 
            setproperty!(m_[1], names[i], tmp_z1*cx.FE*cx.ED*cx.EB) 
            setproperty!(m_[2], names[i], tmp_z2*cx.FE*cx.ED*cx.EB) 
            setproperty!(m_[3], names[i], tmp_z3*cx.FE*cx.ED*cx.EB) 

        elseif i == 6

            setproperty!(mz, names[i], tmp_z*dag(cx.AF)*dag(cx.CF)*dag(cx.FE)) 
            setproperty!(m_[1], names[i], tmp_z1*dag(cx.AF)*dag(cx.CF)*dag(cx.FE)) 
            setproperty!(m_[2], names[i], tmp_z2*dag(cx.AF)*dag(cx.CF)*dag(cx.FE)) 
            setproperty!(m_[3], names[i], tmp_z3*dag(cx.AF)*dag(cx.CF)*dag(cx.FE)) 

        end

    end
  
    # ENERGY

    obs = get_obs_tensors(hamil, tensA, physical_legs, ancilla_legs, cx)
    # @showtime hexx1 = get_obs_tensors(hex1, tensA, physical_legs, ancilla_legs, cx)
    # @showtime hexx2 = get_obs_tensors(hex2, tensA, physical_legs, ancilla_legs, cx)

   
    iza = commonind(tensa.A,tensa.F)
    iya = commonind(tensa.A,tensa.B)
    ixa = commonind(tensa.A,tensa.D)

    ize = commonind(tensa.E,tensa.B)
    ixe = commonind(tensa.E,tensa.F)
    iye = commonind(tensa.E,tensa.D)

    izc = commonind(tensa.C,tensa.D)
    ixc = commonind(tensa.C,tensa.B)
    iyc = commonind(tensa.C,tensa.F)

    
    bd_init = 1

    xa,ya,za,xc,yc,zc,xe,ye,ze = [Index(QN(0)=>bd_init) for _ in 1:9]

    T1 = lattice(randomITensor(xa,(iza),dag(xc)),
    randomITensor(xe,dag(iya),dag(xc)),
    randomITensor(xc,(izc),dag(xe)),
    randomITensor(xc,dag(iye),dag(xa)),
    randomITensor(xe,(ize),dag(xa)),
    randomITensor(xa,dag(iyc),dag(xe)))

    T2 = lattice(randomITensor((ya),(ixa),dag(yc)),
    randomITensor((yc),dag(ize),dag(ya)),
    randomITensor((yc),(ixc),dag(ye)),
    randomITensor((ya),dag(izc),dag(ye)),
    randomITensor((ye),(ixe),dag(ya)),
    randomITensor((ye),dag(iza),dag(yc)))

    T3 = lattice(randomITensor((za),(iya),dag(zc)),
    randomITensor((za),dag(ixc),dag(ze)),
    randomITensor((zc),(iyc),dag(ze)),
    randomITensor((ze),dag(ixa),dag(zc)),
    randomITensor((ze),(iye),dag(za)),
    randomITensor((zc),dag(ixe),dag(za)))

    C1cf = randomITensor(zc,dag(xe));    # C1cf[zc=>1,xe=>1] = 1
    C1ab = randomITensor(za,dag(xc));    # C1ab[za=>1, xc=>1] = 1
    C1ed = randomITensor(ze,dag(xa));    # C1ed[ze=>1, xa=>1] = 1
    C2af = randomITensor((xa),dag(yc));  # C2af[xa=>1, yc=>1] = 1
    C2cd = randomITensor((xc),dag(ye));  # C2cd[xc=>1, ye=>1] = 1
    C2eb = randomITensor((xe),dag(ya));  # C2eb[xe=>1,ya=>1] = 1
    C3da = randomITensor(ya,dag(zc));    #C3da[ya=>1, zc=>1] = 1
    C3fe = randomITensor(ye,dag(za));    #C3fe[ye=>1, za=>1] = 1
    C3bc = randomITensor(yc,dag(ze));    #C3bc[yc=>1, ze=>1] = 1


    renormalise_lattice_to_one!(T1,T2,T3)
    renormalise_tensor_to_one!(C1cf,C1ab,C1ed,C2af,C2cd,C2eb,C3da,C3fe,C3bc)

    chi = min(bond_dimension*bond_dimension + 1, 40)

    ener = 0
    eps = 1e-6
    err = 1
    maz = 2


    while err > eps

        @showtime T1,T2,T3,C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da = evolution_ctm(T1,T2,T3,C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da,tensa,chi)

        enertmp = copy(ener)
        maz_tmp = copy(maz)

        C11 = C1cf*T1.E*T3.D

        C1 = C11*(tensa.A*tensa.B)
        C2 = C2eb*T1.F*T2.A*(tensa.C*tensa.D)
        C3 = C3da*T2.B*T3.C*(tensa.E*tensa.F)
        C23 = C2*C3
        Z = (C1*C23)[1]

        C1af = C11*(obs_crit*tensa.B)
        maz = (C1af*C23)[1]/Z[1]
    
        @show err = norm(maz - maz_tmp)

    end



    
    C1 = C1cf*T1.E*T3.D*(tensa.A*tensa.B)
    C2 = C2eb*T1.F*T2.A*(tensa.C*tensa.D)
    C3 = C3da*T2.B*T3.C*(tensa.E*tensa.F)
    Z = (C1*C2*C3)[]

    C1_e = C1cf*T1.E*T3.D*obs.AB
    C2_e = C2eb*T1.F*T2.A*obs.CD
    C3_e = C3da*T2.B*T3.C*obs.FE
    eab = real((C1_e*C2*C3)[]/Z)
    ecd = real((C1*C2_e*C3)[]/Z)
    efe = real((C1*C2*C3_e)[]/Z)
    ener = eab + ecd + efe

   
    

    # ener
    C1 = C1ab*T1.C*T3.F*(tensa.E*tensa.D)
    C2 = C2cd*T1.B*T2.E*(tensa.A*tensa.F)
    C3 = C3fe*T2.D*T3.A*(tensa.B*tensa.C)
    Z = (C1*C2*C3)[]

    C1_e = C1ab*T1.C*T3.F*obs.ED
    C2_e = C2cd*T1.B*T2.E*obs.AF
    C3_e = C3fe*T2.D*T3.A*obs.BC
    eed = real((C1_e*C2*C3)[]/Z)
    eaf = real((C1*C2_e*C3)[]/Z)
    ecb = real((C1*C2*C3_e)[]/Z)
    ener = ener + eed + eaf + ecb

    
    
    # ener 
    C1 = C1ed*T1.A*T3.B*(tensa.C*tensa.F)
    C2 = C2af*T1.D*T2.C*(tensa.E*tensa.B)
    C3 = C3bc*T2.F*T3.E*(tensa.A*tensa.D)
    Z = (C1*C2*C3)[]

    C1_e = C1ed*T1.A*T3.B*obs.CF
    C2_e = C2af*T1.D*T2.C*obs.EB
    C3_e = C3bc*T2.F*T3.E*obs.DA
    ecf = real((C1_e*C2*C3)[]/Z)
    eeb = real((C1*C2_e*C3)[]/Z)
    ead = real((C1*C2*C3_e)[]/Z)

   
    
    ener = ener + ecf + eeb + ead
    ener = ener/18
    
     
    C1 = C1ab*T1.C*T3.F*(tensa.E*tensa.D)
    C2 = C2cd*T1.B*T2.E*(tensa.A*tensa.F)
    C3 = C3fe*T2.D*T3.A*(tensa.C*tensa.B)
    Z = Array(C1*C2*C3)[1]

    
    

    C1mD = C1ab*T1.C*T3.F*(tensa.E*mz.D)
    C1mE = C1ab*T1.C*T3.F*(mz.E*tensa.D)
    C3mB = C3fe*T2.D*T3.A*(tensa.C*mz.B)
    C3mC = C3fe*T2.D*T3.A*(mz.C*tensa.B)
    C2mF = C2cd*T1.B*T2.E*(tensa.A*mz.F)
    C2mA = C2cd*T1.B*T2.E*(mz.A*tensa.F)

    magne = lattice("tensors")
    magne.A = (C1*C2mA*C3)[1]/Z
    magne.B = (C1*C2*C3mB)[1]/Z
    magne.C = (C1*C2*C3mC)[1]/Z
    magne.D = (C1mD*C2*C3)[1]/Z
    magne.E = (C1mE*C2*C3)[1]/Z
    magne.F = (C1*C2mF*C3)[1]/Z

   
    return ener,T1,T2,T3,C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, magne#, SS1, SS2, magne, all_magne_


end
