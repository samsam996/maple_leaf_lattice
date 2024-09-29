
include("inverse_matrixM.jl")
include("new_qp.jl")

function fullupdate(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate,D)

    tensA, tensa = fullupdateAB(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate)
    tensA, tensa = fullupdateCF(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate)
    tensA, tensa = fullupdateED(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate)

    tensA, tensa = fullupdateAF(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate)
    tensA, tensa = fullupdateCD(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate)
    tensA, tensa = fullupdateEB(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate)

    tensA, tensa = fullupdateDA(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate)
    tensA, tensa = fullupdateBC(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate)
    tensA, tensa = fullupdateFE(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter_gate)

    return tensA, tensa

end


function fullupdateAB(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter)


    # I need 9times what is below
    # AB
    
    # create the environment 
    C1 = C1cf*T1.E*T3.D
    C2 = C2eb*T1.F*T2.A*(tensa.C*tensa.D)
    C3 = C3da*T2.B*T3.C*(tensa.E*tensa.F)

    environment = C1*C2*C3

    iab = commonind(tensA.A, tensA.B)
    ida = commonind(tensA.D, tensA.A)
    iaf = commonind(tensA.A, tensA.F)

    ua,sa,va = svd(tensA.A, (ida,iaf)) 
    ub,sb,vb = svd(tensA.B, (physical_legs.B, iab)) 
    
    p = sa*va;
    q = ub*sb

    uua = (dag(ua)*(prime(ua,noncommoninds(inds(ua),physical_legs.A))))
    uua = uua*cx.DA
    uua = uua*cx.AF
    vvb = (dag(vb)*(prime(vb,noncommoninds(inds(vb),physical_legs.B))))
    vvb = vvb*dag(cx.EB)
    vvb = vvb*dag(cx.BC)
    environment = environment*uua
    environment = environment*vvb

    p, q = new_qp(environment, p, q, trotter.AB, physical_legs.A, physical_legs.B)

    tensA.A = ua*p; 
    tensA.B = q*vb; 
       
    Aprime = prime(tensA.A, noncommoninds(inds(tensA.A), [physical_legs.A]))
    Bprime = prime(tensA.B, noncommoninds(inds(tensA.B), [physical_legs.B]))
            
    tensa.A = (Aprime*dag(tensA.A))*cx.AF*cx.AB*cx.DA
    tensa.B = (Bprime*dag(tensA.B))*dag(cx.AB)*dag(cx.BC)*dag(cx.EB)
    


    return tensA, tensa


end

function fullupdateCF(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter)


    # I need 9times what is below
    # AB
    
    # create the environment 
    C1 = C1ed*T1.A*T3.B
    C2 = C2af*T1.D*T2.C*(tensa.E*tensa.B)
    C3 = C3bc*T2.F*T3.E*(tensa.A*tensa.D)

    environment = C1*C2*C3

    icf = commonind(tensA.C, tensA.F)
    ibc = commonind(tensA.B, tensA.C)
    icd = commonind(tensA.C, tensA.D)

    uc,sc,vc = svd(tensA.C, (ibc,icd)) 
    uf,sf,vf = svd(tensA.F, (physical_legs.F, icf)) 
    
    p = sc*vc;
    q = uf*sf

    uuc = (dag(uc)*(prime(uc,noncommoninds(inds(uc),physical_legs.C))))
    uuc = uuc*cx.BC
    uuc = uuc*cx.CD
    vvf = (dag(vf)*(prime(vf,noncommoninds(inds(vf),physical_legs.F))))
    vvf = vvf*dag(cx.AF)
    vvf = vvf*dag(cx.FE)
    environment = environment*uuc
    environment = environment*vvf

    p, q = new_qp(environment, p, q, trotter.CF, physical_legs.C, physical_legs.F)

    tensA.C = uc*p; 
    tensA.F = q*vf; 
       
    Cprime = prime(tensA.C, noncommoninds(inds(tensA.C), [physical_legs.C]))
    Fprime = prime(tensA.F, noncommoninds(inds(tensA.F), [physical_legs.F]))
            
    tensa.C = (Cprime*dag(tensA.C))*cx.CD*cx.CF*cx.BC
    tensa.F = (Fprime*dag(tensA.F))*dag(cx.CF)*dag(cx.FE)*dag(cx.AF)
    


    return tensA, tensa


end

function fullupdateED(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter)


    # I need 9times what is below
    # ED
    
    # create the environment 
    C1 = C1ab*T1.C*T3.F
    C2 = C2cd*T1.B*T2.E*(tensa.A*tensa.F)
    C3 = C3fe*T2.D*T3.A*(tensa.C*tensa.B)    

    environment = C1*C2*C3

    ied = commonind(tensA.E, tensA.D)
    ife = commonind(tensA.F, tensA.E)
    ieb = commonind(tensA.E, tensA.B)

    ue,se,ve = svd(tensA.E, (ife,ieb)) 
    ud,sd,vd = svd(tensA.D, (physical_legs.D, ied)) 
    
    p = se*ve;
    q = ud*sd

    uue = (dag(ue)*(prime(ue,noncommoninds(inds(ue),physical_legs.E))))
    uue = uue*cx.FE
    uue = uue*cx.EB
    vvd = (dag(vd)*(prime(vd,noncommoninds(inds(vd),physical_legs.D))))
    vvd = vvd*dag(cx.CD)
    vvd = vvd*dag(cx.DA)
    environment = environment*uue
    environment = environment*vvd

    p, q = new_qp(environment, p, q, trotter.ED, physical_legs.E, physical_legs.D)

    tensA.E = ue*p; 
    tensA.D = q*vd; 
       
    Eprime = prime(tensA.E, noncommoninds(inds(tensA.E), [physical_legs.E]))
    Dprime = prime(tensA.D, noncommoninds(inds(tensA.D), [physical_legs.D]))
            
    tensa.E = (Eprime*dag(tensA.E))*cx.EB*cx.ED*cx.FE
    tensa.D = (Dprime*dag(tensA.D))*dag(cx.ED)*dag(cx.DA)*dag(cx.CD)
    


    return tensA, tensa


end

function fullupdateAF(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter)


    # I need 9times what is below
    # AB
    
    # create the environment 
    C1 = C1ab*T1.C*T3.F*(tensa.E*tensa.D)
    C2 = C2cd*T1.B*T2.E
    C3 = C3fe*T2.D*T3.A*(tensa.C*tensa.B)    

    environment = C1*C2*C3

    iaf = commonind(tensA.A, tensA.F)
    iab = commonind(tensA.A, tensA.B)
    ida = commonind(tensA.D, tensA.A)

    ua,sa,va = svd(tensA.A, (ida,iab)) 
    uf,sf,vf = svd(tensA.F, (physical_legs.F, iaf)) 
    
    p = sa*va;
    q = uf*sf

    uua = (dag(ua)*(prime(ua,noncommoninds(inds(ua),physical_legs.A))))
    uua = uua*cx.DA
    uua = uua*cx.AB
    vvf = (dag(vf)*(prime(vf,noncommoninds(inds(vf),physical_legs.F))))
    vvf = vvf*dag(cx.CF)
    vvf = vvf*dag(cx.FE)
    environment = environment*uua
    environment = environment*vvf

    p, q = new_qp(environment, p, q, trotter.AF, physical_legs.A, physical_legs.F)

    tensA.A = ua*p; 
    tensA.F = q*vf; 
       
    Aprime = prime(tensA.A, noncommoninds(inds(tensA.A), [physical_legs.A]))
    Fprime = prime(tensA.F, noncommoninds(inds(tensA.F), [physical_legs.F]))
            
    tensa.A = (Aprime*dag(tensA.A))*cx.AF*cx.AB*cx.DA
    tensa.F = (Fprime*dag(tensA.F))*dag(cx.AF)*dag(cx.CF)*dag(cx.FE)
    


    return tensA, tensa


end

function fullupdateCD(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter)


    # I need 9times what is below
    # CD 
    C1 = C1cf*T1.E*T3.D*(tensa.A*tensa.B)
    C2 = C2eb*T1.F*T2.A
    C3 = C3da*T2.B*T3.C*(tensa.E*tensa.F)
    
    environment = C1*C2*C3
    ## me suis arreté ici
    icd = commonind(tensA.C, tensA.D)
    icf = commonind(tensA.C, tensA.F)
    ibc = commonind(tensA.B, tensA.C)

    uc,sc,vc = svd(tensA.C, (ibc,icf)) 
    ud,sd,vd = svd(tensA.D, (physical_legs.D, icd)) 
    
    p = sc*vc;
    q = ud*sd

    uuc = (dag(uc)*(prime(uc,noncommoninds(inds(uc),physical_legs.C))))
    uuc = uuc*cx.BC
    uuc = uuc*cx.CF
    vvd = (dag(vd)*(prime(vd,noncommoninds(inds(vd),physical_legs.D))))
    vvd = vvd*dag(cx.ED)
    vvd = vvd*dag(cx.DA)
    environment = environment*uuc
    environment = environment*vvd

    p, q = new_qp(environment, p, q, trotter.CD, physical_legs.C, physical_legs.D)

    tensA.C = uc*p; 
    tensA.D = q*vd; 
       
    Cprime = prime(tensA.C, noncommoninds(inds(tensA.C), [physical_legs.C]))
    Dprime = prime(tensA.D, noncommoninds(inds(tensA.D), [physical_legs.D]))
            
    tensa.C = (Cprime*dag(tensA.C))*cx.CD*cx.CF*cx.BC
    tensa.D = (Dprime*dag(tensA.D))*dag(cx.CD)*dag(cx.ED)*dag(cx.DA)
    


    return tensA, tensa


end

function fullupdateEB(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter)


    # I need 9times what is below
    # CD 
    C1 = C1ed*T1.A*T3.B*(tensa.C*tensa.F)
    C2 = C2af*T1.D*T2.C
    C3 = C3bc*T2.F*T3.E*(tensa.A*tensa.D)

    environment = C1*C2*C3

    ieb = commonind(tensA.E, tensA.B)
    ied = commonind(tensA.E, tensA.D)
    ife = commonind(tensA.F, tensA.E)

    ue,se,ve = svd(tensA.E, (ife,ied)) 
    ub,sb,vb = svd(tensA.B, (physical_legs.B, ieb)) 
    
    p = se*ve;
    q = ub*sb

    uue = (dag(ue)*(prime(ue,noncommoninds(inds(ue),physical_legs.E))))
    uue = uue*cx.FE
    uue = uue*cx.ED
    vvb = (dag(vb)*(prime(vb,noncommoninds(inds(vb),physical_legs.B))))
    vvb = vvb*dag(cx.AB)
    vvb = vvb*dag(cx.BC)
    environment = environment*uue
    environment = environment*vvb

    p, q = new_qp(environment, p, q, trotter.EB, physical_legs.E, physical_legs.B)

    tensA.E = ue*p; 
    tensA.B = q*vb; 
       
    Eprime = prime(tensA.E, noncommoninds(inds(tensA.E), [physical_legs.E]))
    Bprime = prime(tensA.B, noncommoninds(inds(tensA.B), [physical_legs.B]))
            
    tensa.E = (Eprime*dag(tensA.E))*cx.EB*cx.ED*cx.FE
    tensa.B = (Bprime*dag(tensA.B))*dag(cx.EB)*dag(cx.AB)*dag(cx.BC)
    


    return tensA, tensa


end

function fullupdateDA(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter)

    # DA
    
    C1 = C1ed*T1.A*T3.B*(tensa.C*tensa.F)
    C2 = C2af*T1.D*T2.C*(tensa.E*tensa.B)
    C3 = C3bc*T2.F*T3.E
    
    environment = C1*C2*C3

    ida = commonind(tensA.D, tensA.A)
    icd = commonind(tensA.C, tensA.D)
    ied = commonind(tensA.E, tensA.D)

    ud,sd,vd = svd(tensA.D, (icd,  ied)) 
    ua,sa,va = svd(tensA.A, (physical_legs.A, ida)) 
    
    p = sd*vd;
    q = ua*sa

    uud = (dag(ud)*(prime(ud,noncommoninds(inds(ud),physical_legs.D))))
    uud = uud*dag(cx.CD)
    uud = uud*dag(cx.ED)
    vva = (dag(va)*(prime(va,noncommoninds(inds(va),physical_legs.A))))
    vva = vva*(cx.AB)
    vva = vva*(cx.AF)
    environment = environment*uud
    environment = environment*vva

    p, q = new_qp2(environment, p, q, trotter.DA, physical_legs.D, physical_legs.A)

    tensA.D = ud*p; 
    tensA.A = q*va; 
       
    Dprime = prime(tensA.D, noncommoninds(inds(tensA.D), [physical_legs.D]))
    Aprime = prime(tensA.A, noncommoninds(inds(tensA.A), [physical_legs.A]))
            
    tensa.D = (Dprime*dag(tensA.D))*dag(cx.CD)*dag(cx.ED)*dag(cx.DA)
    tensa.A = (Aprime*dag(tensA.A))*(cx.DA)*(cx.AF)*(cx.AB)
    


    return tensA, tensa


end

function fullupdateBC(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter)

    # BC
    C1 = C1ab*T1.C*T3.F*(tensa.E*tensa.D)
    C2 = C2cd*T1.B*T2.E*(tensa.A*tensa.F)
    C3 = C3fe*T2.D*T3.A

    environment = C1*C2*C3

    ibc = commonind(tensA.B, tensA.C)
    ieb = commonind(tensA.E, tensA.B)
    iab = commonind(tensA.A, tensA.B)

    ub,sb,vb = svd(tensA.B, (iab,  ieb)) 
    uc,sc,vc = svd(tensA.C, (physical_legs.C, ibc)) 
    
    p = sb*vb;
    q = uc*sc

    uub = (dag(ub)*(prime(ub,noncommoninds(inds(ub),physical_legs.B))))
    uub = uub*dag(cx.AB)
    uub = uub*dag(cx.EB)
    vvc = (dag(vc)*(prime(vc,noncommoninds(inds(vc),physical_legs.C))))
    vvc = vvc*(cx.CF)
    vvc = vvc*(cx.CD)
    environment = environment*uub
    environment = environment*vvc

    p, q = new_qp2(environment, p, q, trotter.BC, physical_legs.B, physical_legs.C)

    tensA.B = ub*p; 
    tensA.C = q*vc; 
       
    Bprime = prime(tensA.B, noncommoninds(inds(tensA.B), [physical_legs.B]))
    Cprime = prime(tensA.C, noncommoninds(inds(tensA.C), [physical_legs.C]))
            
    tensa.B = (Bprime*dag(tensA.B))*dag(cx.AB)*dag(cx.EB)*dag(cx.BC)
    tensa.C = (Cprime*dag(tensA.C))*(cx.BC)*(cx.CD)*(cx.CF)
    


    return tensA, tensa


end

function fullupdateFE(tensa,tensA, cx, physical_legs, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da, T1, T2, T3, trotter)

    # FE
    C1 = C1cf*T1.E*T3.D*(tensa.A*tensa.B)
    C2 = C2eb*T1.F*T2.A*(tensa.C*tensa.D)
    C3 = C3da*T2.B*T3.C

    environment = C1*C2*C3

    ife = commonind(tensA.F, tensA.E)
    iaf = commonind(tensA.A, tensA.F)
    icf = commonind(tensA.C, tensA.F)

    uf,sf,vf = svd(tensA.F, (icf,  iaf)) 
    ue,se,ve = svd(tensA.E, (physical_legs.E, ife)) 
    
    p = sf*vf;
    q = ue*se

    uuf = (dag(uf)*(prime(uf,noncommoninds(inds(uf),physical_legs.F))))
    uuf = uuf*dag(cx.CF)
    uuf = uuf*dag(cx.AF)
    vve = (dag(ve)*(prime(ve,noncommoninds(inds(ve),physical_legs.E))))
    vve = vve*(cx.ED)
    vve = vve*(cx.EB)
    environment = environment*uuf
    environment = environment*vve

    p, q = new_qp2(environment, p, q, trotter.FE, physical_legs.F, physical_legs.E)

    tensA.F = uf*p; 
    tensA.E = q*ve; 
       
    Fprime = prime(tensA.F, noncommoninds(inds(tensA.F), [physical_legs.F]))
    Eprime = prime(tensA.E, noncommoninds(inds(tensA.E), [physical_legs.E]))
            
    tensa.F = (Fprime*dag(tensA.F))*dag(cx.CF)*dag(cx.AF)*dag(cx.FE)
    tensa.E = (Eprime*dag(tensA.E))*(cx.FE)*(cx.EB)*(cx.ED)
    


    return tensA, tensa


end

# function full_update!(tens_a::lattice,tens_A::lattice,cxd::lattice,
#     cyd::lattice,gt::Matrix{Symbol},physical_legs::lattice_ind,
#     gx::lattice,gy::lattice,D::Int64,C::Vector{lattice},T::Vector{lattice})

#     N = size(gt)[1]
#     f(x) = mod(x-1,N) + 1

#     precision_full_update = 1e-15

#     # IN THE X DIRECTION

#     for i = 1:N
#         for j = 1:1

#             ggx = getfield(gx, gt[f(i),f(j)])

#             C1a = getfield(C[1],gt[f(i-1),f(j-1)])
#             T1b = getfield(T[1],gt[f(i),f(j-1)])
#             T1a = getfield(T[1],gt[f(i+1),f(j-1)])
#             C2b = getfield(C[2],gt[f(i+2),f(j-1)])

#             T4c = getfield(T[4],gt[f(i-1),f(j)])
#             A = getfield(tens_A,gt[f(i),f(j)])
#             B = getfield(tens_A,gt[f(i+1),f(j)])
#             T2d = getfield(T[2],gt[f(i+2),f(j)])

#             C4a = getfield(C[4],gt[f(i-1),f(j+1)])
#             T3b = getfield(T[3],gt[f(i),f(j+1)])
#             T3a = getfield(T[3],gt[f(i+1),f(j+1)])
#             C3b = getfield(C[3],gt[f(i+2),f(j+1)])

#             Left = C1a*T4c*C4a
#             Right = C2b*T2d*C3b
        
#             ind_xa = commonind(A,getfield(cxd,gt[f(i),f(j)]))
#             ind_xb = commonind(B,getfield(cxd,gt[f(i+1),f(j)]))

#             physical_legs_ia = getfield(physical_legs,gt[f(i),f(j)])
#             # ancilla_legs_ia = getfield(ancilla_legs,gt[f(i),f(j)])
#             physical_legs_ib = getfield(physical_legs,gt[f(i+1),f(j)])
#             # ancilla_legs_ib = getfield(ancilla_legs,gt[f(i+1),f(j)])

#             indsA = inds(A)
#             indsB = inds(B)
#             all_but_ph_a = noncommoninds(indsA,[ind_xa, physical_legs_ia])
#             all_but_ph_b = noncommoninds(indsB,[ind_xb, physical_legs_ib])

#             u,s,v = svd(A,(all_but_ph_a))
#             u_index_p = commonind(u,s)
#             X = deepcopy(u); p = s*v
#             X_prime = prime(X)
            
#             u,s,v = svd(B,[ind_xa,physical_legs_ib])
#             q = u*s; Y = deepcopy(v); 
#             Y_prime = prime(Y)


#             p0 = deepcopy(p)
#             q0 = deepcopy(q)
#             q0_prime = prime(q0, noncommoninds(inds(q0), [physical_legs_ib]))
#             p0_prime = prime(p0, noncommoninds(inds(p0), [physical_legs_ia]))

#             #=      
#                         | dag(cy22)
#             - dag(cx22) - A  - cx11 -
#                         | cy11
#             =#

#             c1 = getfield(cxd,gt[f(i-1),f(j)])
#             c2 = getfield(cyd,gt[f(i),f(j)])
#             c3 = getfield(cyd,gt[f(i),f(j-1)])
#             XX2 = X_prime*dag(X)
#             XX = XX2*dag(c1)*c2*dag(c3)

#             c1 = getfield(cyd,gt[f(i+1),f(j)])
#             c2 = getfield(cxd,gt[f(i+1),f(j)])
#             c3 = getfield(cyd,gt[f(i+1),f(j-1)])
#             YY2 = Y_prime*dag(Y);
#             YY = YY2*(c1)*c2*dag(c3)

#             LL = (Left*T1b*XX*T3b)
#             RR = (Right*T1a*YY*T3a)
#             environment = LL*RR

#             ida = ITensor(physical_legs_ia, dag(physical_legs_ia)')
#             idb = ITensor(physical_legs_ib, dag(physical_legs_ib)')

#             for ii = 1:4
#                 ida[physical_legs_ia=>ii, dag(physical_legs_ia)' =>ii] = 1;
#                 idb[physical_legs_ib=>ii, dag(physical_legs_ib)' =>ii] = 1;
#             end

#             distance = 1
#             err_fullupdate = 1;

#             q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ib]))
#             p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))

#             bb = ((environment*p0_prime*q0_prime))*ggx;

#             itx = 0; 

#             while distance > precision_full_update && itx < 50


#                 itx = itx + 1

#                 q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ib]))
#                 p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))

#                 M = ((environment*q_prime)*dag(q))*ida
                
#                 b = noprime(bb,physical_legs_ib')
#                 b = noprime(b*dag(q))

#                 indd1 = commoninds(b,M)
#                 c = combiner(indd1, dir = -dir(physical_legs_ia))
#                 indd2 = noncommoninds(inds(M), indd1)
#                 cprime = combiner([indd2[1], indd2[3], indd2[2]])

#                 M_matrix = M*c*cprime
 
#                 inverse_M_matrix = inverse_matrixM(M_matrix)
#                 inverse_M = inverse_M_matrix*(c)*(cprime)

#                 p = inverse_M*b
#                 p = noprime(p)

#                 q_prime = prime(q, noncommoninds(inds(q), [physical_legs_ib]))
#                 p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))

#                 M = (environment*p_prime)*dag(p)*idb
#                 b = noprime(bb,physical_legs_ia')
#                 b = noprime(b*dag(p))

#                 indd1 = commoninds(b,M)
#                 c = combiner(indd1, dir = dir(physical_legs_ib))
#                 indd2 = noncommoninds(inds(M), indd1)
#                 cprime = combiner([indd2[1],indd2[3],indd2[2]])

#                 M_matrix = M*c*cprime

#                 inverse_M_matrix = inverse_matrixM(M_matrix)
#                 inverse_M = inverse_M_matrix*(c)*(cprime)

#                 q = inverse_M*b
#                 q = noprime(q)
            
#                 q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ib]))
#                 p_prime = prime(p, noncommoninds(inds(p),[physical_legs_ia]))

#                 middle_term = p0_prime*q0_prime*ggx
#                 middle_term2 = noprime(middle_term,physical_legs_ia',physical_legs_ib')*ggx

#                 term1 = noprime(environment*middle_term2)*dag(p0)*dag(q0)
#                 term2 = environment*(q_prime*dag(q))*(p_prime*dag(p))
#                 term3 = noprime(environment*middle_term)*dag(p)*dag(q)

#                 Z = norm(array(((environment*q0_prime)*dag(q0))*p0_prime*dag(p0)))

#                 previous_err_fulluptade = deepcopy(err_fullupdate)
#                 err_fullupdate = norm(array(term1) + array(term2) - 2*array(term3))/Z
#                 @show distance = norm(err_fullupdate - previous_err_fulluptade)

#             end

#             p = p/norm(((p)))
#             q = q/norm(((q)))
        
#             A = X*p; 
#             B = q*Y; 
       
#             setproperty!(tens_A, gt[f(i),f(j)], A)
#             setproperty!(tens_A, gt[f(i+1),f(j)], B)
#             Aprime = prime(A, noncommoninds(inds(A), [physical_legs_ia]))
#             Bprime = prime(B, noncommoninds(inds(B), [physical_legs_ib]))
            
#             aaa = (Aprime*dag(A))*getfield(cxd,gt[f(i),f(j)])*dag(getfield(cxd,gt[f(i-1),f(j)]))*
#             getfield(cyd,gt[f(i),f(j)])*dag(getfield(cyd,gt[f(i),f(j-1)]))
#             bbb = Bprime*dag(B)*getfield(cxd,gt[f(i+1),f(j)])*dag(getfield(cxd,gt[f(i),f(j)]))*
#             getfield(cyd,gt[f(i+1),f(j)])*dag(getfield(cyd,gt[f(i+1),f(j-1)]))

#             setproperty!(tens_a,gt[f(i),f(j)], aaa)
#             setproperty!(tens_a,gt[f(i+1),f(j)], bbb)

#         end
#     end

#     # IN THE Y DIRECTION 

#     for i = 1:1
#         for j = 1:N

#             ggy = getfield(gy, gt[f(i),f(j)])

#             C1d = getfield(C[1],gt[f(i-1),f(j-1)])
#             T1c = getfield(T[1],gt[f(i),f(j-1)])
#             C2d = getfield(C[2],gt[f(i+1),f(j-1)])

#             T4b = getfield(T[4],gt[f(i-1),f(j)])
#             A = getfield(tens_A,gt[f(i),f(j)])
#             T2b = getfield(T[2],gt[f(i+1),f(j)])

#             T4d = getfield(T[4],gt[f(i-1),f(j+1)])
#             CC = getfield(tens_A,gt[f(i),f(j+1)])
#             T2d = getfield(T[2],gt[f(i+1),f(j+1)])

#             C4b = getfield(C[4],gt[f(i-1),f(j+2)])
#             T3a = getfield(T[3],gt[f(i),f(j+2)])
#             C3b = getfield(C[3],gt[f(i+1),f(j+2)])

#             Up = C1d*T1c*C2d
#             Down = C4b*T3a*C3b
        
#             ind_ya = commonind(A,getfield(cyd,gt[f(i),f(j)]))
#             ind_yc = commonind(CC,getfield(cyd,gt[f(i),f(j+1)]))

#             physical_legs_ia = getfield(physical_legs,gt[f(i),f(j)])
#             physical_legs_ic = getfield(physical_legs,gt[f(i),f(j+1)])

#             indsA = inds(A)
#             indsC = inds(CC)
#             all_but_ph_a = noncommoninds(indsA,[ind_ya, physical_legs_ia])
#             all_but_ph_c = noncommoninds(indsC,[ind_yc,physical_legs_ic])

#             u,s,v = svd(A,(all_but_ph_a))
#             u_index_p = commonind(u,s)
#             X = deepcopy(u); p = s*v
#             X_prime = prime(X)
            
#             u,s,v = svd(CC,[ind_ya,physical_legs_ic])
#             q = u*s; Y = deepcopy(v); 
#             Y_prime = prime(Y)

#             p0 = deepcopy(p)
#             q0 = deepcopy(q)
#             p0_prime = prime(p0, noncommoninds(inds(p0), [physical_legs_ia]))
#             q0_prime = prime(q0, noncommoninds(inds(q0), [physical_legs_ic]))

#             #=      
#                         | dag(cy22)
#             - dag(cx22) - A  - cx11 -
#                         | cy11
#             =#

#             c1 = getfield(cxd,gt[f(i-1),f(j)])
#             c2 = getfield(cyd,gt[f(i),f(j-1)])
#             c3 = getfield(cxd,gt[f(i),f(j)])
#             XX2 = X_prime*dag(X)
#             XX = XX2*dag(c1)*dag(c2)*(c3)

#             c1 = getfield(cxd,gt[f(i-1),f(j+1)])
#             c2 = getfield(cyd,gt[f(i),f(j+1)])
#             c3 = getfield(cxd,gt[f(i),f(j+1)])
#             YY2 = Y_prime*dag(Y);
#             YY = YY2*dag(c1)*(c2)*(c3)

#             UU = (Up*T4b*XX*T2b)
#             DD = (Down*T4d*YY*T2d)
#             environment = DD*UU

#             ida = ITensor(physical_legs_ia, dag(physical_legs_ia)')
#             idc = ITensor(physical_legs_ic, dag(physical_legs_ic)')

#             for ii = 1:4
#                 ida[physical_legs_ia=>ii, dag(physical_legs_ia)' =>ii] = 1;
#                 idc[physical_legs_ic=>ii, dag(physical_legs_ic)' =>ii] = 1;
#             end

#             distance = 1
#             err_fullupdate = 1;

#             p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))
#             q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ic]))

#             bb = ((environment*p0_prime*q0_prime))*ggy;

#             ity = 0;

#             while distance > precision_full_update && ity < 50

#                 ity = ity + 1; 

#                 p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))
#                 q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ic]))

#                 M = ((environment*q_prime)*dag(q))*ida
                
#                 b = noprime(bb,physical_legs_ic')
#                 b = noprime(b*dag(q))

#                 indd1 = commoninds(b,M)
#                 comb = combiner(indd1, dir = dir(physical_legs_ia))
#                 indd2 = noncommoninds(inds(M), indd1)
#                 cprime = combiner([indd2[1],indd2[3],indd2[2]], dir = -dir(physical_legs_ia))

#                 M_matrix = M*comb*cprime

#                 inverse_M_matrix = inverse_matrixM(M_matrix)
#                 inverse_M = inverse_M_matrix*(comb)*(cprime)

#                 p = inverse_M*b
#                 p = noprime(p)

#                 q_prime = prime(q, noncommoninds(inds(q), [physical_legs_ic]))
#                 p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))

#                 M = (environment*p_prime)*dag(p)*idc
#                 b = noprime(bb,physical_legs_ia')
#                 b = noprime(b*dag(p))

#                 indd1 = commoninds(b,M)
#                 c = combiner(indd1, dir = dir(physical_legs_ic))
#                 indd2 = noncommoninds(inds(M), indd1)
#                 cprime = combiner([indd2[1], indd2[3], indd2[2]])

#                 M_matrix = M*c*cprime

#                 inverse_M_matrix = inverse_matrixM(M_matrix)
#                 inverse_M = inverse_M_matrix*(c)*(cprime)

#                 q = inverse_M*b
#                 q = noprime(q)
            
#                 q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ic]))
#                 p_prime = prime(p, noncommoninds(inds(p),[physical_legs_ia]))

#                 middle_term = p0_prime*q0_prime*ggy
#                 middle_term2 = noprime(middle_term,physical_legs_ia',physical_legs_ic')*ggy

#                 term1 = noprime(environment*middle_term2)*dag(p0)*dag(q0)
#                 term2 = environment*(q_prime*dag(q))*(p_prime*dag(p))
#                 term3 = noprime(environment*middle_term)*dag(p)*dag(q)

#                 Z = norm(array(((environment*q0_prime)*dag(q0))*p0_prime*dag(p0)))

#                 previous_err_fulluptade = deepcopy(err_fullupdate)
#                 err_fullupdate = norm(array(term1) + array(term2) - 2*array(term3))/norm(array(term1))
#                 @show distance = norm(err_fullupdate - previous_err_fulluptade)

#             end

#             p = p/maximum(real(abs.(p)))
#             q = q/maximum(real(abs.(q)))

#             # old_ind = commonind(p,q)
#             # uu,ss,vv = svd(p*q, (u_index_p, physical_legs_ia), maxdim = D)
#             # ss = ss/maximum(ss)
#             # sqrt_s = diag_sqrt(ss)
#             # u_index = commonind(uu,ss)
#             # v_index = commonind(ss,vv)
#             # p = (uu*sqrt_s)*delta(dag(v_index),old_ind); 
#             # q = delta(dag(old_ind), (u_index))*(sqrt_s*vv)


#             A = X*p; 
#             CC = q*Y; 

#             setproperty!(tens_A, gt[f(i),f(j)], A)
#             setproperty!(tens_A, gt[f(i),f(j+1)], CC)
#             Aprime = prime(A, noncommoninds(inds(A), [physical_legs_ia]))
#             Cprime = prime(CC, noncommoninds(inds(CC), [physical_legs_ic]))
            
#             aaa = (Aprime*dag(A))*getfield(cxd,gt[f(i),f(j)])*dag(getfield(cxd,gt[f(i-1),f(j)]))*
#             getfield(cyd,gt[f(i),f(j)])*dag(getfield(cyd,gt[f(i),f(j-1)]))
#             ccc = Cprime*dag(CC)*getfield(cxd,gt[f(i),f(j+1)])*dag(getfield(cxd,gt[f(i-1),f(j+1)]))*
#             getfield(cyd,gt[f(i),f(j+1)])*dag(getfield(cyd,gt[f(i),f(j)]))

#             setproperty!(tens_a,gt[f(i),f(j)], aaa)
#             setproperty!(tens_a,gt[f(i),f(j+1)], ccc)

#         end
#     end
    
#     nothing

# end


