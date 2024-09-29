

function simpleupdate(Gamma, lambda, ix, iy, iz, physical_legs, ancilla_legs, h, D)

    # z : AF CD EB
    # AF
    # AA = Gamma.A*lambda.DA*lambda.AB 
    # FF = (Gamma.F*lambda.AF)*lambda.CF*lambda.FE
    AA = Gamma.A*lambda.DA*lambda.AB 
    FF = Gamma.F*lambda.CF*lambda.FE


    ua,sa,va = svd(AA,[ix.D, iy.B, ancilla_legs.A])
    uf,sf,vf = svd(FF,[iz.F, physical_legs.F])
    index_a = commonind(ua,sa)
    mm = (sa*va*lambda.AF*uf*sf)*h.AF
    mm = noprime(mm)
    u,s,v = svd(mm, [index_a, physical_legs.A], maxdim = D) 
    lambda.AF = s/norm(s)
    Gamma.A = u*ua*inv_lambda(lambda.DA)*inv_lambda(lambda.AB)
    Gamma.F = v*vf*inv_lambda(lambda.CF)*inv_lambda(lambda.FE)
    iz.A = commonind(Gamma.A, lambda.AF)
    iz.F = commonind(lambda.AF, Gamma.F)

    # CD
    CC = Gamma.C*lambda.CF*lambda.BC
    DD = Gamma.D*lambda.ED*lambda.DA
    uc,sc,vc = svd(CC, [ix.B, iy.F, ancilla_legs.C])
    ud,sd,vd = svd(DD, [iz.D, physical_legs.D])
    index_c = commonind(uc,sc)
    mm = (sc*vc*lambda.CD*ud*sd)*h.CD
    mm = noprime(mm)
    u,s,v = svd(mm, [index_c, physical_legs.C], maxdim = D) 
    lambda.CD = s/norm(s)
    Gamma.C = u*uc*inv_lambda(lambda.BC)*inv_lambda(lambda.CF)
    Gamma.D = v*vd*inv_lambda(lambda.ED)*inv_lambda(lambda.DA)
    iz.C = commonind(Gamma.C, lambda.CD)
    iz.D = commonind(lambda.CD, Gamma.D)

    # EB
    EE = Gamma.E*lambda.FE*lambda.ED
    BB = Gamma.B*lambda.AB*lambda.BC
    ue,se,ve = svd(EE, [ix.F, iy.D, ancilla_legs.E])
    ub,sb,vb = svd(BB, [iz.B, physical_legs.B])
    index_e = commonind(ue,se)
    mm = (se*ve*lambda.EB*ub*sb)*h.EB
    mm = noprime(mm)
    u,s,v = svd(mm, [index_e, physical_legs.E], maxdim = D) 
    lambda.EB = s/norm(s)
    Gamma.E = u*ue*inv_lambda(lambda.FE)*inv_lambda(lambda.ED)
    Gamma.B = v*vb*inv_lambda(lambda.AB)*inv_lambda(lambda.BC)
    iz.E = commonind(Gamma.E,lambda.EB)
    iz.B = commonind(lambda.EB, Gamma.B)

    ############ y ############
    # A B
    AA = Gamma.A*lambda.AF*lambda.DA
    BB = Gamma.B*lambda.BC*lambda.EB
    ua,sa,va = svd(AA, [ix.D, iz.F, ancilla_legs.A])
    ub,sb,vb = svd(BB, [iy.B, physical_legs.B])
    index_a = commonind(ua,sa)
    mm = (sa*va*lambda.AB*ub*sb)*h.AB;
    mm = noprime(mm)
    u,s,v = svd(mm, [index_a, physical_legs.A], maxdim = D);
    lambda.AB = s/norm(s)
    Gamma.A = u*ua*inv_lambda(lambda.DA)*inv_lambda(lambda.AF)
    Gamma.B = v*vb*inv_lambda(lambda.EB)*inv_lambda(lambda.BC)
    iy.A = commonind(Gamma.A, lambda.AB)
    iy.B = commonind(lambda.AB, Gamma.B)

    # C F
    CC = Gamma.C*lambda.BC*lambda.CD
    FF = ((Gamma.F)*lambda.AF)*lambda.FE
    uc,sc,vc = svd(CC, [ix.B, iz.D, ancilla_legs.C])
    uf,sf,vf = svd(FF, [iy.F, physical_legs.F])
    index_c = commonind(uc,sc)
    mm = (sc*vc*lambda.CF*uf*sf)*h.CF;
    mm = noprime(mm)
    u,s,v = svd(mm, [index_c, physical_legs.C], maxdim = D);
    lambda.CF = s/norm(s)
    Gamma.C = u*uc*inv_lambda(lambda.BC)*inv_lambda(lambda.CD)
    Gamma.F = v*vf*inv_lambda(lambda.AF)*inv_lambda(lambda.FE)
    iy.C = commonind(Gamma.C, lambda.CF)
    iy.F = commonind(lambda.CF, Gamma.F)

    # E D
    EE = Gamma.E*lambda.FE*lambda.EB
    DD = Gamma.D*lambda.CD*lambda.DA
    ue,se,ve = svd(EE, [ix.F, iz.B, ancilla_legs.E])
    ud,sd,vd = svd(DD, [iy.D, physical_legs.D])
    index_e = commonind(ue,se)
    mm = (se*ve*lambda.ED*ud*sd)*h.ED;
    mm = noprime(mm)
    u,s,v = svd(mm, [index_e, physical_legs.E], maxdim = D);
    lambda.ED = s/norm(s)
    Gamma.E = u*ue*inv_lambda(lambda.FE)*inv_lambda(lambda.EB)
    Gamma.D = v*vd*inv_lambda(lambda.CD)*inv_lambda(lambda.DA)
    iy.E = commonind(Gamma.E, lambda.ED)
    iy.D = commonind(lambda.ED, Gamma.D)

    ############ x ############
    # x : DA BC FE
    # D A
    DD = Gamma.D*lambda.ED*lambda.CD
    AA = Gamma.A*lambda.AB*lambda.AF
    ud,sd,vd = svd(DD, [iy.E, iz.C, ancilla_legs.D])
    ua,sa,va = svd(AA, [ix.A, physical_legs.A])
    index_d = commonind(ud,sd)
    mm = (sd*vd*lambda.DA*ua*sa)*h.DA
    mm = noprime(mm)
    u,s,v = svd(mm, [index_d, physical_legs.D], maxdim = D)
    lambda.DA = s/norm(s)
    Gamma.D = u*ud*inv_lambda(lambda.ED)*inv_lambda(lambda.CD)
    Gamma.A = v*va*inv_lambda(lambda.AB)*inv_lambda(lambda.AF)
    ix.D = commonind(Gamma.D, lambda.DA)
    ix.A = commonind(lambda.DA, Gamma.A)

    # B C
    BB = Gamma.B*lambda.AB*lambda.EB 
    CC = Gamma.C*lambda.CD*lambda.CF
    ub,sb,vb = svd(BB, [iy.A, iz.E, ancilla_legs.B])
    uc,sc,vc = svd(CC, [ix.C, physical_legs.C])
    index_b = commonind(ub,sb)
    mm = (sb*vb*lambda.BC*uc*sc)*h.BC
    mm = noprime(mm)
    u,s,v = svd(mm, [index_b, physical_legs.B], maxdim = D)
    lambda.BC = s/norm(s)
    Gamma.B = u*ub*inv_lambda(lambda.AB)*inv_lambda(lambda.EB)
    Gamma.C = v*vc*inv_lambda(lambda.CD)*inv_lambda(lambda.CF)
    ix.B = commonind(Gamma.B, lambda.BC)
    ix.C = commonind(lambda.BC, Gamma.C)


    # F E
    FF = Gamma.F*lambda.CF*lambda.AF
    EE = Gamma.E*lambda.EB*lambda.ED
    uf,sf,vf = svd(FF, [iy.C, iz.A, ancilla_legs.F])
    ue,se,ve = svd(EE, [ix.E, physical_legs.E])
    index_f = commonind(uf,sf)
    mm = (sf*vf*lambda.FE*ue*se)*h.FE
    mm = noprime(mm)
    u,s,v = svd(mm, [index_f, physical_legs.F], maxdim = D)
    lambda.FE = s/norm(s)
    Gamma.F = u*uf*inv_lambda(lambda.CF)*inv_lambda(lambda.AF)
    Gamma.E = v*ve*inv_lambda(lambda.EB)*inv_lambda(lambda.ED)
    ix.F = commonind(Gamma.F, lambda.FE)
    ix.E = commonind(lambda.FE, Gamma.E)

    return Gamma, lambda, ix, iy, iz

end
