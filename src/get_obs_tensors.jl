


function get_obs_tensors(hex1, tensA, physical_legs, ancilla_legs, cx)

    obs = double_lattice("tensors")

   
    index_da = commonind(tensA.D, tensA.A)
    index_af = commonind(tensA.A, tensA.F)
    index_ab = commonind(tensA.A, tensA.B)
    ua,sa,va = svd(tensA.A, [index_da, index_af, ancilla_legs.A])
    ub,sb,vb = svd(tensA.B, [index_ab, physical_legs.B])
    up = sa*va*hex1.AB*ub*sb
    noprime!(up)
    up = prime(up,noncommoninds(inds(up),[physical_legs.A,physical_legs.B]))
    updown = up*dag(sa*va*ub*sb)
    uap = prime(ua,noncommoninds(inds(ua),[ancilla_legs.A]))
    vbp = prime(vb,noncommoninds(inds(vb),[ancilla_legs.B]))
    obs.AB = (uap*dag(ua))*updown*(vbp*dag(vb))
    obs.AB = obs.AB*cx.AF*cx.DA*dag(cx.BC)*dag(cx.EB)


    
    

    index_bc = commonind(tensA.B, tensA.C)
    index_cf = commonind(tensA.C, tensA.F)
    index_cd = commonind(tensA.C, tensA.D)
    uc,sc,vc = svd(tensA.C, [index_bc, index_cf, ancilla_legs.C])
    ud,sd,vd = svd(tensA.D, [index_cd, physical_legs.D])
    up = sc*vc*hex1.CD*ud*sd
    noprime!(up)
    up = prime(up,noncommoninds(inds(up),[physical_legs.C,physical_legs.D]))
    updown = up*dag(sc*vc*ud*sd)
    ucp = prime(uc,noncommoninds(inds(uc),[ancilla_legs.C]))
    vdp = prime(vd,noncommoninds(inds(vd),[ancilla_legs.D]))
    obs.CD = (ucp*dag(uc))*updown*(vdp*dag(vd))
    obs.CD = obs.CD*cx.BC*cx.CF*dag(cx.DA)*dag(cx.ED)


    
    index_af = commonind(tensA.A, tensA.F)
    index_cf = commonind(tensA.C, tensA.F)
    index_fe = commonind(tensA.F, tensA.E)
    uf,sf,vf = svd(tensA.F, [index_af, index_cf, ancilla_legs.F])
    ue,se,ve = svd(tensA.E, [index_fe, physical_legs.E])
    up = sf*vf*hex1.FE*ue*se
    noprime!(up)
    up = prime(up,noncommoninds(inds(up),[physical_legs.F,physical_legs.E]))
    updown = up*dag(sf*vf*ue*se)
    ufp = prime(uf,noncommoninds(inds(uf),[ancilla_legs.F]))
    vep = prime(ve,noncommoninds(inds(ve),[ancilla_legs.E]))
    obs.FE = (ufp*dag(uf))*updown*(vep*dag(ve))
    obs.FE = obs.FE*dag(cx.AF)*dag(cx.CF)*cx.EB*cx.ED

    

    index_fe = commonind(tensA.F, tensA.E)
    index_eb = commonind(tensA.E, tensA.B)
    index_ed = commonind(tensA.E, tensA.D)
    ue,se,ve = svd(tensA.E, [index_fe, index_eb, ancilla_legs.E])
    ud,sd,vd = svd(tensA.D, [index_ed, physical_legs.D])
    up = se*ve*hex1.ED*ud*sd
    noprime!(up)
    up = prime(up,noncommoninds(inds(up),[physical_legs.E,physical_legs.D]))
    updown = up*dag(se*ve*ud*sd)
    uep = prime(ue,noncommoninds(inds(ue),[ancilla_legs.E]))
    vdp = prime(vd,noncommoninds(inds(vd),[ancilla_legs.D]))
    obs.ED = (uep*dag(ue))*updown*(vdp*dag(vd))
    obs.ED = obs.ED*cx.FE*cx.EB*dag(cx.DA)*dag(cx.CD)


    

    index_da = commonind(tensA.D, tensA.A)
    index_ab = commonind(tensA.A, tensA.B)
    index_af = commonind(tensA.A, tensA.F)
    ua,sa,va = svd(tensA.A, [index_da, index_ab, ancilla_legs.A])
    uf,sf,vf = svd(tensA.F, [index_af, physical_legs.F])
    up = sa*va*hex1.AF*uf*sf
    noprime!(up)
    up = prime(up,noncommoninds(inds(up),[physical_legs.A,physical_legs.F]))
    updown = up*dag(sa*va*uf*sf)
    uap = prime(ua,noncommoninds(inds(ua),[ancilla_legs.A]))
    vfp = prime(vf,noncommoninds(inds(vf),[ancilla_legs.F]))
    obs.AF = (uap*dag(ua))*updown*(vfp*dag(vf))
    obs.AF = obs.AF*cx.DA*cx.AB*dag(cx.CF)*dag(cx.FE)




    index_ab = commonind(tensA.A, tensA.B)
    index_eb = commonind(tensA.E, tensA.B)
    index_bc = commonind(tensA.B, tensA.C)
    ub,sb,vb = svd(tensA.B, [index_ab, index_eb, ancilla_legs.B])
    uc,sc,vc = svd(tensA.C, [index_bc, physical_legs.C])
    up = sb*vb*hex1.BC*uc*sc
    noprime!(up)
    up = prime(up,noncommoninds(inds(up),[physical_legs.B, physical_legs.C]))
    updown = up*dag(sb*vb*uc*sc)
    ubp = prime(ub,noncommoninds(inds(ub),[ancilla_legs.B]))
    vcp = prime(vc,noncommoninds(inds(vc),[ancilla_legs.C]))
    obs.BC = (ubp*dag(ub))*updown*(vcp*dag(vc))
    obs.BC = obs.BC*cx.CF*cx.CD*dag(cx.EB)*dag(cx.AB)



    
    index_bc = commonind(tensA.B, tensA.C)
    index_cd = commonind(tensA.C, tensA.D)
    index_cf = commonind(tensA.C, tensA.F)
    uc,sc,vc = svd(tensA.C, [index_bc, index_cd, ancilla_legs.C])
    uf,sf,vf = svd(tensA.F, [index_cf, physical_legs.F])
    up = sc*vc*hex1.CF*uf*sf
    noprime!(up)
    up = prime(up,noncommoninds(inds(up),[physical_legs.C,physical_legs.F]))
    updown = up*dag(sc*vc*uf*sf)
    ucp = prime(uc,noncommoninds(inds(uc),[ancilla_legs.C]))
    vfp = prime(vf,noncommoninds(inds(vf),[ancilla_legs.F]))
    obs.CF = (ucp*dag(uc))*updown*(vfp*dag(vf))
    obs.CF = obs.CF*cx.CD*cx.BC*dag(cx.AF)*dag(cx.FE)



   
    index_fe = commonind(tensA.F, tensA.E)
    index_ed = commonind(tensA.E, tensA.D)
    index_eb = commonind(tensA.E, tensA.B)
    ue,se,ve = svd(tensA.E, [index_fe, index_ed, ancilla_legs.E])
    ub,sb,vb = svd(tensA.B, [index_eb, physical_legs.B])
    up = se*ve*hex1.EB*ub*sb
    noprime!(up)
    up = prime(up,noncommoninds(inds(up),[physical_legs.E,physical_legs.B]))
    updown = up*dag(se*ve*ub*sb)
    uep = prime(ue,noncommoninds(inds(ue),[ancilla_legs.E]))
    vbp = prime(vb,noncommoninds(inds(vb),[ancilla_legs.B]))
    obs.EB = (uep*dag(ue))*updown*(vbp*dag(vb))
    obs.EB = obs.EB*cx.ED*cx.FE*dag(cx.BC)*dag(cx.AB)




    index_ed = commonind(tensA.E, tensA.D)
    index_cd = commonind(tensA.C, tensA.D)
    index_da = commonind(tensA.D, tensA.A)
    ud,sd,vd = svd(tensA.D, [index_ed, index_cd, ancilla_legs.D])
    ua,sa,va = svd(tensA.A, [index_da, physical_legs.A])
    up = sd*vd*hex1.DA*ua*sa
    noprime!(up)
    up = prime(up,noncommoninds(inds(up),[physical_legs.D,physical_legs.A]))
    updown = up*dag(sd*vd*ua*sa)
    udp = prime(ud,noncommoninds(inds(ud),[ancilla_legs.D]))
    vap = prime(va,noncommoninds(inds(va),[ancilla_legs.A]))
    obs.DA = (udp*dag(ud))*updown*(vap*dag(va))
    obs.DA = obs.DA*cx.AF*cx.AB*dag(cx.CD)*dag(cx.ED)


    return obs


end