

function give_tensaA(Gamma::lattice,lambda,physical_legs,ancilla_legs)


    tensA = lattice("tensors")
    tensa = lattice("tensors")

    tensA.A = Gamma.A*sqrt_lambda(lambda.DA)*sqrt_lambda(lambda.AF)*sqrt_lambda(lambda.AB) # inds : ixd izf iyb
    tensA.C = Gamma.C*sqrt_lambda(lambda.BC)*sqrt_lambda(lambda.CD)*sqrt_lambda(lambda.CF) 
    tensA.E = Gamma.E*sqrt_lambda(lambda.FE)*sqrt_lambda(lambda.EB)*sqrt_lambda(lambda.ED) 

    tensA.D = Gamma.D*sqrt_lambda(lambda.CD)*sqrt_lambda(lambda.DA)*sqrt_lambda(lambda.ED) 
    tensA.D = tensA.D*delta(dag(commonind(lambda.ED,Gamma.E)),dag(commonind(lambda.ED,Gamma.D)))
    tensA.D = tensA.D*delta(dag(commonind(lambda.CD,Gamma.C)),dag(commonind(lambda.CD,Gamma.D)))
    tensA.D = tensA.D*delta(dag(commonind(lambda.DA,Gamma.A)),dag(commonind(lambda.DA,Gamma.D)))

    tensA.B = Gamma.B*sqrt_lambda(lambda.AB)*sqrt_lambda(lambda.BC)*sqrt_lambda(lambda.EB) 
    tensA.B = tensA.B*delta(dag(commonind(lambda.AB,Gamma.A)),dag(commonind(lambda.AB,Gamma.B)))
    tensA.B = tensA.B*delta(dag(commonind(lambda.BC,Gamma.C)),dag(commonind(lambda.BC,Gamma.B)))
    tensA.B = tensA.B*delta(dag(commonind(lambda.EB,Gamma.E)),dag(commonind(lambda.EB,Gamma.B)))

    tensA.F = Gamma.F*sqrt_lambda(lambda.CF)*sqrt_lambda(lambda.FE)*sqrt_lambda(lambda.AF) 
    tensA.F = tensA.F*delta(dag(commonind(lambda.CF,Gamma.C)),dag(commonind(lambda.CF,Gamma.F)))
    tensA.F = tensA.F*delta(dag(commonind(lambda.FE,Gamma.E)),dag(commonind(lambda.FE,Gamma.F)))
    tensA.F = tensA.F*delta(dag(commonind(lambda.AF,Gamma.A)),dag(commonind(lambda.AF,Gamma.F)))


    iaf = commonind(tensA.A,tensA.F)
    iab = commonind(tensA.A,tensA.B)
    iad = commonind(tensA.A,tensA.D)

    icb = commonind(tensA.C,tensA.B)
    icd = commonind(tensA.C,tensA.D)
    icf = commonind(tensA.C,tensA.F)

    ief = commonind(tensA.E,tensA.F)
    ied = commonind(tensA.E,tensA.D)
    ieb = commonind(tensA.E,tensA.B)

    tensa.A = prime(tensA.A,noncommoninds(inds(tensA.A),physical_legs.A, ancilla_legs.A))*dag(tensA.A)
    tensa.B = prime(tensA.B,noncommoninds(inds(tensA.B),physical_legs.B, ancilla_legs.B))*dag(tensA.B)
    tensa.C = prime(tensA.C,noncommoninds(inds(tensA.C),physical_legs.C, ancilla_legs.C))*dag(tensA.C)
    tensa.D = prime(tensA.D,noncommoninds(inds(tensA.D),physical_legs.D, ancilla_legs.D))*dag(tensA.D)
    tensa.E = prime(tensA.E,noncommoninds(inds(tensA.E),physical_legs.E, ancilla_legs.E))*dag(tensA.E)
    tensa.F = prime(tensA.F,noncommoninds(inds(tensA.F),physical_legs.F, ancilla_legs.F))*dag(tensA.F)

    cx = double_lattice("tensors")

    cx.AF = combiner(iaf',dag(iaf); tags = "index AF")
    cx.DA = combiner(iad',dag(iad); tags = "index AD")
    cx.AB = combiner(iab',dag(iab); tags = "index AB")

    cx.BC = combiner(icb',dag(icb); tags = "index CB")
    cx.CD = combiner(icd',dag(icd); tags = "index CD")
    cx.CF = combiner(icf',dag(icf); tags = "index CF")

    cx.FE = combiner(ief',dag(ief); tags = "index EF")
    cx.ED = combiner(ied',dag(ied); tags = "index ED")
    cx.EB = combiner(ieb',dag(ieb); tags = "index EB")

    tensa.A = tensa.A*cx.AF*cx.AB*cx.DA
    tensa.C = tensa.C*cx.BC*cx.CD*cx.CF
    tensa.E = tensa.E*cx.FE*cx.ED*cx.EB

    tensa.F = tensa.F*dag(cx.AF)*dag(cx.CF)*dag(cx.FE)
    tensa.B = tensa.B*dag(cx.AB)*dag(cx.BC)*dag(cx.EB)
    tensa.D = tensa.D*dag(cx.DA)*dag(cx.CD)*dag(cx.ED)

    return tensa, tensA, cx

end
