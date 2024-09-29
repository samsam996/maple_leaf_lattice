

include("diag_tensors.jl")

function isometryCorboz(C1cf,C2eb,chi)

    rho1 = C1cf*C2eb
    u,s,v = svd(rho1, commoninds(C1cf*C2eb,C2eb), maxdim = chi)
    Pbar = (C1cf*dag(v))*inv_sqrt_lambda(s)
    P = (C2eb*dag(u))
    P = P*inv_sqrt_lambda(s)
    indv = commonind(s,v)
    indu = commonind(u,s)
    P = P*delta(indv,dag(indu))
    
    return P, Pbar

end