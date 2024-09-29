


function new_qp(environment, p, q, trotter, physlegA, physlegB)

    p0 = deepcopy(p)
    q0 = deepcopy(q)
    q0_prime = prime(q0, noncommoninds(inds(q0), [physlegB]))
    p0_prime = prime(p0, noncommoninds(inds(p0), [physlegA]))

    ida = ITensor(physlegA, dag(physlegA)')
    idb = ITensor(physlegB, dag(physlegB)')
    for ii = 1:8
        ida[physlegA=>ii, dag(physlegA)' =>ii] = 1;
        idb[physlegB=>ii, dag(physlegB)' =>ii] = 1;
    end

    distance = 1
    err_fullupdate = 1;
    precision_full_update = 1e-7

    q_prime = prime(q, noncommoninds(inds(q), [physlegB]))
    p_prime = prime(p, noncommoninds(inds(p), [physlegA]))

    bb = ((environment*p0_prime*q0_prime))*trotter;

    itx = 0

    while distance > precision_full_update && itx < 50

        itx = itx + 1

        q_prime = prime(q, noncommoninds(inds(q),[physlegB]))
        p_prime = prime(p, noncommoninds(inds(p), [physlegA]))

        M = ((environment*q_prime)*dag(q))*ida
                
        b = noprime(bb,physlegB')
        b = noprime(b*dag(q))

        indd1 = commoninds(b,M) 
        c = combiner(indd1, dir = dir(physlegA))
        indd2 = noncommoninds(inds(M), indd1)
        cprime = combiner([indd2[1], indd2[3], indd2[2]])
        # cprime = combiner([indd2[2], indd2[3], indd2[1]])

        M_matrix = M*c*cprime
        # @show inds(M_matrix)
        inverse_M_matrix = inverse_matrixM(M_matrix)
        inverse_M = inverse_M_matrix*(c)*(cprime)

        p = inverse_M*b
        p = noprime(p)

        q_prime = prime(q, noncommoninds(inds(q), [physlegB]))
        p_prime = prime(p, noncommoninds(inds(p), [physlegA]))

        M = (environment*p_prime)*dag(p)*idb
        b = noprime(bb,physlegA')
        b = noprime(b*dag(p))

        indd1 = commoninds(b,M)
        c = combiner(indd1, dir = -dir(physlegB))
        indd2 = noncommoninds(inds(M), indd1)
        cprime = combiner([indd2[1],indd2[3],indd2[2]])
        # @show indd1
        # @show indd2 # okok !!

        M_matrix = M*c*cprime
        # @show inds(M_matrix)
        inverse_M_matrix = inverse_matrixM(M_matrix)
        inverse_M = inverse_M_matrix*(c)*(cprime)

        q = inverse_M*b
        q = noprime(q)
            
        q_prime = prime(q, noncommoninds(inds(q),[physlegB]))
        p_prime = prime(p, noncommoninds(inds(p),[physlegA]))

        middle_term = p0_prime*q0_prime*trotter
        middle_term2 = noprime(middle_term,physlegA',physlegB')*trotter

        term1 = noprime(environment*middle_term2)*dag(p0)*dag(q0)
        term2 = environment*(q_prime*dag(q))*(p_prime*dag(p))
        term3 = noprime(environment*middle_term)*dag(p)*dag(q)

        Z = norm(array(((environment*q0_prime)*dag(q0))*p0_prime*dag(p0)))

        previous_err_fulluptade = deepcopy(err_fullupdate)
        err_fullupdate = norm(array(term1) + array(term2) - 2*array(term3))/Z
        @show distance = norm(err_fullupdate - previous_err_fulluptade)

    end

    p = p/norm(((p)))
    q = q/norm(((q)))

    return p, q

end

function new_qp2(environment, p, q, trotter, physlegA, physlegB)

    p0 = deepcopy(p)
    q0 = deepcopy(q)
    q0_prime = prime(q0, noncommoninds(inds(q0), [physlegB]))
    p0_prime = prime(p0, noncommoninds(inds(p0), [physlegA]))

    ida = ITensor(physlegA, dag(physlegA)')
    idb = ITensor(physlegB, dag(physlegB)')
    for ii = 1:8
        ida[physlegA=>ii, dag(physlegA)' =>ii] = 1;
        idb[physlegB=>ii, dag(physlegB)' =>ii] = 1;
    end

    distance = 1
    err_fullupdate = 1;
    precision_full_update = 1e-7

    q_prime = prime(q, noncommoninds(inds(q), [physlegB]))
    p_prime = prime(p, noncommoninds(inds(p), [physlegA]))

    bb = ((environment*p0_prime*q0_prime))*trotter;

    itx = 0

    while distance > precision_full_update && itx < 50

        itx = itx + 1

        q_prime = prime(q, noncommoninds(inds(q),[physlegB]))
        p_prime = prime(p, noncommoninds(inds(p), [physlegA]))

        M = ((environment*q_prime)*dag(q))*ida
                
        b = noprime(bb,physlegB')
        b = noprime(b*dag(q))

        indd1 = commoninds(b,M) 
        c = combiner(indd1, dir = -dir(physlegA))
        indd2 = noncommoninds(inds(M), indd1)
        cprime = combiner([indd2[1], indd2[3], indd2[2]])
        # cprime = combiner([indd2[2], indd2[3], indd2[1]])

        M_matrix = M*c*cprime
        # @show inds(M_matrix)
        inverse_M_matrix = inverse_matrixM(M_matrix)
        inverse_M = inverse_M_matrix*(c)*(cprime)

        p = inverse_M*b
        p = noprime(p)

        q_prime = prime(q, noncommoninds(inds(q), [physlegB]))
        p_prime = prime(p, noncommoninds(inds(p), [physlegA]))

        M = (environment*p_prime)*dag(p)*idb
        b = noprime(bb,physlegA')
        b = noprime(b*dag(p))

        indd1 = commoninds(b,M)
        c = combiner(indd1, dir = dir(physlegB))
        indd2 = noncommoninds(inds(M), indd1)
        cprime = combiner([indd2[1],indd2[3],indd2[2]])
        # @show indd1
        # @show indd2 # okok !!

        M_matrix = M*c*cprime
        # @show inds(M_matrix)
        inverse_M_matrix = inverse_matrixM(M_matrix)
        inverse_M = inverse_M_matrix*(c)*(cprime)

        q = inverse_M*b
        q = noprime(q)
            
        q_prime = prime(q, noncommoninds(inds(q),[physlegB]))
        p_prime = prime(p, noncommoninds(inds(p),[physlegA]))

        middle_term = p0_prime*q0_prime*trotter
        middle_term2 = noprime(middle_term,physlegA',physlegB')*trotter

        term1 = noprime(environment*middle_term2)*dag(p0)*dag(q0)
        term2 = environment*(q_prime*dag(q))*(p_prime*dag(p))
        term3 = noprime(environment*middle_term)*dag(p)*dag(q)

        Z = norm(array(((environment*q0_prime)*dag(q0))*p0_prime*dag(p0)))

        previous_err_fulluptade = deepcopy(err_fullupdate)
        err_fullupdate = norm(array(term1) + array(term2) - 2*array(term3))/Z
        @show distance = norm(err_fullupdate - previous_err_fulluptade)

    end

    p = p/norm(((p)))
    q = q/norm(((q)))

    return p, q

end