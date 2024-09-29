

function renormalise_tensor_to_norm!(vec::ITensor...)

    for v in vec
        v=v/norm(v)
    end

    return nothing
end

function renormalise_lattice_to_norm!(Ts::lattice...)

    for T in Ts
        for field in fieldnames(lattice)
            setfield!(T,field, getfield(T,field)/norm(getfield(T,field)))
        end
    end

    return nothing
end

function renormalise_tensor_to_one!(vec::ITensor...)

    for v in vec
        v.=v./v
    end

    return nothing
end

function renormalise_lattice_to_one!(Ts::lattice...)

    for T in Ts
        T.A .= T.A./T.A
        T.B .= T.B./T.B
        T.C .= T.C./T.C
        T.D .= T.D./T.D
        T.E .= T.E./T.E
        T.F .= T.F./T.F
    end

    return nothing
end

