





mutable struct lattice

    A
    B
    C
    D
    E
    F

    function lattice(A,B,C,D,E,F)
        new(A,B,C,D,E,F)
    end

    function lattice(str)
        if str == "tensors"
            new(ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor())

        elseif str == "index"
     
            new(
            Index([QN(-3)=>1,QN(-1)=>2,QN(1)=>1,QN(-1)=>1,QN(1)=>2,QN(3)=>1]; tags="ia"),
            Index([QN(-3)=>1,QN(-1)=>2,QN(1)=>1,QN(-1)=>1,QN(1)=>2,QN(3)=>1]; tags="ib"),
            Index([QN(-3)=>1,QN(-1)=>2,QN(1)=>1,QN(-1)=>1,QN(1)=>2,QN(3)=>1]; tags="ic"),
            Index([QN(-3)=>1,QN(-1)=>2,QN(1)=>1,QN(-1)=>1,QN(1)=>2,QN(3)=>1]; tags="id"),
            Index([QN(-3)=>1,QN(-1)=>2,QN(1)=>1,QN(-1)=>1,QN(1)=>2,QN(3)=>1]; tags="ie"),
            Index([QN(-3)=>1,QN(-1)=>2,QN(1)=>1,QN(-1)=>1,QN(1)=>2,QN(3)=>1]; tags="if"))
            
            
        end
    end        

end

mutable struct double_lattice
    AB
    CF
    ED
    DA
    BC
    FE
    AF
    CD
    EB

    function double_lattice(str)
        if str == "tensors"
        new(ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor())
        end
    end
    function double_lattice(A,B,C,D,E,F,G,H,U)
        new(A,B,C,D,E,F,G,H,U)
    end

    
end


include("simpleupdate.jl")
include("evolution_ctm.jl")
include("get_ener_tensors_J1-5.jl")
include("energy.jl")
include("Initialisation.jl")
include("fullupdate/full_update.jl")
include("get_tensaA.jl")
include("diag_tensors.jl")

function cool_system_down(J1,J2,J3,J4,J5,h,D,temperature)

dbeta = 1e-2

final_temp = temperature[end] - 1e-6

Gamma, ix, iy, iz, lambda, physical_legs, ancilla_legs, temp, energy_array, temperature_array, magne_array  = Initialisation(J1,J2,J3,J4,J5,h,D)
trotter_gate, ener_tens, hex1, hex2  = get_ener_tensors_J15(physical_legs,J1,J2,J3,J4,J5,h,dbeta)

if temperature[1] > temp
    println("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    println("ERROR target temperature is higher then the actual temperature")
    println("TARGET TEMPERATURE : "*string(temperature[1]))
    println("ACTUAL TEMPEARTURE : "*string(temp))
    println("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")

    return 
end

modit = 10
ener = 7; enertmp = 9
it = 0

iteration = []
cx = lattice("tensors")
tensa = lattice("tensors")
tensA = lattice("tensors")

SS1 = double_lattice("tensors")
SS2 = double_lattice("tensors")
magne = lattice("tensors")
all_magne_ = [lattice("tensors"),lattice("tensors"),lattice("tensors")]

beta = 1/temp
relevant = 1

while final_temp < temp

    tempprev = 1/beta
    beta = beta + dbeta; 
    temp = 1/beta;
    it += 1
    println(temp)

    Gamma, lambda, ix, iy, iz = simpleupdate(Gamma, lambda, ix, iy, iz, physical_legs, ancilla_legs, trotter_gate, D)

    if size(temperature)[1] >= relevant && size(temperature)[1] > 0 
        if (temperature[relevant] <= tempprev && temperature[relevant] > temp) 

        relevant = relevant + 1    
        enertmp = ener
        tensa, tensA, cx = give_tensaA(Gamma,lambda,physical_legs,ancilla_legs)

        @showtime ener, T1, T2, T3, C1ab,C1ed,C1cf,C2cd,C2af,C2eb,C3fe,C3bc,C3da,magne = 
        energy(tensa,tensA,cx,D,physical_legs, ancilla_legs, ener_tens, hex1, hex2)

        @show ener

        push!(energy_array, ener)
        push!(temperature_array, temp)
        push!(magne_array, magne)

        end
    end



end

output = "Results"

save(output*"/LocalTensors_J1$(J1)_J2$(J2)_J3$(J3)_J4$(J4)_J5$(J5)_D$(D)_h$(h).jld2",
"Gamma",Gamma,
"lambda",lambda,
"ix",ix,
"iy",iy,
"iz",iz,
"physical_legs",physical_legs, 
"ancilla_legs",ancilla_legs,
"temp", temp,
"dbeta", dbeta,
"trotter_gate",trotter_gate,
"ener",enertmp,
"SS1", SS1,
"SS2", SS2,
"magne", magne,
"all_magne_",all_magne_,
"energy_array", energy_array,
"temperature_array",temperature_array,
"magne_array", magne_array,
"J1",J1,
"J2",J2,
"J3",J3,
"J4",J4,
"J5",J5,
"h",h
)

file = matopen(output*"/LocalTensors_J1$(J1)_J2$(J2)_J3$(J3)_J4$(J4)_J5$(J5)_D$(D)_h$(h).mat", "w")
write(file, "ener", energy_array)
write(file, "temp", temperature_array)
write(file, "dbeta", dbeta)
write(file, "magne", magne_array)
write(file, "D", D)
write(file, "J1", J1)
write(file, "J2", J2)
write(file, "J3", J3)
write(file, "J4", J4)
write(file, "J5", J5)
write(file, "h", h)
close(file)

end


