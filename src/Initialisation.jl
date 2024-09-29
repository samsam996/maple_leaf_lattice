


function Initialisation(J1,J2,J3,J4,J5,h,D)


    energy_array = []
    temperature_array = []
    magne_array = []
    temp = []
    
    Gamma = lattice("tensors")
    lambda = double_lattice("tensors")


    physical_legs, ancilla_legs, ix, iy, iz = [lattice("index") for _ in 1:5]
    
    file_name = "Results/LocalTensors_J1$(J1)_J2$(J2)_J3$(J3)_J4$(J4)_J5$(J5)_D$(D)_h$(h).jld2"

    if isfile(file_name)

        Gamma, ix, iy, iz, lambda, physical_legs, ancilla_legs,temp, energy_array, temperature_array, magne_array = 
        load(file_name,"Gamma","ix","iy","iz","lambda","physical_legs","ancilla_legs","temp",
        "energy_array",
        "temperature_array",
        "magne_array")

        println("ALREADY DONE UNTIL TEMP:"*string(temp))
        
    else
        
        ix.B = Index([QN(0)=>1];tags="ixb"); iy.B = Index([QN(0)=>1];tags ="iyb"); iz.B = Index([QN(0)=>1];tags= "izb"); 
        ix.A = Index([QN(0)=>1];tags="ixa"); iy.A = Index([QN(0)=>1];tags= "iya"); iz.A = Index([QN(0)=>1];tags= "iza");
        ix.D = Index([QN(0)=>1];tags="ixd"); iy.D = Index([QN(0)=>1];tags= "iyd"); iz.D = Index([QN(0)=>1];tags= "izd"); 
        ix.C = Index([QN(0)=>1];tags="ixc"); iy.C = Index([QN(0)=>1];tags= "iyc"); iz.C = Index([QN(0)=>1];tags= "izc"); 
        ix.E = Index([QN(0)=>1];tags="ixe"); iy.E = Index([QN(0)=>1];tags= "iye"); iz.E = Index([QN(0)=>1];tags= "ize"); 
        ix.F = Index([QN(0)=>1];tags="ixf"); iy.F = Index([QN(0)=>1];tags= "iyf"); iz.F = Index([QN(0)=>1];tags= "izf"); 


        Gamma = lattice(
            ITensor(dag(ix.A),iy.A,iz.A,physical_legs.A,dag(ancilla_legs.A)),
            ITensor(ix.B,dag(iy.B),dag(iz.B),physical_legs.B,dag(ancilla_legs.B)),
            ITensor(dag(ix.C),iy.C,iz.C,physical_legs.C,dag(ancilla_legs.C)),
            ITensor(ix.D,dag(iy.D),dag(iz.D),physical_legs.D,dag(ancilla_legs.D)),
            ITensor(dag(ix.E),iy.E,iz.E,physical_legs.E,dag(ancilla_legs.E)),
            ITensor(ix.F,dag(iy.F),dag(iz.F),physical_legs.F,dag(ancilla_legs.F)))
            
        for i = 1:8
            Gamma.A[ix.A=>1,iy.A=>1,iz.A=>1,physical_legs.A=>i,ancilla_legs.A=>i] = 1.
            Gamma.B[ix.B=>1,iy.B=>1,iz.B=>1,physical_legs.B=>i,ancilla_legs.B=>i] = 1.
            Gamma.C[ix.C=>1,iy.C=>1,iz.C=>1,physical_legs.C=>i,ancilla_legs.C=>i] = 1.
            Gamma.D[ix.D=>1,iy.D=>1,iz.D=>1,physical_legs.D=>i,ancilla_legs.D=>i] = 1.
            Gamma.E[ix.E=>1,iy.E=>1,iz.E=>1,physical_legs.E=>i,ancilla_legs.E=>i] = 1.
            Gamma.F[ix.F=>1,iy.F=>1,iz.F=>1,physical_legs.F=>i,ancilla_legs.F=>i] = 1.
        end

        lambda.AB = abs.(ITensor(dag(iy.A),(iy.B))) # y 
        lambda.CF = abs.(ITensor(dag(iy.C),(iy.F))) # y 
        lambda.ED = abs.(ITensor(dag(iy.E),(iy.D))) # y 
        lambda.DA = abs.(ITensor(dag(ix.D),(ix.A))) # x
        lambda.BC = abs.(ITensor(dag(ix.B),(ix.C))) # x
        lambda.FE = abs.(ITensor(dag(ix.F),(ix.E))) # x
        lambda.AF = abs.(ITensor(dag(iz.A),(iz.F))) # z 
        lambda.CD = abs.(ITensor(dag(iz.C),(iz.D))) # z
        lambda.EB = abs.(ITensor(dag(iz.E),(iz.B))) # z

        lambda.AB[iy.A=>1, iy.B=>1] = 1.
        lambda.CF[iy.C=>1, iy.F=>1] = 1.
        lambda.ED[iy.E=>1, iy.D=>1] = 1.
        lambda.DA[ix.D=>1, ix.A=>1] = 1.
        lambda.BC[ix.B=>1, ix.C=>1] = 1.
        lambda.FE[ix.F=>1, ix.E=>1] = 1.
        lambda.AF[iz.A=>1, iz.F=>1] = 1.
        lambda.CD[iz.C=>1, iz.D=>1] = 1.
        lambda.EB[iz.E=>1, iz.B=>1] = 1.

        temp = Inf

    end
    

    return Gamma, ix, iy, iz, lambda, physical_legs, ancilla_legs, temp, energy_array, temperature_array, magne_array
    

end