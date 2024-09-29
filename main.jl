
using ITensors
using Random
using MAT
using JLD2
using FileIO
using MKL
using LinearAlgebra
using MAT


# D = parse(Int, ARGS[1])
# h = parse(Float32, ARGS[2])



include("src/cool_system_down.jl")


BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)
ITensors.disable_threaded_blocksparse()


println("BEGIN OF SIMULATION")



let 
   
    D = 3
    h = 1
    J1 = 1
    J2 = -0.201;
    J3 = -0.204;
    J4 = 0.278;
    J5 = 0.402

    temp_max = 9.99
    temp_min = 8.01
    number_of_points = 10

    temperature = LinRange(temp_max,temp_min,number_of_points)
    cool_system_down(J1,J2,J3,J4,J5,h,D,temperature)

end




