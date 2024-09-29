

function get_ener_tensors_J15(physical_legs,J1,J2,J3,J4,J5,hfield,dbetasu)

    iia = physical_legs.A
    iib = physical_legs.B
    iic = physical_legs.C
    iid = physical_legs.D
    iie = physical_legs.E
    iif = physical_legs.F

    
    h, ener, hex, singlehex = [double_lattice("tensors") for _ in 1:4]

    h.AB = ITensor(dag(iia),dag(iib),iia',iib')
    h.CF = ITensor(dag(iic),dag(iif),iic',iif')
    h.ED = ITensor(dag(iie),dag(iid),iie',iid')
    h.DA = ITensor(dag(iid),dag(iia),iid',iia')
    h.BC = ITensor(dag(iib),dag(iic),iib',iic')
    h.FE = ITensor(dag(iif),dag(iie),iif',iie')
    h.AF = ITensor(dag(iia),dag(iif),iia',iif')
    h.CD = ITensor(dag(iic),dag(iid),iic',iid')
    h.EB = ITensor(dag(iie),dag(iib),iie',iib')
    
    ener.AB = ITensor(dag(iia),dag(iib),iia',iib')
    ener.ED = ITensor(dag(iie),dag(iid),iie',iid')
    ener.CF = ITensor(dag(iic),dag(iif),iic',iif')
    ener.DA = ITensor(dag(iid),dag(iia),iid',iia')
    ener.BC = ITensor(dag(iib),dag(iic),iib',iic')
    ener.FE = ITensor(dag(iif),dag(iie),iif',iie')
    ener.AF = ITensor(dag(iia),dag(iif),iia',iif')
    ener.CD = ITensor(dag(iic),dag(iid),iic',iid')
    ener.EB = ITensor(dag(iie),dag(iib),iie',iib')

    hex.AB = ITensor(dag(iia),dag(iib),iia',iib')
    hex.ED = ITensor(dag(iie),dag(iid),iie',iid')
    hex.CF = ITensor(dag(iic),dag(iif),iic',iif')
    hex.DA = ITensor(dag(iid),dag(iia),iid',iia')
    hex.BC = ITensor(dag(iib),dag(iic),iib',iic')
    hex.FE = ITensor(dag(iif),dag(iie),iif',iie')
    hex.AF = ITensor(dag(iia),dag(iif),iia',iif')
    hex.CD = ITensor(dag(iic),dag(iid),iic',iid')
    hex.EB = ITensor(dag(iie),dag(iib),iie',iib')
    
    singlehex.AB = ITensor(dag(iia),dag(iib),iia',iib')

    hex1 = double_lattice("tensors")
    hex2 = double_lattice("tensors")

    hex1.AB = ITensor(dag(iia),dag(iib),iia',iib')
    hex1.ED = ITensor(dag(iie),dag(iid),iie',iid')
    hex1.CF = ITensor(dag(iic),dag(iif),iic',iif')
    hex1.DA = ITensor(dag(iid),dag(iia),iid',iia')
    hex1.BC = ITensor(dag(iib),dag(iic),iib',iic')
    hex1.FE = ITensor(dag(iif),dag(iie),iif',iie')
    hex1.AF = ITensor(dag(iia),dag(iif),iia',iif')
    hex1.CD = ITensor(dag(iic),dag(iid),iic',iid')
    hex1.EB = ITensor(dag(iie),dag(iib),iie',iib')

    hex2.AB = ITensor(dag(iia),dag(iib),iia',iib')
    hex2.ED = ITensor(dag(iie),dag(iid),iie',iid')
    hex2.CF = ITensor(dag(iic),dag(iif),iic',iif')
    hex2.DA = ITensor(dag(iid),dag(iia),iid',iia')
    hex2.BC = ITensor(dag(iib),dag(iic),iib',iic')
    hex2.FE = ITensor(dag(iif),dag(iie),iif',iie')
    hex2.AF = ITensor(dag(iia),dag(iif),iia',iif')
    hex2.CD = ITensor(dag(iic),dag(iid),iic',iid')
    hex2.EB = ITensor(dag(iie),dag(iib),iie',iib')


    sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
    S = [sx,sy,sz]

    hy, hx, hz = [zeros(64,64) for _ in 1:3]
    
    hex1array = double_lattice(zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64))
    hex2array = double_lattice(zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64),zeros(64,64))
    


    for i = 1:3

        # AB
        hy =  hy + 
        J2*kron(id,S[i],id,S[i],id,id) + # 24 J2
        J3*kron(id,id,S[i],id,id,S[i]) + # 36 J3
        J5*kron(id,S[i],S[i],id,id,id) + # 23 J5
        J4*kron(id,id,id,S[i],id,S[i]) + # 46 J4
        J1*kron(id,S[i],id,id,id,S[i])   # 26 J1

        # BC
        hx = hx + 
        J3*kron(id,S[i],id,id,S[i],id) + # 25 J3
        J2*kron(id,id,S[i],S[i],id,id) + # 34 J2
        J4*kron(id,S[i],S[i],id,id,id) + # 23 J4
        J5*kron(id,id,id,S[i],S[i],id) + # 45 J5
        J1*kron(id,S[i],id,S[i],id,id)  # 24 J1
    
        # CD
        hz = hz + 
        J3*kron(S[i],id,id,S[i],id,id) + # 14 J3
        J2*kron(id,id,S[i],id,S[i],id) + # 35 J2
        J5*kron(S[i],id,S[i],id,id,id) + # 13 J5
        J4*kron(id,id,id,S[i],S[i],id) + # 45 J4
        J1*kron(id,id,S[i],S[i],id,id)   # 34 J1
  
    
    end

    hy =  hy + hfield/2*(
        kron(id,sz,id,id,id,id) + # 2436
        kron(id,id,id,sz,id,id) + 
        kron(id,id,sz,id,id,id) + 
        kron(id,id,id,id,id,sz))

    hx = hx +  hfield/2*(
        kron(id,sz,id,id,id,id) + # 2534
        kron(id,id,id,id,sz,id) + 
        kron(id,id,sz,id,id,id) + 
        kron(id,id,id,sz,id,id))

    hz = hz + hfield/2*(
        kron(sz,id,id,id,id,id) + # 1435
        kron(id,id,id,sz,id,id) + 
        kron(id,id,sz,id,id,id) + 
        kron(id,id,id,id,sz,id))

    gx = exp(-hx*dbetasu/2)
    gx = reshape(gx, (8,8,8,8))
    gy = exp(-hy*dbetasu/2)
    gy = reshape(gy, (8,8,8,8))
    gz = exp(-hz*dbetasu/2)
    gz = reshape(gz, (8,8,8,8))
    
    hx = reshape(hx, (8,8,8,8))
    hy = reshape(hy, (8,8,8,8))
    hz = reshape(hz, (8,8,8,8))
    


    for i = 1:3

        # Spin-Spin correlation whithin the inner and outer honeycomb 

        hex1array.AB = hex1array.AB + kron(id,S[i],id,S[i],id,id) # 24
        hex1array.ED = hex1array.ED + kron(id,S[i],id,S[i],id,id) # 24
        hex1array.CF = hex1array.CF + kron(id,S[i],id,S[i],id,id) # 24

        hex1array.AF = hex1array.AF + kron(S[i],id,id,S[i],id,id) # 14
        hex1array.CD = hex1array.CD + kron(S[i],id,id,S[i],id,id) # 14
        hex1array.EB = hex1array.EB + kron(S[i],id,id,S[i],id,id) # 14
        
        hex1array.DA = hex1array.DA + kron(id,S[i],id,id,S[i],id) # 25
        hex1array.BC = hex1array.BC + kron(id,S[i],id,id,S[i],id) # 25
        hex1array.FE = hex1array.FE + kron(id,S[i],id,id,S[i],id) # 25

        hex2array.AB = hex2array.AB + kron(id,id,S[i],id,id,S[i]) # 36
        hex2array.ED = hex2array.ED + kron(id,id,S[i],id,id,S[i]) # 36
        hex2array.CF = hex2array.CF + kron(id,id,S[i],id,id,S[i]) # 36
        
        hex2array.AF = hex2array.AF + kron(id,id,S[i],id,S[i],id) # 35
        hex2array.CD = hex2array.CD + kron(id,id,S[i],id,S[i],id) # 35
        hex2array.EB = hex2array.EB + kron(id,id,S[i],id,S[i],id) # 35
        
        hex2array.DA = hex2array.DA + kron(id,id,S[i],S[i],id,id) # 34
        hex2array.BC = hex2array.BC + kron(id,id,S[i],S[i],id,id) # 34
        hex2array.FE = hex2array.FE + kron(id,id,S[i],S[i],id,id) # 34

    end

    
    for name in fieldnames(double_lattice)
        setproperty!(hex1array, name, reshape(getproperty(hex1array,name), (8,8,8,8)))
        setproperty!(hex2array, name, reshape(getproperty(hex2array,name), (8,8,8,8)))
    end

    charges = [-3, -1, -1, 1, -1, 1, 1, 3]
   
    
    for i1 = 1:8, i2=1:8, i3=1:8, i4=1:8
    
        if norm(charges[i1] +  charges[i2] -  charges[i3] -  charges[i4]) < 1e-9
                        
            h.AB[iia=>i1, iib =>i2, iia' => i3, iib' =>i4] = gy[i1,i2,i3,i4]
            h.ED[iie=>i1, iid =>i2, iie' => i3, iid' =>i4] = gy[i1,i2,i3,i4]
            h.CF[iic=>i1, iif =>i2, iic' => i3, iif' =>i4] = gy[i1,i2,i3,i4]
        
            h.DA[iid=>i1, iia =>i2, iid' => i3, iia' =>i4] = gx[i1,i2,i3,i4]
            h.BC[iib=>i1, iic =>i2, iib' => i3, iic' =>i4] = gx[i1,i2,i3,i4]
            h.FE[iif=>i1, iie =>i2, iif' => i3, iie' =>i4] = gx[i1,i2,i3,i4]
        
            h.CD[iic=>i1, iid =>i2, iic' => i3, iid' =>i4] = gz[i1,i2,i3,i4]
            h.AF[iia=>i1, iif =>i2, iia' => i3, iif' =>i4] = gz[i1,i2,i3,i4]
            h.EB[iie=>i1, iib =>i2, iie' => i3, iib' =>i4] = gz[i1,i2,i3,i4]
        
            ener.DA[iid=>i1, iia =>i2, iid' => i3, iia' =>i4] = hx[i1,i2,i3,i4]
            ener.BC[iib=>i1, iic =>i2, iib' => i3, iic' =>i4] = hx[i1,i2,i3,i4]
            ener.FE[iif=>i1, iie =>i2, iif' => i3, iie' =>i4] = hx[i1,i2,i3,i4]
        
            ener.CD[iic=>i1, iid =>i2, iic' => i3, iid' =>i4] = hz[i1,i2,i3,i4]
            ener.AF[iia=>i1, iif =>i2, iia' => i3, iif' =>i4] = hz[i1,i2,i3,i4]
            ener.EB[iie=>i1, iib =>i2, iie' => i3, iib' =>i4] = hz[i1,i2,i3,i4]

            ener.AB[iia=>i1, iib =>i2, iia' => i3, iib' =>i4] = hy[i1,i2,i3,i4]
            ener.ED[iie=>i1, iid =>i2, iie' => i3, iid' =>i4] = hy[i1,i2,i3,i4]
            ener.CF[iic=>i1, iif =>i2, iic' => i3, iif' =>i4] = hy[i1,i2,i3,i4]
        
 
                        
            hex1.DA[iid=>i1, iia =>i2, iid' => i3, iia' =>i4] = hex1array.DA[i1,i2,i3,i4]
            hex1.BC[iib=>i1, iic =>i2, iib' => i3, iic' =>i4] = hex1array.BC[i1,i2,i3,i4]
            hex1.FE[iif=>i1, iie =>i2, iif' => i3, iie' =>i4] = hex1array.FE[i1,i2,i3,i4]
        
            hex1.CD[iic=>i1, iid =>i2, iic' => i3, iid' =>i4] = hex1array.CD[i1,i2,i3,i4]
            hex1.AF[iia=>i1, iif =>i2, iia' => i3, iif' =>i4] = hex1array.AF[i1,i2,i3,i4]
            hex1.EB[iie=>i1, iib =>i2, iie' => i3, iib' =>i4] = hex1array.EB[i1,i2,i3,i4]

            hex1.AB[iia=>i1, iib =>i2, iia' => i3, iib' =>i4] = hex1array.AB[i1,i2,i3,i4]
            hex1.ED[iie=>i1, iid =>i2, iie' => i3, iid' =>i4] = hex1array.ED[i1,i2,i3,i4]
            hex1.CF[iic=>i1, iif =>i2, iic' => i3, iif' =>i4] = hex1array.CF[i1,i2,i3,i4]

            hex2.DA[iid=>i1, iia =>i2, iid' => i3, iia' =>i4] = hex2array.DA[i1,i2,i3,i4]
            hex2.BC[iib=>i1, iic =>i2, iib' => i3, iic' =>i4] = hex2array.BC[i1,i2,i3,i4]
            hex2.FE[iif=>i1, iie =>i2, iif' => i3, iie' =>i4] = hex2array.FE[i1,i2,i3,i4]
    
            hex2.CD[iic=>i1, iid =>i2, iic' => i3, iid' =>i4] = hex2array.CD[i1,i2,i3,i4]
            hex2.AF[iia=>i1, iif =>i2, iia' => i3, iif' =>i4] = hex2array.AF[i1,i2,i3,i4]
            hex2.EB[iie=>i1, iib =>i2, iie' => i3, iib' =>i4] = hex2array.EB[i1,i2,i3,i4]

            hex2.AB[iia=>i1, iib =>i2, iia' => i3, iib' =>i4] = hex2array.AB[i1,i2,i3,i4]
            hex2.ED[iie=>i1, iid =>i2, iie' => i3, iid' =>i4] = hex2array.ED[i1,i2,i3,i4]
            hex2.CF[iic=>i1, iif =>i2, iic' => i3, iif' =>i4] = hex2array.CF[i1,i2,i3,i4]

        end
    
    end




    return h, ener, hex1, hex2

end

