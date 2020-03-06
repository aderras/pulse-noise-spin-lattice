# This module contains functions that call the numerical solver. Useful
# function is runComputation(), which runs the runge kutta solver and 
# saves the resulting data. 
module computeSolution

    using modifyFiles,calculateEnergy,calculateTopologicalCharge,LLequation,
    calculateMagnetization,rungeKutta,LLequation,normalize,noiseRotation,Dates
    export evaluateLL!,runRelaxation!,pulseNoiseStep!

    # The following function evaluates the LLG equation. In the comments you'll find
    # two options: pulseNoise or rk4. For zero temperature relaxation, it doesn't matter
    # which you use. For nonzero temperature, must use pulse noise. 
    #
    # inputs: mat = (N,N,3) initial condition, topChargeArray = list of topological
    # charges during computation, energyArray = list of energy values, reldir = the
    # directory where data is stored, relative to pwd(), Hext = current external field,
    # maxLoop = max iterations
    # outputs: nothing. The resulting spin field is stored in mat
    function evaluateLL!(mat::Array{Float64,3}, matParams::Array{Float64,1},
        evalParams::Array{Float64,1})

        reldir = "/data"
        timestampString = string(Dates.format(Dates.now(),"HH_MM_SS"),
                trunc(Int,rand()*10000))

        j, h, dmi, pma, ed, pbc = matParams
        tRK, nSteps, damping, t = evalParams

        maxLoop = 100000   
        
        enArray = zeros(maxLoop)
        topChargeArray = zeros(maxLoop)
        reffArray = zeros(maxLoop)


        for i in 1:maxLoop
        
            enArray[i]          = calcEnergy(mat, matParams)
	        topChargeArray[i]	= calcQ(mat)
	        reffArray[i]	    = sqrt(calcM2(mat)*m*n/(4*pi))

            pulseNoiseStep!(mat, matParams, evalParams)

            # println(topChargeArray[i])
            
            # If the topological charge drops, quit
            if i > 10 && abs(topChargeArray[i]) < 0.1
                break 
            end

        end

        writeDataH5(string(reldir,"/spin_field_after_dynamics_T=",t,"_H=",h,"_DMI=",dmi,
        "_DDI=",ed,"_PMA=",pma,"_",timestampString,"_.h5"),mat)
        writeDataH5(string(reldir,"/topological_charge_during_dynamics_T=",t,"_",
        "H=",h,"_DMI=",dmi,"_DDI=",ed,"_PMA=",pma,"_",timestampString,".h5"),filter(x->x!=0.,topChargeArray))
        writeDataH5(string(reldir,"/energy_during_dynamics_T=",t,"_","H=",h,"_DDI=",
        ed,"_DMI=",dmi,"_PMA=",pma,"_",timestampString,".h5"),filter(x->x!=0.,enArray))
	    writeDataH5(string(reldir,"/effective_size_during_dynamics_T=",t,"_","H=",h,"_DDI=",
        ed,"_DMI=",dmi,"_PMA=",pma,"_",timestampString,".h5"),filter(x->x!=0.,reffArray))

    end

    # runs a single pulse noise step
    function pulseNoiseStep!(s::Array{Float64,3}, matParams::Array{Float64,1}, 
        evalParams::Array{Float64,1})

        # This is the pulse-noise algorithm
        rk4!(s,RHS!, matParams, evalParams)
        rotateSpins!(s, evalParams)
        rk4!(s,RHS!, matParams, evalParams)
        normalizeSpins!(s)

    end

    # first run relaxation on the spin field
    function runRelaxation!(mat, matParams, evalParams)

        maxLoop = 100000

        j, h, dmi, pma, ed, pbc = matParams
        tRK, nSteps, damping, t = evalParams
        p, m, n = size(mat)

        # Only use the following 3 lines when saving data
        enArray = zeros(maxLoop)
        topChargeArray = zeros(maxLoop)
        reffArray = zeros(maxLoop)

        prevEnergy = 0.0
        currEnergy = 0.0
        diff = 10.0

        @time for i in 1:maxLoop

            rk4!(mat,RHS!, matParams, evalParams)
            normalizeSpins!(mat)
            currEnergy = calcEnergy(mat, matParams);
            diff = currEnergy -prevEnergy

            enArray[i]          = currEnergy
            topChargeArray[i]   = calcQ(mat)
            reffArray[i]        = sqrt(calcM2(mat)*m*n/(4*pi))

            #println(string("i = ",i," E = ",enArray[i],", Q = ",topChargeArray[i],", Reff = ",reffArray[i]))
            #sleep(0.5)

            if abs(diff) < 0.00004 || topChargeArray[i] < 0
                break
            else
                prevEnergy = currEnergy
            end
        end

        reldir = "/data"
        writeDataH5(string(reldir,"/spin_field_after_relaxation_T=",t,"_H=",h,
        "_DMI=",dmi,"_DDI=",ed,"_PMA=",pma,"_.h5"),mat)
        writeDataH5(string(reldir,"/topological_charge_during_relaxation_T=",t,"_",
        "H=",h,"_DMI=",dmi,"_DDI=",ed,"_PMA=",pma,"_.h5"),filter(x->x!=0.,topChargeArray))
        writeDataH5(string(reldir,"/energy_during_relaxation_T=",t,"_","H=",h,
        "_DMI=",dmi,"_DDI=",ed,"_PMA=",pma,"_.h5"),filter(x->x!=0.,enArray))
        writeDataH5(string(reldir,"/effective_size_during_relaxation_T=",t,"_","H=",h,
        "_DMI=",dmi,"_DDI=",ed,"_PMA=",pma,"_.h5"),filter(x->x!=0.,reffArray))

    end

end
