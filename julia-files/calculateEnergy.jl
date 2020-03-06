#!/usr/bin/env julia
push!(LOAD_PATH, pwd())
# This module contains functions that calculate the energy of a 
# spin lattice. 
module calculateEnergy
    using ddiFunctions
    export calcEnergy

    # inputs: spin array (N,N,3), matParams = [H, DMI, PMA, ED, PBC]
    # outputs: Float64 of energy 
    function calcEnergy(mat::Array{Float64,3}, params)

        J, H, dmi, Anis, ed, pbc = params

        totalenergy = exchangeenergy(mat, J, pbc==1.0)

        if H != 0.0
            totalenergy += zeemanenergy(mat, H) 
        end
        if dmi != 0.0
            totalenergy += dmienergy(mat, dmi, pbc==1.0)
        end
        if Anis != 0.0
            totalenergy += anisotropyenergy(mat, Anis)
        end
        if ed != 0.0
            totalenergy += ddiEnergy(mat, ed)
        end
        
        return totalenergy
    end

    # inputs: mat = spin matrix
    # outputs: the zeeman contribution to the energy
    function zeemanenergy(mat::Array{Float64,3}, Hext::Float64)

        p, m, n = size(mat)

        en = 0.0
        
        for i in 1:n 
            for j in 1:m
                en += -Hext * mat[3,i,j]
            end
        end
        
        return en
    
    end 

    # inputs: mat = spin matrix
    # outputs: the dmi contribution of the energy
    function dmienergy(mat::Array{Float64,3}, DMI::Float64, pbc::Bool)

        p, m, n = size(mat)
        en = 0.0

        # The following is for nx > 1
        for i in 2:m
            for j in 1:n
                en += mat[2,i,j] * mat[3,i-1,j] - mat[3,i,j]*mat[2,i-1,j]
            end
        end
        
        # The following is for ny > 1
        for i in 1:m
            for j in 2:n
                en+=mat[3,i,j]*mat[1,i,j-1] - mat[1,i,j]*mat[3,i,j-1]
            end
        end
        
        # Deal with the edge of the matrix. If periodic boundary 
        # conditions, run the following. Otherwise do nothing
        if pbc
            for i in 1:n 
                en += (mat[2,1,i]*mat[3,m,i] - mat[3,1,i]*mat[2,m,i])
            end
            
            for i in 1:m 
                en += (mat[3,i,1]*mat[1,i,n] - mat[1,i,1]*mat[3,i,n])
            end
        end

        en = -DMI * en
        
        return en
    end

    # inputs: mat = spin Array
    # outputs: en = exchange energy of the spin array, mat
    function exchangeenergy(mat::Array{Float64,3}, J::Float64, pbc::Bool)

        p, m, n = size(mat)
        en = 0.0
        
        for j in 1:n       
            for i in 1:m-1
                for k in 1:p
                    en += mat[k,i,j] * mat[k,i+1,j]
                end
            end
        end

        for j in 1:n-1    
            for i in 1:m   
                for k in 1:p
                    en += mat[k,i,j] * mat[k,i,j+1]
                end
            end
        end

        if pbc
            for j in 1:n
                for k in 1:p
                    en += mat[k,1,j] * mat[k,m,j]
                end
            end
            for i in 1:m
                for k in 1:p
                    en += mat[k,i,1] * mat[k,i,n]
                end
            end

            en -= (2 * m * n)
        else 
            en -= (2 * m * n - m - n)
        end

        return -J * en
    end

    # inputs: mat = spin Array
    # outputs: en = exchange energy of the spin array, mat
    function anisotropyenergy(mat::Array{Float64,3}, AnisZ::Float64)
        
        en = 0.0
        p, m, n = size(mat)

        for i in 1:m 
            for j in 1:n
                en -= 0.5 * AnisZ * (mat[3,i,j] * mat[3,i,j] - 1)
            end
        end

        return en
    end

    function ddiEnergy(mat::Array{Float64,3}, ed::Float64)

        p, m, n = size(mat)
        ddiEnArray = Array{Float64}(undef,p,m,n)

        ddiEnArray = FHD(mat,vddMatrices(m,n))
        tot = 0.0

        for i in eachindex(view(mat,:,:,:))
            tot += mat[i] * 0.5 * ddiEnArray[i]
        end

        return -ed * tot

    end

end
