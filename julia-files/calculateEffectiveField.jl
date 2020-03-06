# This module contains functions that calculate the effective field at a certain
# point of the spin lattice, and the function that builds the entire effective 
# field matrix.
#
# Useful functions: effectivefield!(),getFullEffField!()
# Note: effectivefield!() does not compute ddi, since we need to compute the
# field in the whole lattice to calculate DDI
module calculateEffectiveField
    using ddiFunctions,FFTW
    export effectivefield!,getFullEffField!,ddifield

    # inputs: mat = spin matrix, effField = (3,1) output float array, nx,ny = 
    # the position of the spin component we're interested in 
    # outputs: nothing. Ressults writte to effField
    #
    # This function modifies the (3,1) array to equal the effective field 
    # at some point nx,ny. Prealocating this way improves speed. 
    function effectivefield!(mat::Array{Float64,3}, effField::Array{Float64,1},
		nx::Int64, ny::Int64, params)

        j, h, dmi, anis, ed, pbc = params

        # The following accounts for zeeman contribution
        effField[1] = 0.
        effField[2] = 0.
        effField[3] = h

        exchangefield!(mat, effField, nx, ny, j, pbc==1.0)

        if dmi != 0.0
            dmiblochfield!(mat, effField, nx, ny, dmi, pbc==1.0)
        end
        if anis != 0.0
            anisotropyfield!(mat, effField, nx, ny, anis)
        end

    end

    # This function modifies effField to add the zeeman contribution to 
    # the effective field. 
    #
    # inputs: mat = spin matrix, effField = (3,1) effective field array,
    # nx,ny = positions where you want to find effective field
    # outputs: nothing. Results are stored in effField input 
    function exchangefield!(mat::Array{Float64,3}, effField::Array{Float64,1},
        nx::Int, ny::Int, J::Float64, pbc::Bool)
        
        p,m,n = size(mat)

        if nx > 1 
            for k in 1:3 effField[k] = effField[k] + mat[k,nx-1,ny] end
        else
            if pbc
                for k in 1:3 effField[k] = effField[k] + mat[k,m,ny] end
            end
        end

        if nx < m
            for k in 1:3 effField[k] = effField[k] + mat[k,nx+1,ny] end
        else
            if pbc
                for k in 1:3 effField[k] = effField[k] + mat[k,1,ny] end
            end
        end

        if ny > 1
            for k in 1:3 effField[k] = effField[k] + mat[k,nx,ny-1] end
        else
            if pbc
                for k in 1:3 effField[k] = effField[k] + mat[k,nx,n] end
            end 
        end

        if ny < n
            for k in 1:3 effField[k] = effField[k] + mat[k,nx,ny+1] end
        else
            if pbc
                for k in 1:3 effField[k] = effField[k] + mat[k,nx,1] end
            end 
        end

        for i in 1:3 effField[i] = J*effField[i] end
        
        nothing

    end

    # This function modifies effField to add the dmi contribution to the 
    # effective field at some point nx, ny.
    #
    # inputs: mat = spin matrix, effField = (3,1) effective field, 
    # nx,ny = pos, dmi = dmi const, pbc = periodic boundary conditions
    # outputs: nothing. Results are stored in effField input
    function dmiblochfield!(mat::Array{Float64,3}, effField::Array{Float64,1}, 
        nx::Int, ny::Int, dmi::Float64, pbc::Bool)

        p,m,n = size(mat)

        if ny==1
            
            if pbc
                effField[1]-=dmi*mat[3,nx,n];
                effField[3]+=dmi*mat[1,nx,n];
            end

            effField[1]+=dmi*mat[3,nx,ny+1];
            effField[3]-=dmi*mat[1,nx,ny+1];

        elseif ny==n

            if pbc
                effField[1]+=dmi*mat[3,nx,1];
                effField[3]-=dmi*mat[1,nx,1];
            end

            effField[1]-=dmi*mat[3,nx,ny-1];
            effField[3]+=dmi*mat[1,nx,ny-1];

        else

            effField[1]+=dmi*(mat[3,nx,ny+1]-mat[3,nx,ny-1]);
            effField[3]+=dmi*(mat[1,nx,ny-1]-mat[1,nx,ny+1]);

        end

        if nx==1

            if pbc

                effField[2]+=dmi*mat[3,m,ny];
                effField[3]-=dmi*mat[2,m,ny];

            end

            effField[2]-=dmi*mat[3,nx+1,ny];
            effField[3]+=dmi*mat[2,nx+1,ny];

        elseif nx==m

            if pbc
            
                effField[2]-=dmi*mat[3,1,ny];
                effField[3]+=dmi*mat[2,1,ny];
            
            end
            
            effField[2]+=dmi*mat[3,nx-1,ny];
            effField[3]-=dmi*mat[2,nx-1,ny];
        
        else
        
            effField[2]+=dmi*(mat[3,nx-1,ny]-mat[3,nx+1,ny]);
            effField[3]+=dmi*(mat[2,nx+1,ny]-mat[2,nx-1,ny]);
        
        end

    end

    # inputs: mat = spin matrix, effField = effective field at some nx, ny,
    # nx & ny are coordinattes of in the spin matrix.
    # outputs: nothing. This function changes effField to the effective field.
    function anisotropyfield!(mat::Array{Float64,3},effField::Array{Float64,1},
        nx::Int, ny::Int, pma::Float64)

        effField[3] += pma*mat[3,nx,ny]

    end

  
    # Gets the effective field for the entire spin array in matrix form
    #
    # inputs: mat = (N,N,3) spin array, Heff = (N,N,3) effective field matrix
    # effField = (3,1) effective field at a particular point
    # outputs: nothing
    function getFullEffField!(mat::Array{Float64,3},
        effField::Array{Float64,1}, params)
        
        p,m,n = size(mat)
        j, h, dmi, anis, ed, pbc = params

        dipField = Array{Float64}(undef,p,m,n)
        Heff = Array{Float64}(undef,p,m,n)

        for j in 1:n
            for i in 1:m
                effectivefield!(mat,effField,i,j,params)
                for k in 1:p Heff[k,i,j] = effField[k] end
            end 
        end

        if ed != 0.0
            dipField = ddifield(mat,ed)
            return Heff + dipField
        end

        return Heff
    end
  
    function ddifield(mat::Array{Float64,3}, ed::Float64)

        p, m, n = size(mat)
        field = Array{Float64}(undef,p,m,n)

        field = ed*FHD(mat,vddMatrices(m,n))

        return field

    end

    function fieldFspace(mat::Array{Float64,3})

        p,m,n = size(mat)
        j, h, dmi, anis, ed, pbc = params

        dipField = Array{Float64}(undef,p,m,n)
        Heff = Array{Float64}(undef,p,m,n)

        for j in 1:n
            for i in 1:m
                effectivefield!(mat,effField,i,j,params)
                for k in 1:p Heff[k,i,j] = effField[k] end
            end 
        end

        # fourier xform of the field so far
        Heff = fft(Heff)

        if ed != 0.0
            dipField = ddifield(mat,ed)
            return Heff + dipField
        end

        return Heff

    end
    
end
