# This module contains functions used to calculate the topological charge 
# of a spin field
#
#
module calculateTopologicalCharge

    export calcQ

    # The following function calculates the topological charge of
    # the entire lattice.
    #
    # inputs: mat = NxNx3 spin array
    # outputs: q = float
    function calcQ(mat::Array{Float64,3})
        
        p,m,n = size(mat)
        q = 0.0

        dxsQ = Array{Float64}(undef,p,m,n)
        dysQ = Array{Float64}(undef,p,m,n)
        # dxsQ = zeros(p,m,n)
        # dysQ = zeros(p,m,n)


        dxsQ = (-rotateRight!(mat,[2,0,0])+8*rotateRight!(mat,[1,0,0])-
               8*rotateRight!(mat,[-1,0,0])+rotateRight!(mat,[-2,-0,0]))/12
        dysQ = (-rotateRight!(mat,[0,2,0])+8*rotateRight!(mat,[0,1,0])-
                8*rotateRight!(mat,[0,-1,0])+rotateRight!(mat,[0,-2,0]))/12

        for i in 1:m, j in 1:n

            q += mat[1,i,j]*dxsQ[2,i,j]*dysQ[3,i,j] + 
            mat[2,i,j]*dxsQ[3,i,j]*dysQ[1,i,j] + 
            mat[3,i,j]*dxsQ[1,i,j]*dysQ[2,i,j] - 
            mat[3,i,j]*dxsQ[2,i,j]*dysQ[1,i,j] - 
            mat[2,i,j]*dxsQ[1,i,j]*dysQ[3,i,j] - 
            mat[1,i,j]*dxsQ[3,i,j]*dysQ[2,i,j]  
            
        end

        q = q/(4*pi)

        return q

    end

    # This function shifts the rows and columns of mat a certain amount 
    # defined in rotNum.
    #
    # inputs: mat = NxNx3 spin array, rotNum = 3x1 array of integers 
    # containing the amount of times to shift rows and columns.
    # outputs: matNew = rotated matix.
    function rotateRight!(mat::Array{Float64,3},rotNum::Array{Int64,1})

        p,m,n = size(mat)
        
        matNew = Array{Float64}(undef,p,m,n)
        matTemp = Array{Float64}(undef,p,m,n)

        matTemp .= mat

        for i in 1:m, j in 1:n
            shiftCols!(matTemp,matNew,i,j,rotNum[2])
        end

        matTemp .= matNew

        for i in 1:m, j in 1:n
            shiftRows!(matTemp,matNew,i,j,rotNum[1])
        end

        return matNew

    end

    # The following function shifts the element in a row in
    # a periodic way such that the last element becomes the first.
    # 
    # inputs: mat = NxNx3 spin matrix, matNew = NxNx3 shifted mat,
    # i,j = index of the element to shift, ishift = amount to shift by
    # outputs: nothing. The result is stored in matNew.
    function shiftRows!(mat::Array{Float64,3}, matNew::Array{Float64,3}, 
        i::Int64, j::Int64, ishift::Int64)

        p,m,n = size(mat)

        if abs(ishift) > m
            println("Error in shiftRows! function. Shift is larger than array.")
        end

        if ishift > 0
            if i > ishift
                for k in 1:p matNew[k,i,j] = mat[k,i-ishift,j] end 
            else
                for k in 1:p matNew[k,i,j] = mat[k,m-ishift+i, j] end
            end
        elseif ishift < 0
            if i <= m + ishift
                for k in 1:p matNew[k,i,j] = mat[k,i-ishift,j] end
            else
                for k in 1:p matNew[k,i,j] = mat[k,(i-ishift)-m,j] end
            end
        else
            for k in 1:p matNew[k,i,j] = mat[k,i,j] end
        end

    end

    # The following function shifts the element in an column in
    # a periodic way such that the last element becomes the first.
    # 
    # inputs: mat = NxNx3 spin matrix, matNew = NxNx3 shifted mat,
    # i,j = index of the element to shift, jshift = amount to shift by
    # outputs: nothing. The result is stored in matNew.
    function shiftCols!(mat::Array{Float64,3}, matNew::Array{Float64,3}, 
        i::Int64, j::Int64, jshift::Int64)

        p,m,n = size(mat)

        if jshift > 0
            if j > jshift
                for k in 1:p matNew[k,i,j] = mat[k,i,j-jshift] end 
            else
                for k in 1:p matNew[k,i,j] = mat[k,i, n-jshift+j] end
            end
        elseif jshift < 0
            if j <= n + jshift
                for k in 1:p matNew[k,i,j] = mat[k,i,j-jshift] end 
            else
                for k in 1:p matNew[k,i,j] = mat[k,i, (j-jshift)-n] end
            end
        else
            for k in 1:p matNew[k,i,j] = mat[k,i,j] end
        end

    end

end