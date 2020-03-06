module LLequation

    using calculateEffectiveField,modifyFiles
    export RHS!,saveRHS

    # inputs: t = time of current step, mat = spin array, matRHS = right hand
    # side of the LL equation. This array is modified when calling the function
    # outputs: nothing
    function RHS!(t::Float64,mat::Array{Float64,3},matRHS::Array{Float64,3},
        lambda::Float64,matParams::Array{Float64,1})
       
        p, m, n = size(mat)

        Heff = Array{Float64}(undef,p,m,n)
        heff = Array{Float64}(undef,p)

        # calculate effective field 
        Heff = getFullEffField!(mat,heff,matParams)

        # # calculate S dot H. 
        SDotH = zeros(1,m,n)       
        SDotH .= sum(mat.*Heff,dims=1)
        
        fillRHS!(mat,Heff,SDotH,matRHS,lambda)

    end
    
    # Right side of LL equation
    #
    # inputs: mat = (N,N,3) spin array, Heff = (N,N,3) effective field matrix, 
    # SDotH = (N,N,1) array of SDotH, matRHS = (N,N,3) right hand side of LL equation. Updates this value
    # outputs: nothing
    function fillRHS!(mat::Array{Float64,3},Heff::Array{Float64,3},SDotH::Array{Float64,3},matRHS::Array{Float64,3},lambda::Float64)
    
        p, m, n = size(mat)

        for i in 1:m, j in 1:n
            matRHS[1,i,j] = mat[2,i,j]*Heff[3,i,j] - mat[3,i,j]*Heff[2,i,j] + lambda*(Heff[1,i,j]-mat[1,i,j]*SDotH[1,i,j])
            matRHS[2,i,j] = mat[3,i,j]*Heff[1,i,j] - mat[1,i,j]*Heff[3,i,j] + lambda*(Heff[2,i,j]-mat[2,i,j]*SDotH[1,i,j])
            matRHS[3,i,j] = mat[1,i,j]*Heff[2,i,j] - mat[2,i,j]*Heff[1,i,j] + lambda*(Heff[3,i,j]-mat[3,i,j]*SDotH[1,i,j])
        end
	
    end

end
