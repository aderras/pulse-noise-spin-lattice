module calculateMagnetization

export calcM, calcM2

	function calcM(mat)

		p,m,n = size(mat)
		mz = 0.

		for i in eachindex(view(mat,1,1:m,1:n))

				mz = mz + mat[3,i[1],i[2]]

		end

		return mz/(m*n)

	end
	
	function calcM2(mat)

		p,m,n = size(mat)
		mz = 0.

		for i in eachindex(view(mat,1,1:m,1:n))

				mz = mz + (mat[3,i[1],i[2]]+1)^2

		end

		return mz/(n*m)

	end


end 

