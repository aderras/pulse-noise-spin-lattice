
# This file contains several functions used to determine the dipole-
# dipole field in a spin lattice.
module ddiFunctions

    using main,FFTW,PaddedViews
    export FHD,convfft,buildVD,convolve_explicit,convolve_transpose,
    vddMatrices

    function buildVD(stype::String, Nx, Ny)

        if (Nx % 2 == Ny % 2 == 0)
            VD = zeros(2*Nx,2*Ny)
        else
            VD = zeros(2*Nx-1,2*Ny-1)
        end

        for i in -Nx+1:Nx-1
            for j in -Ny+1:Ny-1
                VD[j+Ny,i+Nx] = phi(i,j,stype)
            end
        end

        return VD
    end

    # Returns the xth and yth matrix element of phi
    function phi(nx::Int,ny::Int,stype::String)

        dim = 10
        dimFloat = float(dim)
        sum = 0.

        jx = jy = Float64

        nxfloat = nyfloat = Float64
        nxfloat = float(nx)
        nyfloat = float(ny)
        
        #implementing the so called smart formulas
        if stype == "xy"
            sum += ddixy(nxfloat,nyfloat) + 2/dim*sumddixy(nxfloat,nyfloat,dim)
        elseif stype == "xx" 
            sum += ddiSmart(nyfloat,nxfloat) + (2/dim)*sumddi(nyfloat,nxfloat,dim) 
        elseif stype == "yy" 
            sum += ddiSmart(nxfloat,nyfloat) + 2/dim*sumddi(nxfloat,nyfloat,dim)
        elseif stype == "zz" 
            sum += ddizz(nxfloat,nyfloat) + 2/dim*sumddizz(nxfloat,nyfloat,dim)
        end
        
        return sum
    end

    function sumddi(a,b,N)
        tot = 0.
        for i in 1:N-1
            tot += (N-i)*(2. *a^2-b^2-i^2)/denominator(a,b,float(i)) 
        end
        return tot
    end

    function sumddizz(a,b,N)
        tot = 0.
        for i in 1:N-1
            tot += (N-i)*(2. *i^2-b^2-a^2)/denominator(a,b,float(i)) 
        end
        return tot
    end

    function sumddixy(a,b,N)
        tot = 0.
        for i in 1:N-1
            tot += (N-i)*(3. *a*b)/denominator(a,b,float(i)) 
        end
        return tot
    end

    @inline ddiSmart(a,b) = if a==b==0 0. else (2. * a^2-b^2)/denominator(a,b,0.) end
    @inline ddizz(a,b) = if a==b==0 0. else (-1)/(sqrt((a^2 + b^2)^2 * (a^2 + b^2))) end
    @inline ddixy(a,b) = if a==b==0 0. else (3. * a * b)/denominator(a,b,0.) end
    # I have to write the denominator in this awful way because ^n for n>2 is slow in julia
    @inline denominator(a::Float64,b::Float64,c::Float64) =  sqrt((a^2 + b^2 + c^2)^2*(a^2 + b^2 + c^2)^2*(a^2 + b^2 + c^2))

    function vddMatrices(Nx::Int,Ny::Int) 
    
        phixx = Array{Float32}(undef,2*Nx, 2*Ny)
        phiyy = Array{Float32}(undef,2*Nx, 2*Ny)
        phizz = Array{Float32}(undef,2*Nx, 2*Ny)
        phixy = Array{Float32}(undef,2*Nx, 2*Ny)

        phixx = buildVD("xx",Nx,Ny)
        phiyy = buildVD("yy",Nx,Ny)
        phizz = buildVD("zz",Nx,Ny)
        phixy = buildVD("xy",Nx,Ny)
                
        return [phixx,phiyy,phizz,phixy]

    end

    function FHD(mat,phi)

        p,m,n = size(mat)

        matx = Array{Float64}(undef,m,n)
        maty = Array{Float64}(undef,m,n)
        matz = Array{Float64}(undef,m,n)
        matx,maty,matz = slicematrix(mat)
     
        return permutedims(cat(dims=3, convfft(matx,phi[1])+convfft(maty,phi[4]),
        convfft(maty,phi[2])+convfft(matx,phi[4]),convfft(matz,phi[3])),[3,1,2])

    end

    function FHDfourier(matF,phiF)

        p,m,n = size(matF)
        mm,nn = size(phiF)

        matFx = Array{Float64}(undef,mm,nn)
        matFy = Array{Float64}(undef,mm,nn)
        matFz = Array{Float64}(undef,mm,nn)
        matFx,matFy,matFz = slicematrixF(matF,mm,nn)
     
        return permutedims(cat(dims=3, matFx.*phiF[1] + matFy.*phiF[4],
        matFy.*phiF[2] + matFx.*phiF[4], matFz.*phiF[3]),[3,1,2])

    end

    function ddiconv(matslice,dimstr::String)

        m,n = size(matslice)
        conv = Array{Float64}(undef,m,n)
             
        conv = convfft(matslice,phi)

        return conv

    end
    
    # Convolution of two arrays using FFT
    function convfft(a,b)#,pfor,planifft)

        # FFTW.set_num_threads(2)

        ax,ay = size(a)
        bx,by = size(b)

        # pfor = plan_fft(Array{Float64,2}(undef,199,199));
        # pireal = plan_ifft!(Array{Float64,2}(undef,199,199));

        res     = Array{Float64}(undef,bx,by)

        apadded = PaddedView(0.0, a, (Base.OneTo(bx), Base.OneTo(by)))

        # res = ((planfft * apadded) .* (planfft * b))
        # planifft * res

        res = real(ifft(fft(apadded).*fft(b)))

        return res[bx-ax:bx-1,by-ay:by-1]

    end
    
    function slicematrix(A::Array{Float64,3})
        
        p,m,n = size(A)
        
        B1 = Array{Float64}(undef,m,n)
        B2 = Array{Float64}(undef,m,n)
        B3 = Array{Float64}(undef,m,n)
        
        for i in 1:m
            for j in 1:n
                B1[i,j] = A[1,i,j]
                B2[i,j] = A[2,i,j]
                B3[i,j] = A[3,i,j]
            end
        end
        
        return [B1,B2,B3]
    end

    function slicematrixF(A::Array{Float64,3},m,n)
                
        B1 = zeros(m,n)
        B2 = zeros(m,n)
        B3 = zeros(m,n)
        
        for i in 1:m
            for j in 1:n
                B1[i,j] = A[1,i,j]
                B2[i,j] = A[2,i,j]
                B3[i,j] = A[3,i,j]
            end
        end
        
        return [B1,B2,B3]
    end

    function convolve_transpose(K::Array{Float64,2},M::Array{Float64,2})

        (n,m) = size(K)     # get the size of the array
        (nn,mm) = size(M)

        res = zeros(Float64,nn-n+1,mm-m+1) 

        @inbounds @simd for ii in eachindex(view(res,1:nn-n+1,1:mm-n+1))
            res[ii] = convolve_elem_t(K,M,ii[1],ii[2])
        end

        return res[1:n,1:m]
    end

    function convolve_elem_t(K,M,xshift,yshift)

        tot = 0.0
        n = N()
        m = N()

        @simd for i in eachindex(view(K, 1:n, 1:m))
            @inbounds tot -= prod(K[i], M[i[1]+xshift-1, i[2]+yshift-1])
        end

        return tot

    end

    function convolve_elem(K,M,xshift,yshift)

        tot = 0.0
        (n,m) = size(K)

        for i in eachindex(view(K,1:n,1:m))
            tot += prod(K[i],M[ xshift + n - i[1], yshift + m - i[2]])
        end

        return tot

    end

    function convolve_explicit(K::Array{Float64,2},M::Array{Float64,2})

        (n,m) = size(K)     # get the size of the array
        (nn,mm) = size(M)

        res = zeros(Float64,nn-n+1,mm-m+1) # make a vector with the same number of rows as the matrix

        @inbounds for k=1:m # k from 1 to 100
            @inbounds @simd for l=1:n # l from 1 to 199

                res[l,k] = convolve_elem(K,M,l,k)
            end
        end
        return res[1:n,1:m]
    end


end