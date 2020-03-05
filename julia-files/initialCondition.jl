module initialCondition

    using main
    export buildInitial

    @inline xcomp(z::Complex) = if isnan(real(z)/(1+(abs(z)^2)/4)) 0. else real(z)/(1+(abs(z)^2)/4) end
    @inline ycomp(z::Complex) = if isnan(imag(z)/(1+(abs(z)^2)/4)) 0. else imag(z)/(1+(abs(z)^2)/4) end
    @inline zcomp(z::Complex) = if isnan((-1+(abs(z)^2)/4)/(1+(abs(z)^2)/4)) 1. else -(-1+(abs(z)^2)/4)/(1+(abs(z)^2)/4) end
    
    @inline skyr(x,y,z0,l) = (x+y*1im-z0)/l
    @inline antiskyr(x,y,z0,l) = l/(x+y*1im-z0)
    @inline skyrAnti(x,y,z1,z2,l) = (x+y*1im-z1)*(x-y*1im-z2)/l^2

    function buildInitial(t::String,l::Float64,z0,Nx::Int64,Ny::Int64)

        # mat = Array{Float64}(undef,Nx,Ny,3)
        mat = Array{Float64}(undef,3,Nx,Ny)
        xx = 0.
        yy = 0.

        xstart = Nx/2
        ystart = Ny/2
        dx = 1
        dy = 1

        vec = Array{Float64}(undef,3)
        if t == "skyrmion"
            for i in eachindex(view(mat,1,1:Nx,1:Ny))

                skyrmioncomponent!(vec,i[1]-(Nx+1)/2,i[2]-(Ny+1)/2,l*2.)
            
                mat[1,i] = -vec[1]
                mat[2,i] = -vec[2]
                mat[3,i] = vec[3]


            end
        elseif t == "skyrmionAntiskyrmion"
            for i in eachindex(view(mat,1,1:Nx,1:Ny))

                xx = -xstart + dx*(i[1]) 
                yy = -ystart + dy*(i[2])

                z = skyrAnti(xx,yy,z0,-z0,l)
                mat[1,i] = xcomp(z)
                mat[2,i] = ycomp(z)
                mat[3,i] = zcomp(z)
            
            end
        elseif t == "antiskyrmion"
            for i in eachindex(view(mat,1,1:Nx,1:Ny))

                xx = -xstart + dx*(i[1] + 1) 
                yy = -ystart + dy*(i[2] + 1)

                z = antiskyr(xx,yy,z0,l)
                mat[1,i] = xcomp(z)
                mat[2,i] = ycomp(z)
                mat[3,i] = zcomp(z)
            
            end
        elseif t == "AFMskyrmion"
            for i in 1:Nx
                for j in 1:Ny
                    xx = -xstart + dx*(i + 1) 
                    yy = -ystart + dy*(j + 1)

                    z = skyr(xx,yy,z0,l)

                    if i%2==0 && j%2==1 ||
                        j%2==0 && i%2==1
                        mat[1,i,j] = -xcomp(z)
                        mat[2,i,j] = -ycomp(z)
                        mat[3,i,j] = -zcomp(z)
                    else
                        mat[1,i,j] = xcomp(z)
                        mat[2,i,j] = ycomp(z)
                        mat[3,i,j] = zcomp(z)
                    end
                end
            end
        elseif t == "AFMskyrmionAntiskyrmion"
            for i in 1:Nx
                for j in 1:Ny
                    xx = -xstart + dx*(i + 1) 
                    yy = -ystart + dy*(j + 1)

                    z = skyrAnti(xx,yy,z0,-z0,l)

                    if i%2==0 && j%2==1 ||
                        j%2==0 && i%2==1
                        mat[1,i,j] = -xcomp(z)
                        mat[2,i,j] = -ycomp(z)
                        mat[3,i,j] = -zcomp(z)
                    else
                        mat[1,i,j] = xcomp(z)
                        mat[2,i,j] = ycomp(z)
                        mat[3,i,j] = zcomp(z)
                    end
                end
            end
        end

        return mat

    end

    # inputs: x,y = position on the plane of the centrosymmetic skyrmion,
    # λ = the radius, Q = topological charge, ϕ = chirality
    # output: {x,y,z} vector representing skyrmion components at x,y
    function skyrmioncomponent!(val::Array{Float64},x::Float64,
        y::Float64,λ::Float64)

        Q = 1.
        ϕ = pi/2

        rsquared = sz = sperp = ang = 0.

        rsquared = x^2+y^2
        sz = -(rsquared-λ^2)/(rsquared+λ^2)
        sperp = sqrt(1-sz^2)

        if y>0
            ang = acos(x/sqrt(rsquared))
        else
            ang = 2*pi-acos(x/sqrt(rsquared))
        end

        # for some reason allocating this way leads to smaller allocation time
        val[1] = sperp*cos(Q*ang+ϕ)
        val[2] = sperp*sin(Q*ang+ϕ)
        val[3] = sz
        
        nothing    
    end


end
