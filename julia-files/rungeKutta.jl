# This module contains the Runge-Kutta stepping function.
module rungeKutta

    using calculateEffectiveField,LLequation
    export rk4!

    # Calling this function once leads to 'n' iterations of the 
    # fourth order runge-kutta solver. (See vars in main.jl for what
    # the values mean.) 
    # Rather than return an array of the past n solutions, this modifies
    # the X0 array to constantly contain most recent solution. 
    function rk4!(X0::Array{Float64,3},f::Function, matParams::Array{Float64,1},
         evalParams::Array{Float64,1})
        
        p, m, n = size(X0)

        Δt, nn, lambda, T = evalParams

        # Initialize arrays for the RK steps
        K1 = Array{Float64}(undef,p,m,n)
        K2 = Array{Float64}(undef,p,m,n)
        K3 = Array{Float64}(undef,p,m,n)
        K4 = Array{Float64}(undef,p,m,n)
        X = Array{Float64}(undef,p,m,n)

        temp = zeros(p,m,n)

        # # # Initialize vars
        h = t = Float64
        h = Δt/nn       # h = 1/5 = 0.2 
        t = 0.0

        X .= X0 # starting with specified

        for k in 1:nn
                   
            # In this particular instance the equation were solving
            # has no t in RHS, so the t value is arbitrary.
            # t = t0 + k*h

            f(t, X, K1, lambda, matParams) # function modifies last argument to equal RK Coeff
            K1 .= h*K1

            temp .= X .+ 0.5 .* K1
            f(t + 0.5*h, temp , K2, lambda, matParams)
            K2 .= h*K2

            temp .= X .+ 0.5 .* K2
            f(t + 0.5*h, temp, K3, lambda, matParams)
            K3 .= h*K3

            temp .= X .+ K3
            f(t + h, temp, K4, lambda, matParams)
            K4 .= h*K4

            X .= X .+ (1/6).*K1 .+ (1/3).*K2 .+ (1/3).*K3 .+ (1/6).*K4
        end

        X0 .= X # rewrite input variable with the result

    end

end
