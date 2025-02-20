module Sampling

export MC, find_Δx, metropolis

function MC(x0::Float64, βΦ::Function, nsteps::Int64, Δx::Float64; criterion::Function=metropolis)
    x = Vector{Float64}(undef, nsteps)#Trajectory
    x[1] = x0
    acceptance_rate = 0.0
    for t in 1:(nsteps-1)
        xₜ = x[t]
        y = xₜ + (rand() - 0.5)*Δx
        # Metropolis criterion
        ΔβΦ = βΦ(y) - βΦ(xₜ)
        accept = criterion(ΔβΦ)
        if accept
            x[t+1] = y
            acceptance_rate += 1.0
        else
            x[t+1] = xₜ
        end
    end
    acceptance_rate /= (nsteps - 1)
    return acceptance_rate, x
end

function metropolis(ΔβΦ::Float64)
    if ΔβΦ <= 0.0
        return true
    end
    if rand() <= exp(-ΔβΦ)
        return true
    end
    return false
end

function find_Δx(x0::Float64, βΦ::Function, nequil::Int64;  # Required arguments first
                criterion::Function=metropolis,               # Optional arguments after ;
                acceptance_target::Float64=0.3, 
                window_size::Int64=10, 
                tol::Float64=1e-3)
    #Calibrate the step size
    Δx = 0.001
    for _ in 1:nequil
        sliding_window = 0.0
        for _ in 1:window_size
            acceptance_rate, _ = MC(x0, βΦ, nequil, Δx, criterion=criterion)
            sliding_window += acceptance_rate
        end
        sliding_window /= window_size
        factor = sliding_window/acceptance_target
        Δx *= sqrt(factor)

        if abs(sliding_window - acceptance_target) < tol
            println("Converged to Δx = $Δx with acceptance rate = $sliding_window")
            return Δx  
        end
    end
    @warn "Step size did not converge"
    return Δx
end

end