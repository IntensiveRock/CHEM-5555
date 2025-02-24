include("HardSphereMC.jl")
using .HardSphereMC
using Plots

#Problem 3.4: Calculate the Mean Squared Displacement to optimize the acceptance rate. Choose show_rate=true 
# as a quick and dirty way to judge if you're equilibrated before you start collecting data.
function compute_D(; T::Int64=100,target_accept::Float64=0.5)
    sweeps=50000
    equil_time = 10000
    calc_freq = 10

    hs = hard_spheres("parameters.toml")
    hs.acceptance = target_accept #Adjust the target acceptance rate in the function.
    c = msd(T=T)
    calcs = [c]

    #Run the simulation
    mc_simulation!(sys=hs,
        sweeps=sweeps,
        calculations=calcs,
        equil_time=equil_time,
        show_rate=false,
        calculate_freq=calc_freq)

    #Extract the data
    r² = calcs[1].r²
    t = collect(0:(T-1))
    
    mean(x) = sum(x)/length(x)
    std(x) = sqrt(sum((x .- mean(x)).^2)/length(x))
    Dt = calc_freq*r²[2:end]./t[2:end]/6
    D = mean(Dt)
    σ = std(Dt)/sqrt(length(Dt))#Use the standard error estimate.
    
    p = plot(t*calc_freq, r², label="MSD", xlabel="Time", ylabel="MSD")
    return p, D, σ, hs
end

#Now calculate the Diffusion constant as a function of the target acceptance rate.

rates = collect(0.2:0.05:0.55)
combined_plot = plot(xlabel="Time", ylabel="MSD", title="MSD vs Time for Different Acceptance Rates")


println("Calculating the Diffusion constant for different acceptance rates")
println("D = ⟨r²(t)/6t⟩")
println("Errors estimated from the standard error of the mean.")

# Initialize arrays outside the loop
D_values = Float64[]
σ_values = Float64[]

for rate in rates
    p, D, σ, sys = compute_D(T=50, target_accept=rate)
    plot!(combined_plot, p.series_list[1].plotattributes[:x], 
          p.series_list[1].plotattributes[:y], 
          label="acceptance rate = $(round(rate, digits=2))")
    println("D = $D ± $σ for rate = $rate")
    #Write the values to a file
    open("Diffusion Constant vs Acceptance for N=$(number_of_particles(sys)).dat", "a") do file
        if rate == rates[1]
            write(file, "rate D error (σ) \n")
        end
        write(file, "$rate, $D, $σ\n")
    end
    
    # Collect values
    push!(D_values, D)
    push!(σ_values, σ)
    
    # Update display after each calculation
    display(combined_plot)
end

# Create diffusion constant plot after collecting all data
p_diff = plot(rates, D_values, yerror=σ_values, 
              xlabel="Acceptance Rate", 
              ylabel="Diffusion Constant D",
              label="D vs. Acceptance Rate",
              marker=:circle)
println("\nPlot of the Diffusion constant vs acceptance rate.")
display(p_diff)

println("\nPlot of the MSD vs time for different acceptance rates, verifying linear diffusion.")
display(combined_plot)

savefig(p_diff, "Diffusion_vs_Acceptance_Rate.png")
savefig(combined_plot, "MSD_vs_Time.png")