module nBodySolver

struct State
    positions::Matrix{Float64}  # 3×N array for N bodies in 3D space
    velocities::Matrix{Float64} # 3×N array for N bodies in 3D space
    masses::Vector{Float64}     # N-element array for masses
    time::Float64               # Current time
end

end # module
