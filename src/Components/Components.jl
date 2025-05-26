module Components

using LinearAlgebra
export Bus, Branch, Generator, Load
export BusType, GeneratorType

# Bus types
@enum BusType begin
    PQ_BUS = 1
    PV_BUS = 2
    SLACK_BUS = 3
end

# Generator types
@enum GeneratorType begin
    THERMAL
    HYDRO
    WIND
    SOLAR
    NUCLEAR
end

include("bus.jl")
include("branch.jl")
include("generator.jl")
include("load.jl")

end # module Components
