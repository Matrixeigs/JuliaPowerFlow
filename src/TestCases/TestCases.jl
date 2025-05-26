module TestCases

using ..Components, ..Systems

export create_ieee9_system, case9_legacy

include("ieee9.jl")
include("case9_legacy.jl")

end # module TestCases
