using Test
using JuliaPowerFlow

@testset "JuliaPowerFlow.jl" begin
    @testset "Component Tests" begin
        # Test Bus creation
        bus = Bus(1, bus_type=PQ_BUS, voltage_magnitude=1.05)
        @test bus.id == 1
        @test bus.voltage_magnitude == 1.05
        
        # Test Generator creation
        gen = Generator(1, 1, p_max=100.0, q_max=50.0)
        @test gen.id == 1
        @test gen.bus_id == 1
        @test gen.p_max == 100.0
        
        # Test capability constraints
        valid, msg = check_capability_constraints(gen, 50.0, 25.0)
        @test valid == true
    end
    
    @testset "Power System Tests" begin
        # Test power system creation
        sys = PowerSystem(100.0)
        @test sys.base_mva == 100.0
        @test length(sys.buses) == 0
        
        # Test adding components
        bus = Bus(1)
        add_component!(sys, bus)
        @test haskey(sys.buses, 1)
    end
    
    @testset "IEEE 9-Bus System" begin
        # Test IEEE 9-bus system creation
        sys = create_ieee9_system()
        @test get_bus_count(sys) == 9
        @test length(sys.generators) == 3
        @test length(sys.branches) == 9
        
        # Test power flow solution
        V, converged, iterations, S = newton_raphson_power_flow(sys, verbose=false)
        @test converged == true
        @test iterations < 10
        @test length(V) == 9
    end
end
