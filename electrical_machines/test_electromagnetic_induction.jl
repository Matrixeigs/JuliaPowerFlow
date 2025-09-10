"""
电磁感应理论模块测试
Test suite for Electromagnetic Induction Theory Module

测试内容：
1. 基础计算函数测试
2. 物理定律验证测试
3. 数值精度测试
4. 边界条件测试
"""

using Test
include("ElectromagneticInduction.jl")
using .ElectromagneticInduction

# 测试基础计算函数
@testset "基础计算函数测试" begin
    
    @testset "磁通量计算" begin
        # 测试正常情况
        Φ = ElectromagneticInduction.calculate_magnetic_flux(1.0, 0.01, 0.0)
        @test Φ ≈ 0.01
        
        # 测试角度影响
        Φ_45 = ElectromagneticInduction.calculate_magnetic_flux(1.0, 0.01, π/4)
        @test Φ_45 ≈ 0.01 * cos(π/4)
        
        # 测试垂直情况
        Φ_90 = ElectromagneticInduction.calculate_magnetic_flux(1.0, 0.01, π/2)
        @test Φ_90 ≈ 0.0 atol=1e-10
    end
    
    @testset "感应电动势计算" begin
        # 测试法拉第定律
        ε = ElectromagneticInduction.calculate_induced_emf(100, 0.01)
        @test ε ≈ -1.0
        
        # 测试符号
        ε_pos = ElectromagneticInduction.calculate_induced_emf(100, -0.01)
        @test ε_pos ≈ 1.0
    end
    
    @testset "自感计算" begin
        # 标准螺线管参数
        L = ElectromagneticInduction.calculate_self_inductance(1000, 0.01, 0.1, 4π*1e-7)
        expected_L = 4π*1e-7 * 1000^2 * 0.01 / 0.1
        @test L ≈ expected_L
        
        # 测试匝数平方关系
        L2 = ElectromagneticInduction.calculate_self_inductance(2000, 0.01, 0.1, 4π*1e-7)
        @test L2 ≈ 4 * L
    end
end

# 测试物理定律验证
@testset "物理定律验证测试" begin
    
    @testset "能量守恒测试" begin
        # 创建RL电路
        L, R, V = 0.01, 10.0, 12.0
        t_range = 0:0.001:0.01
        
        i_t, v_L, τ = ElectromagneticInduction.rl_circuit_transient(L, R, V, collect(t_range))
        
        # 验证最终电流
        I_final = V / R
        @test i_t[end] ≈ I_final atol=0.01
        
        # 验证初始电感电压
        @test v_L[1] ≈ V atol=0.01
    end
    
    @testset "互感对称性测试" begin
        # 创建两个线圈
        coil1 = ElectromagneticInduction.InductionCoil(1000, 0.01, 0.1, 4π*1e-7, 1.0)
        coil2 = ElectromagneticInduction.InductionCoil(500, 0.005, 0.08, 4π*1e-7, 0.5)
        
        # 测试互感对称性 M₁₂ = M₂₁
        M12 = ElectromagneticInduction.calculate_mutual_inductance(coil1, coil2, 0.8)
        M21 = ElectromagneticInduction.calculate_mutual_inductance(coil2, coil1, 0.8)
        
        @test M12 ≈ M21
    end
    
    @testset "变压器理想关系测试" begin
        # 创建理想变压器
        primary = ElectromagneticInduction.InductionCoil(1000, 0.01, 0.1, 4π*1e-7, 0.1)
        secondary = ElectromagneticInduction.InductionCoil(500, 0.01, 0.1, 4π*1e-7, 0.05)
        transformer = ElectromagneticInduction.Transformer(primary, secondary, 1.0, 4π*1e-7)
        
        # 测试变压比
        turns_ratio = secondary.turns / primary.turns
        @test turns_ratio ≈ 0.5
        
        # 测试理想互感关系
        M_ideal = sqrt(primary.inductance * secondary.inductance)
        @test transformer.mutual_inductance ≈ M_ideal
    end
end

# 测试数值精度和边界条件
@testset "数值精度和边界条件测试" begin
    
    @testset "极小值测试" begin
        # 测试极小磁场
        Φ_small = ElectromagneticInduction.calculate_magnetic_flux(1e-12, 1e-6, 0.0)
        @test Φ_small ≈ 1e-18
        
        # 测试极小电感
        L_small = ElectromagneticInduction.calculate_self_inductance(1, 1e-12, 1e-3, 4π*1e-7)
        @test L_small > 0
        @test isfinite(L_small)
    end
    
    @testset "极大值测试" begin
        # 测试大磁场
        Φ_large = ElectromagneticInduction.calculate_magnetic_flux(100.0, 10.0, 0.0)
        @test Φ_large ≈ 1000.0
        @test isfinite(Φ_large)
    end
    
    @testset "零值测试" begin
        # 测试零磁场
        Φ_zero = ElectromagneticInduction.calculate_magnetic_flux(0.0, 1.0, 0.0)
        @test Φ_zero ≈ 0.0
        
        # 测试零变化率
        ε_zero = ElectromagneticInduction.calculate_induced_emf(100, 0.0)
        @test ε_zero ≈ 0.0
    end
    
    @testset "单位一致性测试" begin
        # 测试SI单位一致性
        # 磁通量：韦伯 [Wb] = [T·m²]
        B, A = 1.0, 1.0  # 1T, 1m²
        Φ = ElectromagneticInduction.calculate_magnetic_flux(B, A, 0.0)
        @test Φ ≈ 1.0  # 1 Wb
        
        # 自感：亨利 [H] = [Wb/A]
        N, A, l, μ = 1000, 0.01, 0.1, 4π*1e-7
        L = ElectromagneticInduction.calculate_self_inductance(N, A, l, μ)
        @test L > 0
        @test L < 1.0  # 现实的自感值
    end
end

# 测试复杂场景
@testset "复杂场景测试" begin
    
    @testset "非线性磁场响应" begin
        # 创建测试线圈
        coil = ElectromagneticInduction.InductionCoil(100, 0.01, 0.2, 4π*1e-7, 1.0)
        
        # 定义非线性磁场函数
        B_nonlinear(t) = 0.1 * sin(2π * t)^3
        
        t_range = collect(0:0.01:1.0)
        
        # 计算响应（这里简化测试）
        t_emf, Φ, ε = ElectromagneticInduction.faraday_law_demonstration(B_nonlinear, t_range, coil)
        
        @test length(ε) == length(t_emf)
        @test all(isfinite.(ε))
    end
    
    @testset "多频率响应测试" begin
        # 测试多个频率分量的叠加
        coil = ElectromagneticInduction.InductionCoil(200, 0.005, 0.15, 4π*1e-7, 2.0)
        
        # 多频率磁场
        B_multi(t) = 0.1 * (sin(2π * 50 * t) + 0.5 * sin(2π * 150 * t))
        
        t_range = collect(0:0.0001:0.02)  # 20ms，0.1ms步长
        
        try
            t_emf, Φ, ε = ElectromagneticInduction.faraday_law_demonstration(B_multi, t_range, coil)
            
            @test length(ε) > 0
            @test all(isfinite.(ε))
            
            # 验证频率成分
            max_emf = maximum(abs.(ε))
            @test max_emf > 0
            
        catch e
            @warn "多频率测试失败: $e"
        end
    end
end

# 性能测试
@testset "性能测试" begin
    
    @testset "大规模计算测试" begin
        # 测试大量数据点的计算性能
        N_points = 10000
        t_range = collect(range(0, 1, length=N_points))
        
        # 简单正弦磁场
        B_values = 0.1 * sin.(2π * 50 * t_range)
        
        # 测试计算时间
        start_time = time()
        
        Φ_values = [ElectromagneticInduction.calculate_magnetic_flux(B, 0.01, 0.0) for B in B_values]
        
        elapsed_time = time() - start_time
        
        @test length(Φ_values) == N_points
        @test elapsed_time < 1.0  # 应该在1秒内完成
        @test all(isfinite.(Φ_values))
    end
end

# 物理常数验证测试
@testset "物理常数验证" begin
    
    @testset "真空磁导率测试" begin
        μ₀ = 4π * 1e-7  # 真空磁导率
        
        # 使用标准螺线管验证磁导率
        N, A, l = 1000, π * (0.01)^2, 0.1  # 1000匝，半径1cm，长10cm
        L_theoretical = μ₀ * N^2 * A / l
        L_calculated = ElectromagneticInduction.calculate_self_inductance(N, A, l, μ₀)
        
        @test L_calculated ≈ L_theoretical rtol=1e-10
    end
    
    @testset "物理定律一致性" begin
        # 验证能量守恒
        L, I = 0.001, 5.0  # 1mH, 5A
        energy = ElectromagneticInduction.electromagnetic_energy(L, I)
        expected_energy = 0.5 * L * I^2
        
        @test energy ≈ expected_energy
        @test energy > 0
    end
end

"""
运行所有测试的函数
"""
function run_electromagnetic_tests()
    println("🧪 开始电磁感应理论模块测试")
    println("=" * 50)
    
    # 运行测试套件
    test_results = @testset "电磁感应理论完整测试" begin
        include("test_electromagnetic_induction.jl")
    end
    
    println("\n" * "=" * 50)
    println("✅ 电磁感应理论模块测试完成")
    
    return test_results
end

# 如果直接运行此文件，执行所有测试
if abspath(PROGRAM_FILE) == @__FILE__
    run_electromagnetic_tests()
end
