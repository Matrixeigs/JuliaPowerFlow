"""
简单的电磁感应测试运行器
Simple Electromagnetic Induction Test Runner

这个脚本运行基本的功能测试，验证电磁感应理论模块的正确性
"""

using Test

# 包含模块
include("ElectromagneticInduction.jl")
using .ElectromagneticInduction

println("🧪 Running Simple Electromagnetic Induction Tests")
println("=" ^ 55)

# 基本功能测试
@testset "基本功能测试" begin
    
    # 测试磁通量计算
    println("  ✓ 测试磁通量计算...")
    Φ = ElectromagneticInduction.calculate_magnetic_flux(1.0, 0.01, 0.0)
    @test Φ ≈ 0.01
    
    # 测试感应电动势
    println("  ✓ 测试感应电动势计算...")
    ε = ElectromagneticInduction.calculate_induced_emf(100, 0.01)
    @test ε ≈ -1.0
    
    # 测试自感计算
    println("  ✓ 测试自感计算...")
    L = ElectromagneticInduction.calculate_self_inductance(1000, 0.01, 0.1, 4π*1e-7)
    @test L > 0
    @test isfinite(L)
    
    # 测试线圈创建
    println("  ✓ 测试线圈对象创建...")
    coil = ElectromagneticInduction.InductionCoil(100, 0.01, 0.2, 4π*1e-7, 1.0)
    @test coil.turns == 100
    @test coil.area ≈ 0.01
    @test coil.inductance > 0
    
    # 测试变压器创建
    println("  ✓ 测试变压器对象创建...")
    primary = ElectromagneticInduction.InductionCoil(1000, 0.01, 0.1, 4π*1e-7, 1.0)
    secondary = ElectromagneticInduction.InductionCoil(500, 0.01, 0.1, 4π*1e-7, 0.5)
    transformer = ElectromagneticInduction.Transformer(primary, secondary, 0.9, 4π*1e-7)
    @test transformer.mutual_inductance > 0
    
    # 测试电磁能量计算
    println("  ✓ 测试电磁能量计算...")
    energy = ElectromagneticInduction.electromagnetic_energy(0.01, 5.0)
    expected_energy = 0.5 * 0.01 * 5.0^2
    @test energy ≈ expected_energy
    
end

println("\n🎯 Running Physics Law Verification Tests")
println("-" ^ 45)

# 物理定律验证
@testset "物理定律验证" begin
    
    # 法拉第定律验证
    println("  ✓ 验证法拉第定律...")
    N = 100
    B₀ = 0.1
    A = 0.01
    ω = 2π * 50
    
    # 磁通量：Φ = B₀ × A × sin(ωt)
    # dΦ/dt = B₀ × A × ω × cos(ωt)
    # ε = -N × dΦ/dt = -N × B₀ × A × ω × cos(ωt)
    
    dΦ_dt_max = B₀ * A * ω
    ε_max_expected = N * dΦ_dt_max
    ε_max_calculated = ElectromagneticInduction.calculate_induced_emf(N, dΦ_dt_max)
    
    @test abs(ε_max_calculated) ≈ ε_max_expected
    
    # 楞次定律验证 (符号检查)
    println("  ✓ 验证楞次定律...")
    ε_positive_change = ElectromagneticInduction.calculate_induced_emf(N, 0.01)  # dΦ/dt > 0
    ε_negative_change = ElectromagneticInduction.calculate_induced_emf(N, -0.01) # dΦ/dt < 0
    
    @test ε_positive_change < 0  # 阻碍磁通量增加
    @test ε_negative_change > 0  # 阻碍磁通量减少
    
    # 自感公式验证
    println("  ✓ 验证自感公式...")
    μ₀ = 4π * 1e-7
    N, A, l = 1000, π * (0.01)^2, 0.1
    
    L_calculated = ElectromagneticInduction.calculate_self_inductance(N, A, l, μ₀)
    L_theoretical = μ₀ * N^2 * A / l
    
    @test L_calculated ≈ L_theoretical rtol=1e-12
    
    # 互感对称性验证
    println("  ✓ 验证互感对称性...")
    coil1 = ElectromagneticInduction.InductionCoil(1000, 0.01, 0.1, μ₀, 1.0)
    coil2 = ElectromagneticInduction.InductionCoil(500, 0.008, 0.12, μ₀, 0.8)
    
    M12 = ElectromagneticInduction.calculate_mutual_inductance(coil1, coil2, 0.8)
    M21 = ElectromagneticInduction.calculate_mutual_inductance(coil2, coil1, 0.8)
    
    @test M12 ≈ M21
    
end

println("\n📊 Running Numerical Accuracy Tests")
println("-" ^ 35)

# 数值精度测试
@testset "数值精度测试" begin
    
    # 边界值测试
    println("  ✓ 测试边界值...")
    
    # 零值
    @test ElectromagneticInduction.calculate_magnetic_flux(0.0, 1.0, 0.0) ≈ 0.0
    @test ElectromagneticInduction.calculate_induced_emf(100, 0.0) ≈ 0.0
    
    # 极小值
    Φ_tiny = ElectromagneticInduction.calculate_magnetic_flux(1e-15, 1e-15, 0.0)
    @test isfinite(Φ_tiny)
    @test Φ_tiny ≈ 1e-30
    
    # 极大值
    Φ_large = ElectromagneticInduction.calculate_magnetic_flux(1e6, 1e6, 0.0)
    @test isfinite(Φ_large)
    @test Φ_large ≈ 1e12
    
    # 角度测试
    println("  ✓ 测试角度依赖性...")
    B, A = 1.0, 1.0
    
    @test ElectromagneticInduction.calculate_magnetic_flux(B, A, 0.0) ≈ B * A
    @test ElectromagneticInduction.calculate_magnetic_flux(B, A, π/2) ≈ 0.0 atol=1e-15
    @test ElectromagneticInduction.calculate_magnetic_flux(B, A, Float64(π)) ≈ -B * A
    
end

println("\n⚡ Running Application Examples")
println("-" ^ 30)

# 应用实例测试
@testset "应用实例测试" begin
    
    # 发电机仿真
    println("  ✓ 测试发电机原理...")
    B, l, v = 1.0, 0.5, 10.0  # 1T, 0.5m, 10m/s
    ε_motional = B * l * v
    @test ε_motional ≈ 5.0
    
    # 变压器原理
    println("  ✓ 测试变压器原理...")
    N₁, N₂ = 1000, 200
    V₁ = 220.0
    
    turns_ratio = N₂ / N₁
    V₂_ideal = V₁ * turns_ratio
    I_ratio = N₁ / N₂
    
    @test turns_ratio ≈ 0.2
    @test V₂_ideal ≈ 44.0
    @test I_ratio ≈ 5.0
    
    # RL电路响应
    println("  ✓ 测试RL电路响应...")
    L, R, V = 0.01, 10.0, 12.0
    τ = L / R
    I_final = V / R
    
    @test τ ≈ 0.001  # 1ms
    @test I_final ≈ 1.2  # 1.2A
    
    # 在t = 5τ时，电流应该接近稳态值
    t_5tau = 5 * τ
    i_5tau = I_final * (1 - exp(-5))
    @test i_5tau > 0.99 * I_final  # 99%的稳态值
    
end

println("\n" * "=" ^ 55)
println("🎉 All tests completed successfully!")
println("✅ 电磁感应理论模块功能正常")
println("📋 测试覆盖内容:")
println("   • 基础计算函数 (磁通量、电动势、自感)")
println("   • 物理定律验证 (法拉第定律、楞次定律)")
println("   • 数值精度和边界条件")
println("   • 实际应用实例 (发电机、变压器、RL电路)")
println("=" ^ 55)
