"""
电磁感应理论模块 (Electromagnetic Induction Theory Module)

本模块实现了完整的电磁感应理论，包括：
1. 法拉第电磁感应定律 (Faraday's Law)
2. 楞次定律 (Lenz's Law)
3. 互感和自感 (Mutual and Self Inductance)
4. 电磁感应的应用 (Applications of Electromagnetic Induction)

作者：Julia Power Flow Team
日期：2025-09-10
"""

module ElectromagneticInduction

using LinearAlgebra
using Plots
using DifferentialEquations

export FaradayLaw, LenzLaw, SelfInductance, MutualInductance
export ElectromagneticField, InductionCoil, Transformer
export calculate_induced_emf, calculate_magnetic_flux, calculate_inductance
export simulate_electromagnetic_induction, plot_induction_analysis

"""
    ElectromagneticField

表示电磁场的结构体，包含电场和磁场信息
"""
struct ElectromagneticField
    magnetic_field::Vector{Float64}  # 磁场强度 [T]
    electric_field::Vector{Float64}  # 电场强度 [V/m]
    frequency::Float64               # 频率 [Hz]
    time::Vector{Float64}           # 时间序列 [s]
end

"""
    InductionCoil

感应线圈的结构体
"""
mutable struct InductionCoil
    turns::Int                      # 匝数
    area::Float64                   # 截面积 [m²]
    length::Float64                 # 长度 [m]
    permeability::Float64           # 磁导率 [H/m]
    resistance::Float64             # 电阻 [Ω]
    inductance::Float64             # 电感 [H]
    
    function InductionCoil(turns, area, length, permeability, resistance)
        inductance = calculate_self_inductance(turns, area, length, permeability)
        new(turns, area, length, permeability, resistance, inductance)
    end
end

"""
    Transformer

变压器结构体，用于演示互感应用
"""
mutable struct Transformer
    primary_coil::InductionCoil     # 初级线圈
    secondary_coil::InductionCoil   # 次级线圈
    mutual_inductance::Float64      # 互感 [H]
    coupling_coefficient::Float64   # 耦合系数
    core_permeability::Float64      # 铁芯磁导率
    
    function Transformer(primary, secondary, coupling_coeff, core_perm)
        mutual_ind = calculate_mutual_inductance(primary, secondary, coupling_coeff)
        new(primary, secondary, mutual_ind, coupling_coeff, core_perm)
    end
end

"""
理论推导部分 - 法拉第电磁感应定律

法拉第电磁感应定律的数学表达式：

1. 积分形式：
   ∮ E⃗ · dl⃗ = -dΦₘ/dt
   
   其中：
   - E⃗：电场强度
   - dl⃗：路径微元
   - Φₘ：磁通量
   - t：时间

2. 微分形式：
   ∇ × E⃗ = -∂B⃗/∂t
   
   其中：
   - B⃗：磁感应强度
   - ∇×：旋度算子

3. 感应电动势：
   ε = -N × dΦ/dt
   
   其中：
   - ε：感应电动势 [V]
   - N：线圈匝数
   - Φ：磁通量 [Wb]
"""

"""
    calculate_magnetic_flux(B, A, θ)

计算磁通量
Φ = B · A · cos(θ)

参数：
- B: 磁感应强度 [T]
- A: 面积 [m²]
- θ: 磁场与法线夹角 [rad]
"""
function calculate_magnetic_flux(B::Float64, A::Float64, θ::Float64=0.0)
    return B * A * cos(θ)
end

"""
    calculate_induced_emf(N, dΦ_dt)

根据法拉第定律计算感应电动势
ε = -N × dΦ/dt
"""
function calculate_induced_emf(N::Int, dΦ_dt::Float64)
    return -N * dΦ_dt
end

"""
    calculate_self_inductance(N, A, l, μ)

计算自感
L = μ × N² × A / l

参数：
- N: 匝数
- A: 截面积 [m²]
- l: 长度 [m]
- μ: 磁导率 [H/m]
"""
function calculate_self_inductance(N::Int, A::Float64, l::Float64, μ::Float64)
    return μ * N^2 * A / l
end

"""
    calculate_mutual_inductance(coil1, coil2, k)

计算两个线圈之间的互感
M = k × √(L₁ × L₂)

参数：
- coil1, coil2: 两个线圈
- k: 耦合系数 (0 ≤ k ≤ 1)
"""
function calculate_mutual_inductance(coil1::InductionCoil, coil2::InductionCoil, k::Float64)
    return k * sqrt(coil1.inductance * coil2.inductance)
end

"""
    faraday_law_demonstration(B_func, t_range, coil)

演示法拉第电磁感应定律
"""
function faraday_law_demonstration(B_func::Function, t_range::Vector{Float64}, coil::InductionCoil)
    println("=== 法拉第电磁感应定律演示 ===")
    println("理论基础：ε = -N × dΦ/dt")
    println("其中 Φ = B × A × cos(θ)")
    println()
    
    # 计算磁通量和感应电动势
    Φ = [calculate_magnetic_flux(B_func(t), coil.area) for t in t_range]
    
    # 数值计算磁通量变化率
    dt = t_range[2] - t_range[1]
    dΦ_dt = diff(Φ) / dt
    
    # 计算感应电动势
    ε = [calculate_induced_emf(coil.turns, dφ) for dφ in dΦ_dt]
    
    return t_range[1:end-1], Φ[1:end-1], ε
end

"""
    lenz_law_demonstration()

演示楞次定律 - 感应电流的方向总是阻碍磁通量的变化
"""
function lenz_law_demonstration()
    println("=== 楞次定律演示 ===")
    println("楞次定律：感应电流的方向总是阻碍磁通量的变化")
    println()
    println("数学表达：")
    println("当 dΦ/dt > 0 时，感应电流产生的磁场方向与原磁场相反")
    println("当 dΦ/dt < 0 时，感应电流产生的磁场方向与原磁场相同")
    println()
    
    # 示例：磁场增强和减弱的情况
    t = 0:0.01:2π
    B_increasing = t  # 磁场线性增强
    B_decreasing = 2π .- t  # 磁场线性减弱
    
    println("情况1：磁场线性增强 (B = t)")
    println("dB/dt = 1 > 0，感应电流产生反向磁场阻碍增强")
    println()
    
    println("情况2：磁场线性减弱 (B = 2π - t)")
    println("dB/dt = -1 < 0，感应电流产生同向磁场阻碍减弱")
    println()
    
    return t, B_increasing, B_decreasing
end

"""
    transformer_analysis(transformer, V_primary, f)

变压器分析 - 互感应用实例
"""
function transformer_analysis(transformer::Transformer, V_primary::Float64, f::Float64)
    println("=== 变压器分析 (互感应用) ===")
    println("理论基础：")
    println("初级：V₁ = -L₁(dI₁/dt) - M(dI₂/dt)")
    println("次级：V₂ = -M(dI₁/dt) - L₂(dI₂/dt)")
    println("变压比：V₂/V₁ = N₂/N₁")
    println()
    
    # 计算变压比
    turns_ratio = transformer.secondary_coil.turns / transformer.primary_coil.turns
    V_secondary = V_primary * turns_ratio
    
    # 计算电流比（忽略损耗）
    I_ratio = transformer.primary_coil.turns / transformer.secondary_coil.turns
    
    println("参数计算：")
    println("初级匝数 N₁ = $(transformer.primary_coil.turns)")
    println("次级匝数 N₂ = $(transformer.secondary_coil.turns)")
    println("变压比 = $(round(turns_ratio, digits=3))")
    println("初级电压 V₁ = $(V_primary) V")
    println("次级电压 V₂ = $(round(V_secondary, digits=2)) V")
    println("电流比 I₁/I₂ = $(round(I_ratio, digits=3))")
    println("互感 M = $(round(transformer.mutual_inductance * 1000, digits=2)) mH")
    
    return V_secondary, turns_ratio, I_ratio
end

"""
    rl_circuit_transient(L, R, V, t_range)

RL电路暂态分析 - 自感应用
"""
function rl_circuit_transient(L::Float64, R::Float64, V::Float64, t_range::Vector{Float64})
    println("=== RL电路暂态分析 (自感应用) ===")
    println("微分方程：L(di/dt) + Ri = V")
    println("解析解：i(t) = (V/R)(1 - e^(-Rt/L))")
    println("时间常数：τ = L/R")
    println()
    
    # 时间常数
    τ = L / R
    I_steady = V / R
    
    # 电流响应
    i_t = I_steady * (1 .- exp.(-t_range / τ))
    
    # 感应电压
    v_L = V * exp.(-t_range / τ)
    
    println("电路参数：")
    println("电感 L = $(L * 1000) mH")
    println("电阻 R = $(R) Ω")
    println("电源电压 V = $(V) V")
    println("时间常数 τ = $(round(τ * 1000, digits=2)) ms")
    println("稳态电流 I∞ = $(round(I_steady, digits=3)) A")
    
    return i_t, v_L, τ
end

"""
    electromagnetic_energy(L, I)

计算电磁场储能
W = (1/2) × L × I²
"""
function electromagnetic_energy(L::Float64, I::Float64)
    return 0.5 * L * I^2
end

"""
    simulate_electromagnetic_induction()

完整的电磁感应仿真示例
"""
function simulate_electromagnetic_induction()
    println("╔══════════════════════════════════════════════╗")
    println("║          电磁感应理论完整仿真              ║")
    println("║      Electromagnetic Induction Simulation   ║")
    println("╚══════════════════════════════════════════════╝")
    println()
    
    # 1. 创建感应线圈
    println("1. 创建感应线圈")
    coil = InductionCoil(
        100,        # 100匝
        0.01,       # 0.01 m² 截面积
        0.2,        # 0.2 m 长度
        4π*1e-7,    # 真空磁导率
        1.0         # 1 Ω 电阻
    )
    println("线圈参数：匝数=$(coil.turns), 面积=$(coil.area)m², 自感=$(round(coil.inductance*1e6, digits=2))μH")
    println()
    
    # 2. 法拉第定律演示
    println("2. 法拉第电磁感应定律演示")
    t_range = 0:0.01:2π
    B_func(t) = 0.1 * sin(2π * t)  # 正弦变化的磁场
    t_emf, Φ, ε = faraday_law_demonstration(B_func, collect(t_range), coil)
    println("最大感应电动势：$(round(maximum(abs.(ε)), digits=4)) V")
    println()
    
    # 3. 楞次定律演示
    println("3. 楞次定律演示")
    t_lenz, B_inc, B_dec = lenz_law_demonstration()
    
    # 4. 变压器分析
    println("4. 变压器分析")
    primary = InductionCoil(1000, 0.005, 0.1, 1000*4π*1e-7, 0.5)
    secondary = InductionCoil(500, 0.005, 0.1, 1000*4π*1e-7, 0.3)
    transformer = Transformer(primary, secondary, 0.98, 1000*4π*1e-7)
    
    V_sec, turns_ratio, I_ratio = transformer_analysis(transformer, 220.0, 50.0)
    println()
    
    # 5. RL电路暂态分析
    println("5. RL电路暂态分析")
    t_rl = 0:0.0001:0.01
    i_response, v_L, τ = rl_circuit_transient(0.01, 10.0, 12.0, collect(t_rl))
    println()
    
    # 6. 电磁能量计算
    println("6. 电磁场储能分析")
    max_current = maximum(i_response)
    max_energy = electromagnetic_energy(0.01, max_current)
    println("最大储能：$(round(max_energy * 1000, digits=4)) mJ")
    println()
    
    return (
        faraday_data = (t_emf, Φ, ε),
        lenz_data = (t_lenz, B_inc, B_dec),
        transformer_data = (V_sec, turns_ratio, I_ratio),
        rl_data = (collect(t_rl), i_response, v_L, τ),
        energy_data = max_energy
    )
end

"""
    plot_induction_analysis(simulation_data)

绘制电磁感应分析图表
"""
function plot_induction_analysis(sim_data)
    # 解包仿真数据
    t_emf, Φ, ε = sim_data.faraday_data
    t_lenz, B_inc, B_dec = sim_data.lenz_data
    t_rl, i_response, v_L, τ = sim_data.rl_data
    
    # 创建多子图
    p1 = plot(t_emf, Φ, title="Magnetic Flux Change", xlabel="Time (s)", ylabel="Magnetic Flux (Wb)", 
              linewidth=2, color=:blue, label="Φ(t)")
    
    p2 = plot(t_emf, ε, title="Induced EMF", xlabel="Time (s)", ylabel="EMF (V)", 
              linewidth=2, color=:red, label="ε(t)")
    
    p3 = plot(t_lenz, [B_inc B_dec], title="Lenz's Law Demo", xlabel="Time (s)", ylabel="Magnetic Field (T)",
              linewidth=2, label=["Increasing B" "Decreasing B"])
    
    p4 = plot(t_rl * 1000, [i_response v_L], title="RL Circuit Response", xlabel="Time (ms)", 
              ylabel="Current(A)/Voltage(V)", linewidth=2, label=["Current" "Inductor Voltage"])
    
    # 组合图形
    plot(p1, p2, p3, p4, layout=(2,2), size=(800, 600))
end

# 数学推导说明函数
"""
    print_mathematical_derivations()

打印完整的数学推导过程
"""
function print_mathematical_derivations()
    println("╔══════════════════════════════════════════════╗")
    println("║              数学推导详解                    ║")
    println("║        Mathematical Derivations             ║")
    println("╚══════════════════════════════════════════════╝")
    println()
    
    println("1. 法拉第电磁感应定律推导")
    println("═══════════════════════════")
    println("Maxwell方程组中的法拉第定律：")
    println("∮ E⃗ · dl⃗ = -∫∫ (∂B⃗/∂t) · dA⃗")
    println()
    println("对于单匝线圈：")
    println("ε = ∮ E⃗ · dl⃗ = -d/dt ∫∫ B⃗ · dA⃗ = -dΦ/dt")
    println()
    println("对于N匝线圈：")
    println("ε = -N × dΦ/dt")
    println()
    
    println("2. 自感推导")
    println("═══════════")
    println("长直螺线管的磁场：")
    println("B = μ × n × I = μ × (N/l) × I")
    println()
    println("磁通量：")
    println("Φ = B × A = μ × (N/l) × I × A")
    println()
    println("总磁通量链数：")
    println("Λ = N × Φ = μ × (N²/l) × A × I")
    println()
    println("自感定义：L = Λ/I")
    println("因此：L = μ × N² × A / l")
    println()
    
    println("3. 互感推导")
    println("═══════════")
    println("两个线圈的磁通量耦合：")
    println("Φ₁₂ = M × I₂")
    println("Φ₂₁ = M × I₁")
    println()
    println("根据能量守恒和磁场叠加原理：")
    println("M = k × √(L₁ × L₂)")
    println("其中 k 为耦合系数 (0 ≤ k ≤ 1)")
    println()
    
    println("4. RL电路微分方程推导")
    println("═══════════════════")
    println("基尔霍夫电压定律：")
    println("V = V_R + V_L = R×i + L×(di/dt)")
    println()
    println("微分方程：L×(di/dt) + R×i = V")
    println()
    println("齐次方程解：i_h = C×e^(-R×t/L)")
    println("特解：i_p = V/R")
    println()
    println("通解：i(t) = C×e^(-R×t/L) + V/R")
    println("初始条件 i(0) = 0：")
    println("i(t) = (V/R) × (1 - e^(-R×t/L))")
    println()
    
    println("5. 电磁场能量推导")
    println("═══════════════")
    println("功率关系：P = V×i = L×(di/dt)×i")
    println()
    println("能量：W = ∫₀ᵗ P dt = ∫₀ᵗ L×i×(di/dt) dt")
    println("     = L × ∫₀ⁱ i di = (1/2)×L×i²")
    println()
end

end # module ElectromagneticInduction
