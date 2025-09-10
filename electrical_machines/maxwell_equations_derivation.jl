# 麦克斯韦方程组推导与数值验证 (修复版)
# Maxwell Equations Derivation and Numerical Verification (Fixed Version)

using LinearAlgebra, Plots, DifferentialEquations
import Base.MathConstants: π  # 明确导入π避免冲突

# 设置绘图后端以避免中文字体问题
gr()  # 使用GR后端，但用英文标签

"""
物理常数定义模块
Physical Constants Module
"""
module PhysicalConstants
    # 基本物理常数 / Fundamental Physical Constants
    const μ₀ = 4π * 10.0^(-7)        # 真空磁导率 [H/m] - Vacuum permeability
    const ε₀ = 8.854187817 * 10.0^(-12)  # 真空介电常数 [F/m] - Vacuum permittivity  
    const c = 299792458         # 光速 [m/s] - Speed of light
    const e_charge = 1.602176634 * 10.0^(-19)   # 基本电荷 [C] - Elementary charge
    const k_e = 1/(4π*ε₀)       # 库仑常数 [N⋅m²/C²] - Coulomb constant
    
    # 验证光速关系 / Verify speed of light relation
    const c_calculated = 1/sqrt(μ₀ * ε₀)
    
    function verify_constants()
        println("=== 物理常数验证 / Physical Constants Verification ===")
        println("真空磁导率 mu_0 = $μ₀ H/m")
        println("真空介电常数 epsilon_0 = $ε₀ F/m")
        println("光速定义值 c = $c m/s")
        println("光速计算值 c = 1/√(mu_0*epsilon_0) = $c_calculated m/s")
        println("相对误差 = $(abs(c - c_calculated)/c * 100)%")
        println("库仑常数 k_e = $(k_e) N⋅m²/C²")
        println()
    end
end

"""
库仑定律与高斯定律推导
Coulomb's Law and Gauss's Law Derivation
"""
function coulomb_to_gauss_derivation()
    println("=== 从库仑定律到高斯定律推导 ===")
    println("=== Coulomb's Law to Gauss's Law Derivation ===")
    
    # 库仑定律：F = k_e * q1*q2/r²
    println("1. 库仑定律 / Coulomb's Law:")
    println("   F = (1/4πε₀) × (q₁q₂/r²) × r̂")
    println("   F = k_e × (q₁q₂/r²) × r̂")
    
    # 电场定义
    println("\n2. 电场定义 / Electric Field Definition:")
    println("   E = F/q = k_e × (Q/r²) × r̂")
    
    # 数值示例：点电荷电场
    Q = 1.0 * 10.0^(-6)  # 1微库仑
    r_values = [0.1, 0.2, 0.5, 1.0, 2.0]  # 距离 [m]
    
    println("\n3. 数值验证 - 点电荷电场强度 / Numerical Verification:")
    println("   电荷 Q = $Q C")
    for r in r_values
        E_field = PhysicalConstants.k_e * Q / r^2
        println("   r = $(r) m → E = $(round(E_field, sigdigits=4)) V/m")
    end
    
    # 高斯定律推导
    println("\n4. 高斯定律推导 / Gauss's Law Derivation:")
    println("   通过半径r的球面，电通量 Phi_E = ∮ E · dA")
    println("   球面面积 = 4πr²")
    println("   Phi_E = E × 4πr² = k_e × Q/r² × 4πr² = Q/ε₀")
    println("   因此：∮ E · dA = Q/ε₀ (高斯定律)")
    
    # 微分形式
    println("\n5. 微分形式 / Differential Form:")
    println("   ∇ · E = ρ/ε₀")
    
    return r_values, [PhysicalConstants.k_e * Q / r^2 for r in r_values]
end

"""
安培定律与位移电流
Ampère's Law and Displacement Current
"""
function ampere_to_maxwell_derivation()
    println("\n=== 安培定律到安培-麦克斯韦定律推导 ===")
    println("=== Ampère's Law to Ampère-Maxwell Law Derivation ===")
    
    # 原始安培定律
    println("1. 原始安培定律 / Original Ampère's Law:")
    println("   ∮ B · dl = μ₀I")
    println("   微分形式：∇ × B = μ₀J")
    
    # 电荷守恒问题
    println("\n2. 电荷守恒问题 / Charge Conservation Problem:")
    println("   电荷守恒：∂ρ/∂t + ∇ · J = 0")
    println("   对安培定律取散度：∇ · (∇ × B) = μ₀∇ · J")
    println("   但 ∇ · (∇ × B) = 0 (恒等式)")
    println("   这要求 ∇ · J = 0，与电荷守恒矛盾！")
    
    # 麦克斯韦修正
    println("\n3. 麦克斯韦修正 / Maxwell's Correction:")
    println("   结合高斯定律：∇ · E = ρ/ε₀")
    println("   电荷守恒：∂ρ/∂t + ∇ · J = 0")
    println("   得到：∇ · J + ε₀∂(∇ · E)/∂t = 0")
    println("   即：∇ · (J + ε₀∂E/∂t) = 0")
    
    # 位移电流
    println("\n4. 位移电流 / Displacement Current:")
    println("   定义位移电流密度：J_D = ε₀∂E/∂t")
    println("   修正的安培定律：∇ × B = μ₀(J + J_D)")
    println("   即：∇ × B = μ₀J + μ₀ε₀∂E/∂t")
    
    # 数值示例：平行板电容器
    println("\n5. 数值示例 - 平行板电容器充电:")
    C = 1.0 * 10.0^(-6)  # 电容 [F]
    V = 100.0   # 电压 [V] 
    t_charge = 1.0 * 10.0^(-3)  # 充电时间 [s]
    
    I_conduction = C * V / t_charge  # 传导电流
    println("   电容 C = $C F")
    println("   充电电压 V = $V V")
    println("   充电时间 t = $t_charge s")
    println("   传导电流 I = CV/t = $(round(I_conduction, sigdigits=4)) A")
    
    # 位移电流计算
    A = 0.01  # 电容器板面积 [m²]
    d = 0.001 # 板间距离 [m]
    dE_dt = V / (d * t_charge)  # 电场变化率
    I_displacement = PhysicalConstants.ε₀ * A * dE_dt
    
    println("   板面积 A = $A m²")
    println("   板间距 d = $d m")
    println("   电场变化率 dE/dt = $(round(dE_dt, sigdigits=4)) V/(m⋅s)")
    println("   位移电流 I_D = ε₀A(dE/dt) = $(round(I_displacement, sigdigits=6)) A")
    
    # 电流连续性检查
    current_ratio = I_displacement / I_conduction
    println("   电流比值 I_D/I_conduction = $(round(current_ratio, sigdigits=4))")
    println("   注：在此例中位移电流远小于传导电流")
    
    return I_conduction, I_displacement
end

"""
法拉第电磁感应定律
Faraday's Law of Electromagnetic Induction
"""
function faraday_law_derivation()
    println("\n=== 法拉第电磁感应定律推导 ===")
    println("=== Faraday's Law Derivation ===")
    
    # 法拉第实验
    println("1. 法拉第实验 (1831年) / Faraday's Experiments:")
    println("   - 移动磁铁靠近线圈 → 产生电流")
    println("   - 改变线圈中的电流 → 邻近线圈产生电流")  
    println("   - 旋转线圈在磁场中 → 产生交流电")
    
    # 数学表达
    println("\n2. 数学表达 / Mathematical Expression:")
    println("   感应电动势：ε = -dΦ/dt")
    println("   磁通量：Φ = ∫∫ B · dA")
    println("   积分形式：∮ E · dl = -∫∫ ∂B/∂t · dA")
    println("   微分形式：∇ × E = -∂B/∂t")
    
    # 数值示例：正弦磁场中的线圈
    println("\n3. 数值示例 - 正弦磁场中的线圈:")
    
    # 参数设置
    B₀ = 0.1      # 磁场幅值 [T]
    f = 50        # 频率 [Hz]
    ω = 2π * f    # 角频率 [rad/s]
    A = 0.01      # 线圈面积 [m²]
    N = 100       # 线圈匝数
    
    println("   磁场幅值 B₀ = $B₀ T")
    println("   频率 f = $f Hz")
    println("   角频率 ω = $(round(ω, sigdigits=4)) rad/s")
    println("   线圈面积 A = $A m²")
    println("   线圈匝数 N = $N")
    
    # 时间函数
    t_range = 0:0.0001:0.04  # 2个周期
    
    # 磁通量和感应电动势
    Φ(t) = B₀ * A * cos(ω * t)
    ε(t) = -N * (-B₀ * A * ω * sin(ω * t))  # -dΦ/dt
    
    println("\n   磁通量：Φ(t) = B₀A cos(ωt)")
    println("   感应电动势：ε(t) = NB₀Aω sin(ωt)")
    println("   最大电动势：ε_max = $(round(N * B₀ * A * ω, sigdigits=4)) V")
    
    return collect(t_range), [Φ(t) for t in t_range], [ε(t) for t in t_range]
end

"""
电磁波推导
Electromagnetic Wave Derivation
"""
function electromagnetic_wave_derivation()
    println("\n=== 电磁波推导 ===")
    println("=== Electromagnetic Wave Derivation ===")
    
    # 真空中的麦克斯韦方程组
    println("1. 真空中的麦克斯韦方程组 (ρ=0, J=0):")
    println("   (1) ∇ · E = 0")
    println("   (2) ∇ × B = μ₀ε₀ ∂E/∂t")
    println("   (3) ∇ × E = -∂B/∂t")
    println("   (4) ∇ · B = 0")
    
    # 波动方程推导
    println("\n2. 波动方程推导:")
    println("   对方程(3)取旋度：∇ × (∇ × E) = ∇ × (-∂B/∂t)")
    println("   = -∂(∇ × B)/∂t = -μ₀ε₀ ∂²E/∂t²")
    println("   使用矢量恒等式：∇ × (∇ × E) = ∇(∇ · E) - ∇²E")
    println("   在真空中 ∇ · E = 0，所以：")
    println("   ∇²E = μ₀ε₀ ∂²E/∂t²")
    
    println("\n   这是标准波动方程：∇²f = (1/v²) ∂²f/∂t²")
    println("   波速：v = 1/√(μ₀ε₀) = c")
    
    # 数值验证
    println("\n3. 数值验证:")
    c_theory = 1/sqrt(PhysicalConstants.μ₀ * PhysicalConstants.ε₀)
    println("   理论计算：c = 1/√(μ₀ε₀) = $(round(c_theory, sigdigits=8)) m/s")
    println("   实验测量：c = $(PhysicalConstants.c) m/s")
    error_percent = abs(c_theory - PhysicalConstants.c)/PhysicalConstants.c * 100
    println("   相对误差：$(round(error_percent, sigdigits=4))%")
    
    # 平面波解
    println("\n4. 平面波解 (z方向传播):")
    
    # 参数设置
    E₀ = 1.0      # 电场幅值 [V/m]
    f = 1.0 * 10.0^9       # 频率 1GHz [Hz]
    ω = 2π * f    # 角频率
    k = ω / PhysicalConstants.c  # 波数
    λ = 2π / k    # 波长
    
    println("   频率 f = $(f) Hz")
    println("   角频率 ω = $(round(ω, sigdigits=4)) rad/s")
    println("   波数 k = $(round(k, sigdigits=4)) rad/m")
    println("   波长 λ = $(round(λ, sigdigits=4)) m")
    
    # 电场和磁场
    z_range = 0:λ/100:2λ
    t = 0  # t=0时刻
    
    Ex(z,t) = E₀ * cos(k*z - ω*t)  # x方向电场
    By(z,t) = (E₀/PhysicalConstants.c) * cos(k*z - ω*t)  # y方向磁场
    
    println("\n   电场：Ex(z,t) = E₀ cos(kz - ωt)")
    println("   磁场：By(z,t) = (E₀/c) cos(kz - ωt)")
    println("   E/B比值 = c = $(PhysicalConstants.c) m/s")
    
    # 能量密度
    ε₀ = PhysicalConstants.ε₀
    μ₀ = PhysicalConstants.μ₀
    
    u_E = 0.5 * ε₀ * E₀^2  # 电场能量密度
    u_B = 0.5 * (E₀/PhysicalConstants.c)^2 / μ₀  # 磁场能量密度
    
    println("\n5. 能量分析:")
    println("   电场能量密度：u_E = ½ε₀E² = $(round(u_E, sigdigits=4)) J/m³")
    println("   磁场能量密度：u_B = ½B²/μ₀ = $(round(u_B, sigdigits=4)) J/m³")
    println("   总能量密度：u = u_E + u_B = $(round(u_E + u_B, sigdigits=4)) J/m³")
    println("   能量密度相等验证：u_E ≈ u_B ? $(abs(u_E - u_B) < 1.0 * 10.0^(-15))")
    
    return collect(z_range), [Ex(z,t) for z in z_range], [By(z,t) for z in z_range]
end

"""
可视化函数 (使用英文标签避免字体问题)
Visualization Functions (Using English labels to avoid font issues)
"""
function plot_maxwell_demonstrations()
    println("\n=== 生成麦克斯韦方程组演示图表 (英文标签版) ===")
    
    # 1. 库仑定律到高斯定律
    r_vals, E_vals = coulomb_to_gauss_derivation()
    
    p1 = plot(r_vals, E_vals, 
              title="Point Charge Electric Field vs Distance",
              xlabel="Distance r [m]", ylabel="Electric Field E [V/m]",
              linewidth=2, marker=:circle, markersize=6,
              yscale=:log10, label="E ∝ 1/r²",
              legend=:topright)
    
    # 2. 安培定律修正
    I_cond, I_disp = ampere_to_maxwell_derivation()
    
    p2 = bar(["Conduction\nCurrent", "Displacement\nCurrent"], 
             [I_cond, I_disp],
             title="Capacitor Charging Current Comparison",
             ylabel="Current [A]", 
             color=[:blue, :red],
             alpha=0.7,
             legend=false)
    
    # 3. 法拉第定律
    t_vals, Φ_vals, ε_vals = faraday_law_derivation()
    
    p3 = plot(t_vals.*1000, [Φ_vals.*1000, ε_vals], 
              title="Faraday's Electromagnetic Induction",
              xlabel="Time [ms]", ylabel="Flux [mWb] / EMF [V]",
              linewidth=2, 
              label=["Magnetic Flux Φ(t)" "Induced EMF ε(t)"],
              color=[:blue :red])
    
    # 4. 电磁波
    z_vals, Ex_vals, By_vals = electromagnetic_wave_derivation()
    
    p4 = plot(z_vals.*1e9, [Ex_vals, By_vals.*PhysicalConstants.c], 
              title="Electromagnetic Wave (1 GHz)",
              xlabel="Position [nm]", ylabel="Field Strength [V/m]",
              linewidth=2,
              label=["Electric Field Ex" "Magnetic Field By×c"],
              color=[:blue :red])
    
    # 组合图表
    layout = @layout [a b; c d]
    combined_plot = plot(p1, p2, p3, p4, 
                        layout=layout,
                        size=(1000, 700),
                        plot_title="Maxwell's Equations Numerical Verification")
    
    return combined_plot
end

"""
完整的麦克斯韦方程组数值验证
Complete Maxwell Equations Numerical Verification
"""
function complete_maxwell_verification()
    println("="^70)
    println("麦克斯韦方程组完整推导与数值验证 (修复版)")
    println("Maxwell's Equations Complete Derivation and Numerical Verification (Fixed)")
    println("="^70)
    
    # 验证物理常数
    PhysicalConstants.verify_constants()
    
    # 各个定律的推导
    coulomb_to_gauss_derivation()
    ampere_to_maxwell_derivation()
    faraday_law_derivation()
    electromagnetic_wave_derivation()
    
    # 生成可视化
    plots = plot_maxwell_demonstrations()
    
    println("\n=== 总结 / Summary ===")
    println("✅ 验证了真空常数之间的关系：c = 1/√(μ₀ε₀)")
    println("✅ 推导了从库仑定律到高斯定律")
    println("✅ 展示了麦克斯韦对安培定律的修正")
    println("✅ 验证了法拉第电磁感应定律")
    println("✅ 推导了电磁波方程并验证波速")
    println("✅ 生成了数值验证图表 (英文标签)")
    
    println("\n麦克斯韦方程组 (微分形式):")
    println("(1) ∇ · E = ρ/ε₀             [Gauss's Law]")
    println("(2) ∇ × B = μ₀J + μ₀ε₀∂E/∂t  [Ampère-Maxwell Law]")
    println("(3) ∇ × E = -∂B/∂t           [Faraday's Law]")
    println("(4) ∇ · B = 0                [Gauss's Law for Magnetism]")
    
    return plots
end

# 主程序执行
if abspath(PROGRAM_FILE) == @__FILE__
    plots = complete_maxwell_verification()
    
    # 保存图表
    savefig(plots, "maxwell_equations_verification_fixed.png")
    println("\n图表已保存为: maxwell_equations_verification_fixed.png")
    
    # 显示图表
    display(plots)
end
