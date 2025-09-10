"""
电磁感应理论示例程序
Electromagnetic Induction Theory Examples

本程序演示了电磁感应理论的各种应用场景，包括：
1. 基础电磁感应现象
2. 发电机原理
3. 电动机原理
4. 变压器工作原理
5. 感应加热原理

运行此程序可以看到完整的理论分析和可视化结果
"""

using Pkg

# 检查并安装必要的包
required_packages = ["Plots", "DifferentialEquations", "LinearAlgebra"]
for package in required_packages
    try
        eval(Meta.parse("using $package"))
    catch
        println("正在安装 $package...")
        Pkg.add(package)
        eval(Meta.parse("using $package"))
    end
end

# 导入我们的电磁感应模块
include("ElectromagneticInduction.jl")
using .ElectromagneticInduction

using Plots
using LinearAlgebra
using DifferentialEquations

# 配置绘图后端以支持中文或使用英文标签
function setup_plotting()
    try
        # 尝试使用支持中文的字体
        Plots.default(fontfamily="SimHei")
    catch
        try
            # 如果SimHei不可用，尝试其他中文字体
            Plots.default(fontfamily="Microsoft YaHei")
        catch
            # 如果都不可用，使用默认字体并输出警告
            println("⚠️  中文字体不可用，将使用英文标签")
            Plots.default(fontfamily="sans-serif")
        end
    end
end

"""
主演示函数
"""
function main_electromagnetic_induction_demo()
    println("🔌 Electromagnetic Induction Theory Complete Demo")
    println("🔌 电磁感应理论完整演示程序")
    println("=" ^ 50)
    
    # 设置绘图配置
    setup_plotting()
    
    # 1. 打印数学推导
    ElectromagneticInduction.print_mathematical_derivations()
    
    println("\n" * "=" ^ 50)
    println("开始仿真演示...")
    println("=" ^ 50)
    
    # 2. 运行完整仿真
    simulation_results = ElectromagneticInduction.simulate_electromagnetic_induction()
    
    # 3. 生成可视化图表
    try
        println("生成分析图表...")
        plots = ElectromagneticInduction.plot_induction_analysis(simulation_results)
        display(plots)
        
        # 保存图表
        savefig(plots, "electromagnetic_induction_analysis.png")
        println("✅ 图表已保存为 electromagnetic_induction_analysis.png")
    catch e
        println("⚠️  图表生成失败: $e")
        println("请确保已安装 Plots 包")
    end
    
    return simulation_results
end

"""
发电机原理演示
"""
function generator_principle_demo()
    println("\n🔄 发电机原理演示")
    println("=" * 30)
    
    println("理论基础：")
    println("当导体在磁场中运动时，根据洛伦兹力定律：")
    println("F⃗ = q(v⃗ × B⃗)")
    println("对于导体中的自由电荷，会产生电动势：")
    println("ε = B × l × v")
    println("其中 B=磁感应强度, l=导体长度, v=运动速度")
    println()
    
    # 参数设置
    B = 1.0        # 磁感应强度 1T
    l = 0.5        # 导体长度 0.5m
    ω = 2π * 50    # 角频率 50Hz
    t = 0:0.001:0.04  # 时间 0-40ms
    
    # 旋转发电机
    θ = ω * t
    v = ω * 0.1    # 线速度（假设半径0.1m）
    ε = B * l * v * sin.(θ)  # 感应电动势
    
    println("仿真参数：")
    println("磁感应强度 B = $B T")
    println("导体长度 l = $l m")
    println("旋转频率 f = 50 Hz")
    println("最大电动势 = $(round(maximum(abs.(ε)), digits=3)) V")
    
    # 绘图
    try
        p1 = plot(t*1000, ε, title="Generator Induced EMF", xlabel="Time (ms)", 
                 ylabel="EMF (V)", linewidth=2, color=:red, label="ε(t)")
        
        p2 = plot(t*1000, θ*180/π, title="Rotor Angle", xlabel="Time (ms)", 
                 ylabel="Angle (°)", linewidth=2, color=:blue, label="θ(t)")
        
        generator_plot = plot(p1, p2, layout=(2,1), size=(600, 400))
        display(generator_plot)
        savefig(generator_plot, "generator_principle.png")
        println("✅ 发电机原理图已保存")
    catch e
        println("⚠️  发电机图表生成失败: $e")
    end
    
    return t, ε, θ
end

"""
电动机原理演示
"""
function motor_principle_demo()
    println("\n⚡ 电动机原理演示")
    println("=" * 30)
    
    println("理论基础：")
    println("载流导体在磁场中受到安培力：")
    println("F⃗ = I × l⃗ × B⃗")
    println("力的大小：F = B × I × l × sin(θ)")
    println("产生的转矩：τ = F × r = B × I × l × r")
    println()
    
    # 参数设置
    B = 0.5      # 磁感应强度
    I = 10.0     # 电流
    l = 0.3      # 导体长度
    r = 0.05     # 转子半径
    
    # 计算力和转矩
    F_max = B * I * l
    τ_max = F_max * r
    
    println("计算结果：")
    println("最大安培力 F_max = $(F_max) N")
    println("最大转矩 τ_max = $(round(τ_max, digits=4)) N·m")
    
    # 转矩随角度变化
    θ_motor = 0:0.1:2π
    τ = τ_max * sin.(θ_motor)
    
    try
        motor_plot = plot(θ_motor*180/π, τ, title="Motor Torque Characteristics", 
                         xlabel="Rotor Angle (°)", ylabel="Torque (N·m)", 
                         linewidth=2, color=:green, label="τ(θ)")
        display(motor_plot)
        savefig(motor_plot, "motor_principle.png")
        println("✅ 电动机原理图已保存")
    catch e
        println("⚠️  电动机图表生成失败: $e")
    end
    
    return θ_motor, τ, F_max, τ_max
end

"""
感应加热原理演示
"""
function induction_heating_demo()
    println("\n🔥 感应加热原理演示")
    println("=" * 30)
    
    println("理论基础：")
    println("涡流损耗功率：P = σ × (∂A⃗/∂t)² × V")
    println("其中：σ=电导率, A⃗=矢量磁势, V=体积")
    println("简化形式：P = σ × E² × V")
    println("焦耳热：Q = P × t = I²Rt")
    println()
    
    # 参数设置
    σ = 5.96e7     # 铜的电导率 S/m
    f = 100e3      # 感应频率 100kHz
    B0 = 0.1       # 磁感应强度幅值
    V_material = 1e-6  # 材料体积 1cm³
    
    # 感应电场强度（麦克斯韦方程）
    E = 2π * f * B0 * 0.01  # 简化计算，假设特征长度1cm
    
    # 功率密度
    P_density = σ * E^2
    P_total = P_density * V_material
    
    println("计算参数：")
    println("工作频率 f = $(f/1000) kHz")
    println("材料电导率 σ = $(σ/1e6) MS/m")
    println("感应电场强度 E = $(round(E, digits=2)) V/m")
    println("功率密度 = $(round(P_density/1e6, digits=2)) MW/m³")
    println("总功率 = $(round(P_total, digits=4)) W")
    
    # 温度上升模拟（简化）
    c_p = 385      # 铜的比热容 J/(kg·K)
    ρ = 8960       # 铜的密度 kg/m³
    mass = ρ * V_material
    
    t_heat = 0:0.1:10  # 加热时间 10秒
    T_rise = (P_total .* t_heat) ./ (mass * c_p)
    
    try
        heating_plot = plot(t_heat, T_rise, title="Induction Heating Temperature Rise", 
                           xlabel="Time (s)", ylabel="Temperature Rise (K)", 
                           linewidth=2, color=:orange, label="ΔT(t)")
        display(heating_plot)
        savefig(heating_plot, "induction_heating.png")
        println("✅ 感应加热图已保存")
    catch e
        println("⚠️  感应加热图表生成失败: $e")
    end
    
    return t_heat, T_rise, P_total
end

"""
变压器详细分析
"""
function detailed_transformer_analysis()
    println("\n🔌 变压器详细分析")
    println("=" * 30)
    
    # 创建理想变压器模型
    primary = ElectromagneticInduction.InductionCoil(2000, 0.01, 0.2, 1000*4π*1e-7, 2.0)
    secondary = ElectromagneticInduction.InductionCoil(400, 0.01, 0.2, 1000*4π*1e-7, 0.4)
    transformer = ElectromagneticInduction.Transformer(primary, secondary, 0.95, 1000*4π*1e-7)
    
    println("变压器参数：")
    println("初级：$(primary.turns)匝，电阻=$(primary.resistance)Ω，自感=$(round(primary.inductance*1000,digits=2))mH")
    println("次级：$(secondary.turns)匝，电阻=$(secondary.resistance)Ω，自感=$(round(secondary.inductance*1000,digits=2))mH")
    println("互感：$(round(transformer.mutual_inductance*1000,digits=2))mH")
    println("耦合系数：$(transformer.coupling_coefficient)")
    println()
    
    # 频率响应分析
    f_range = 10.0:10.0:1000.0  # 10Hz到1kHz
    ω_range = 2π * f_range
    
    # 计算阻抗
    Z1 = sqrt.(primary.resistance^2 .+ (ω_range * primary.inductance).^2)
    Z2 = sqrt.(secondary.resistance^2 .+ (ω_range * secondary.inductance).^2)
    
    # 变压器传输特性
    V_primary = 220.0  # 初级电压220V
    turns_ratio = secondary.turns / primary.turns
    
    # 理想变压器次级电压
    V_secondary_ideal = V_primary * turns_ratio
    
    # 考虑阻抗的实际次级电压（简化计算）
    efficiency = 0.95 * exp.(-f_range / 1000)  # 频率越高效率越低
    V_secondary_actual = V_secondary_ideal * efficiency
    
    println("变压器性能：")
    println("变压比：$(round(turns_ratio, digits=3))")
    println("理想次级电压：$(round(V_secondary_ideal, digits=1))V")
    println("50Hz时实际次级电压：$(round(V_secondary_actual[5], digits=1))V")
    println("50Hz时效率：$(round(efficiency[5]*100, digits=1))%")
    
    try
        p1 = plot(f_range, V_secondary_actual, title="Transformer Frequency Response", 
                 xlabel="Frequency (Hz)", ylabel="Secondary Voltage (V)", 
                 linewidth=2, color=:purple, label="V₂(f)")
        
        p2 = plot(f_range, efficiency*100, title="Transformer Efficiency", 
                 xlabel="Frequency (Hz)", ylabel="Efficiency (%)", 
                 linewidth=2, color=:blue, label="η(f)")
        
        transformer_plot = plot(p1, p2, layout=(2,1), size=(600, 400))
        display(transformer_plot)
        savefig(transformer_plot, "transformer_analysis.png")
        println("✅ 变压器分析图已保存")
    catch e
        println("⚠️  变压器图表生成失败: $e")
    end
    
    return f_range, V_secondary_actual, efficiency
end

"""
电磁感应应用总结
"""
function electromagnetic_applications_summary()
    println("\n📋 电磁感应应用总结")
    println("=" * 40)
    
    println("1. 发电机 (Generator)")
    println("   - 原理：机械能 → 电能")
    println("   - 应用：火电、水电、风电等")
    println("   - 关键：ε = BLv或ε = -dΦ/dt")
    println()
    
    println("2. 电动机 (Motor)")
    println("   - 原理：电能 → 机械能")
    println("   - 应用：工业设备、交通工具")
    println("   - 关键：F = BIL，τ = Fr")
    println()
    
    println("3. 变压器 (Transformer)")
    println("   - 原理：互感耦合")
    println("   - 应用：电力传输、电压变换")
    println("   - 关键：V₂/V₁ = N₂/N₁")
    println()
    
    println("4. 感应加热 (Induction Heating)")
    println("   - 原理：涡流损耗")
    println("   - 应用：金属熔炼、感应炉")
    println("   - 关键：P = σE²V")
    println()
    
    println("5. 电磁制动 (Electromagnetic Braking)")
    println("   - 原理：反向感应力矩")
    println("   - 应用：电动汽车、电梯")
    println("   - 关键：楞次定律阻碍运动")
    println()
    
    println("6. 无线充电 (Wireless Charging)")
    println("   - 原理：磁耦合共振")
    println("   - 应用：手机、电动汽车充电")
    println("   - 关键：高频交变磁场")
    println()
end

"""
主程序入口
"""
function run_complete_demonstration()
    println("🚀 启动电磁感应理论完整演示程序")
    println("=" ^ 60)
    
    try
        # 运行各个演示模块
        println("1️⃣  基础理论仿真")
        sim_results = main_electromagnetic_induction_demo()
        
        println("\n2️⃣  发电机原理")
        gen_results = generator_principle_demo()
        
        println("\n3️⃣  电动机原理")
        motor_results = motor_principle_demo()
        
        println("\n4️⃣  感应加热原理")
        heating_results = induction_heating_demo()
        
        println("\n5️⃣  变压器详细分析")
        transformer_results = detailed_transformer_analysis()
        
        println("\n6️⃣  应用总结")
        electromagnetic_applications_summary()
        
        println("\n" * "=" * 60)
        println("🎉 电磁感应理论演示完成！")
        println("所有图表已保存到当前目录")
        println("=" * 60)
        
        return (
            simulation = sim_results,
            generator = gen_results,
            motor = motor_results,
            heating = heating_results,
            transformer = transformer_results
        )
        
    catch e
        println("❌ 演示过程中发生错误：$e")
        println("请检查依赖包是否正确安装")
        return nothing
    end
end

# 如果直接运行此文件，执行完整演示
if abspath(PROGRAM_FILE) == @__FILE__
    results = run_complete_demonstration()
end
