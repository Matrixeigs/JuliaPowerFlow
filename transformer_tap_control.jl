using Plots, LaTeXStrings, LinearAlgebra
plotlyjs()

"""
三绕组变压器模型 - 实际工程配置
"""
mutable struct ThreeWindingTransformer
    # 基本参数
    S_H::Float64        # 高压侧额定容量 (MVA)
    S_M::Float64        # 中压侧额定容量 (MVA)  
    S_L::Float64        # 低压侧额定容量 (MVA)
    V_H::Float64        # 高压侧额定电压 (kV)
    V_M::Float64        # 中压侧额定电压 (kV)
    V_L::Float64        # 低压侧额定电压 (kV)
    
    # 阻抗参数 (标幺值)
    R_H::Float64        # 高压侧电阻
    X_H::Float64        # 高压侧电抗
    R_M::Float64        # 中压侧电阻
    X_M::Float64        # 中压侧电抗
    R_L::Float64        # 低压侧电阻
    X_L::Float64        # 低压侧电抗
    
    # 分接头参数 - 实际工程配置
    tap_H::Float64      # 高压侧分接头位置 (有载分接头)
    tap_M::Float64      # 中压侧分接头位置 (可选配置)
    # 注意：低压侧通常无分接头，固定为1.0
    
    tap_H_range::Tuple{Float64,Float64}     # 高压侧分接头范围
    tap_M_range::Tuple{Float64,Float64}     # 中压侧分接头范围 (如果配置)
    tap_H_step::Float64                     # 高压侧分接头步长
    tap_M_step::Float64                     # 中压侧分接头步长
    
    # 分接头开关类型
    tap_H_type::String  # "OLTC" 有载分接头开关
    tap_M_type::String  # "Off-load" 无载分接头开关或 "None" 无分接头
    
    # 控制参数
    tap_H_enabled::Bool # 高压侧分接头是否启用
    tap_M_enabled::Bool # 中压侧分接头是否启用
end

"""
创建典型三绕组变压器 - 实际工程配置
"""
function create_typical_transformer()
    return ThreeWindingTransformer(
        180.0, 120.0, 60.0,     # 容量 MVA
        220.0, 110.0, 35.0,     # 电压 kV
        0.005, 0.08,            # 高压侧阻抗
        0.008, 0.12,            # 中压侧阻抗
        0.012, 0.15,            # 低压侧阻抗
        1.0, 1.0,               # 分接头位置
        (0.9, 1.1),             # 高压侧分接头范围 ±10%
        (0.95, 1.05),           # 中压侧分接头范围 ±5% (较小范围)
        0.0125,                 # 高压侧分接头步长 1.25%
        0.025,                  # 中压侧分接头步长 2.5%
        "OLTC",                 # 高压侧有载分接头开关
        "Off-load",             # 中压侧无载分接头开关
        true,                   # 高压侧分接头启用
        false                   # 中压侧分接头通常不启用自动控制
    )
end

"""
计算变压器导纳矩阵 - 考虑实际分接头配置
"""
function calculate_transformer_admittance(transformer::ThreeWindingTransformer)
    # 计算修正后的阻抗
    Z_H = (transformer.R_H + 1im * transformer.X_H) * transformer.tap_H^2
    Z_M = (transformer.R_M + 1im * transformer.X_M) * transformer.tap_M^2
    Z_L = (transformer.R_L + 1im * transformer.X_L)  # 低压侧无分接头，固定阻抗
    
    # 星形等效电路导纳
    Y_H = 1.0 / Z_H
    Y_M = 1.0 / Z_M
    Y_L = 1.0 / Z_L
    
    # 构建导纳矩阵
    Y = zeros(Complex{Float64}, 3, 3)
    Y[1,1] = Y_H + Y_M + Y_L  # 高压母线
    Y[2,2] = Y_M * transformer.tap_M^2  # 中压母线
    Y[3,3] = Y_L  # 低压母线 (固定)
    
    Y[1,2] = Y[2,1] = -Y_M * transformer.tap_M
    Y[1,3] = Y[3,1] = -Y_L  # 低压侧固定
    Y[2,3] = Y[3,2] = 0.0
    
    return Y
end

"""
分接头对电压的影响分析 - 实际工程配置
"""
function analyze_tap_voltage_effect()
    transformer = create_typical_transformer()
    
    # 分接头变化范围
    tap_H_values = transformer.tap_H_range[1]:transformer.tap_H_step:transformer.tap_H_range[2]
    tap_M_values = transformer.tap_M_range[1]:transformer.tap_M_step:transformer.tap_M_range[2]
    
    # 固定运行条件
    V_H = 1.05  # 高压侧电压 (标幺值)
    P_M = 0.8   # 中压侧有功负荷 (标幺值)
    Q_M = 0.3   # 中压侧无功负荷 (标幺值)
    P_L = 0.5   # 低压侧有功负荷 (标幺值)
    Q_L = 0.2   # 低压侧无功负荷 (标幺值)
    
    # 分析高压侧分接头影响 (主要控制手段)
    V_M_h = Float64[]
    V_L_h = Float64[]
    
    for tap in tap_H_values
        transformer.tap_H = tap
        
        # 简化潮流计算
        V_M = V_H / tap * transformer.V_M / transformer.V_H
        V_L = V_H * transformer.V_L / transformer.V_H  # 低压侧主要受高压侧影响
        
        push!(V_M_h, V_M)
        push!(V_L_h, V_L)
    end
    
    # 分析中压侧分接头影响 (辅助调节，通常为无励磁调节)
    transformer.tap_H = 1.0  # 恢复默认值
    V_M_m = Float64[]
    V_L_m = Float64[]
    
    for tap in tap_M_values
        transformer.tap_M = tap
        
        # 中压侧分接头主要影响中压侧电压
        V_M = V_H * transformer.V_M / transformer.V_H / tap
        V_L = V_H * transformer.V_L / transformer.V_H  # 低压侧受影响较小
        
        push!(V_M_m, V_M)
        push!(V_L_m, V_L)
    end
    
    p1 = plot(tap_H_values, V_M_h,
              title="高压侧OLTC对电压的影响",
              xlabel="高压侧分接头位置",
              ylabel="电压 (标幺值)",
              label="中压侧电压",
              linewidth=2,
              color=:blue)
    
    plot!(p1, tap_H_values, V_L_h,
          label="低压侧电压",
          linewidth=2,
          color=:red)
    
    hline!(p1, [1.05, 0.95], 
           label="电压限制", 
           color=:gray, 
           linestyle=:dash)
    
    p2 = plot(tap_M_values, V_M_m,
              title="中压侧分接头对电压的影响 (无励磁调节)",
              xlabel="中压侧分接头位置", 
              ylabel="电压 (标幺值)",
              label="中压侧电压",
              linewidth=2,
              color=:green)
    
    plot!(p2, tap_M_values, V_L_m,
          label="低压侧电压",
          linewidth=2,
          color=:orange)
    
    hline!(p2, [1.05, 0.95],
           label="电压限制",
           color=:gray,
           linestyle=:dash)
    
    return plot(p1, p2, layout=(1,2), size=(1200,400))
end

"""
实际工程电压控制策略 - 仅使用高压侧OLTC
"""
function practical_voltage_control_simulation()
    transformer = create_typical_transformer()
    
    # 时间序列（24小时）
    hours = 1:24
    
    # 负荷曲线（标幺值）
    load_curve = [0.6, 0.55, 0.5, 0.48, 0.5, 0.55, 0.65, 0.8, 
                  0.9, 0.95, 1.0, 1.05, 1.0, 0.95, 0.9, 0.85,
                  0.8, 0.85, 0.9, 0.95, 0.9, 0.8, 0.7, 0.65]
    
    # 高压侧电压变化
    V_H_curve = 1.05 .+ 0.02 * sin.(2π * hours / 24)
    
    # 目标电压 (中压侧)
    V_M_target = 1.0
    V_deadband = 0.02
    
    # 仿真变量
    V_M_uncontrolled = Float64[]
    V_M_controlled = Float64[]
    tap_H_positions = Float64[]
    current_tap_H = 1.0
    
    # OLTC控制参数
    tap_change_delay = 0  # 延时计数器
    min_delay_steps = 2   # 最小延时步数，防止频繁动作
    
    for (i, hour) in enumerate(hours)
        V_H = V_H_curve[i]
        load = load_curve[i]
        
        # 无控制情况下的中压侧电压
        V_M_uncon = V_H * transformer.V_M / transformer.V_H
        push!(V_M_uncontrolled, V_M_uncon)
        
        # 实际OLTC控制 - 仅使用高压侧分接头
        V_M_actual = V_H * transformer.V_M / transformer.V_H / current_tap_H
        
        # 电压偏差检测
        voltage_error = V_M_actual - V_M_target
        
        # OLTC控制逻辑
        tap_change_needed = false
        if abs(voltage_error) > V_deadband
            if tap_change_delay <= 0
                # 计算所需分接头调整
                if voltage_error > V_deadband
                    # 电压过高，增加分接头（降低二次侧电压）
                    new_tap = min(current_tap_H + transformer.tap_H_step, 
                                transformer.tap_H_range[2])
                    if new_tap != current_tap_H
                        current_tap_H = new_tap
                        tap_change_needed = true
                    end
                elseif voltage_error < -V_deadband
                    # 电压过低，减少分接头（提高二次侧电压）
                    new_tap = max(current_tap_H - transformer.tap_H_step,
                                transformer.tap_H_range[1])
                    if new_tap != current_tap_H
                        current_tap_H = new_tap
                        tap_change_needed = true
                    end
                end
                
                if tap_change_needed
                    tap_change_delay = min_delay_steps
                end
            else
                tap_change_delay -= 1
            end
        else
            tap_change_delay = 0  # 在死区内，重置延时
        end
        
        # 计算实际控制后的电压
        V_M_controlled_actual = V_H * transformer.V_M / transformer.V_H / current_tap_H
        push!(V_M_controlled, V_M_controlled_actual)
        push!(tap_H_positions, current_tap_H)
    end
    
    p1 = plot(hours, V_M_uncontrolled,
              title="实际工程电压控制效果 (仅高压侧OLTC)",
              xlabel="时间 (小时)",
              ylabel="中压侧电压 (标幺值)",
              label="无OLTC控制",
              linewidth=2,
              color=:red)
    
    plot!(p1, hours, V_M_controlled,
          label="OLTC控制",
          linewidth=2,
          color=:blue)
    
    hline!(p1, [V_M_target + V_deadband, V_M_target - V_deadband],
           label="死区边界",
           color=:gray,
           linestyle=:dash)
    
    hline!(p1, [V_M_target],
           label="目标电压", 
           color=:green,
           linestyle=:dot)
    
    p2 = plot(hours, tap_H_positions,
              title="高压侧OLTC分接头位置变化",
              xlabel="时间 (小时)",
              ylabel="分接头位置",
              label="高压侧OLTC",
              linewidth=2,
              color=:purple,
              marker=:circle)
    
    hline!(p2, [transformer.tap_H_range[1], transformer.tap_H_range[2]],
           label="分接头限制",
           color=:red,
           linestyle=:dash)
    
    return plot(p1, p2, layout=(2,1), size=(800,600))
end

"""
分接头调节对网损的影响
"""
function tap_loss_analysis()
    transformer = create_typical_transformer()
    
    tap_values = 0.9:0.005:1.1
    losses = Float64[]
    voltages = Float64[]
    
    # 固定负荷条件
    S_load = 0.8 + 0.3im  # 负荷功率 (标幺值)
    V_H = 1.05           # 高压侧电压
    
    for tap in tap_values
        transformer.tap_M = tap
        
        # 计算中压侧电压
        V_M = V_H * transformer.V_M / transformer.V_H / tap
        
        # 计算电流
        I_M = abs(S_load) / V_M
        
        # 计算损耗
        R_total = transformer.R_M * tap^2
        P_loss = I_M^2 * R_total
        
        push!(losses, P_loss)
        push!(voltages, V_M)
    end
    
    p1 = plot(tap_values, losses,
              title="分接头调节对变压器损耗的影响",
              xlabel="中压侧分接头位置",
              ylabel="变压器损耗 (标幺值)",
              linewidth=2,
              color=:red,
              label="变压器损耗")
    
    p2 = plot(tap_values, voltages,
              title="分接头调节对电压的影响",
              xlabel="中压侧分接头位置",
              ylabel="中压侧电压 (标幺值)",
              linewidth=2,
              color=:blue,
              label="中压侧电压")
    
    hline!(p2, [1.05, 0.95],
           label="电压限制",
           color=:gray,
           linestyle=:dash)
    
    return plot(p1, p2, layout=(1,2), size=(1000,400))
end

"""
多目标优化：电压质量与网损平衡
"""
function multi_objective_optimization()
    transformer = create_typical_transformer()
    
    tap_values = 0.9:0.01:1.1
    V_target = 1.0
    
    # 权重系数
    w_voltage = 0.7
    w_loss = 0.3
    
    objective_values = Float64[]
    voltage_deviations = Float64[]
    loss_values = Float64[]
    
    for tap in tap_values
        transformer.tap_M = tap
        
        # 电压偏差
        V_M = 1.05 * transformer.V_M / transformer.V_H / tap
        voltage_dev = abs(V_M - V_target)
        
        # 损耗
        I_M = 0.8 / V_M  # 简化电流计算
        loss = I_M^2 * transformer.R_M * tap^2
        
        # 多目标函数
        objective = w_voltage * voltage_dev + w_loss * loss
        
        push!(objective_values, objective)
        push!(voltage_deviations, voltage_dev)
        push!(loss_values, loss)
    end
    
    # 找到最优分接头位置
    optimal_idx = argmin(objective_values)
    optimal_tap = tap_values[optimal_idx]
    
    p1 = plot(tap_values, objective_values,
              title="多目标优化",
              xlabel="分接头位置",
              ylabel="目标函数值",
              linewidth=2,
              color=:purple,
              label="综合目标函数")
    
    scatter!(p1, [optimal_tap], [objective_values[optimal_idx]],
             markersize=8,
             color=:red,
             label="最优解")
    
    p2 = plot(tap_values, voltage_deviations,
              title="目标分量",
              xlabel="分接头位置",
              ylabel="标幺值",
              linewidth=2,
              color=:blue,
              label="电压偏差")
    
    plot!(p2, tap_values, loss_values,
          linewidth=2,
          color=:red,
          label="变压器损耗")
    
    return plot(p1, p2, layout=(1,2), size=(1000,400))
end

"""
运行所有变压器分接头控制分析
"""
function run_transformer_tap_analysis()
    println("="^60)
    println("三绕组变压器分接头电压控制分析 (实际工程配置)")
    println("="^60)
    
    # 1. 创建示例变压器
    transformer = create_typical_transformer()
    println("\n1. 变压器参数 (实际工程配置):")
    println("   容量: $(transformer.S_H)/$(transformer.S_M)/$(transformer.S_L) MVA")
    println("   电压: $(transformer.V_H)/$(transformer.V_M)/$(transformer.V_L) kV")
    println("   高压侧分接头: $(transformer.tap_H_type), 范围: $(transformer.tap_H_range)")
    println("   中压侧分接头: $(transformer.tap_M_type), 范围: $(transformer.tap_M_range)")
    println("   低压侧: 固定变比，无分接头")
    println("   高压侧步长: $(transformer.tap_H_step), 中压侧步长: $(transformer.tap_M_step)")
    
    # 2. 生成分析图形
    analysis_plots = [
        ("分接头电压影响 (实际配置)", analyze_tap_voltage_effect()),
        ("实际工程电压控制", practical_voltage_control_simulation()),
        ("分接头损耗分析", tap_loss_analysis()),
        ("多目标优化", multi_objective_optimization())
    ]
    
    # 显示所有图形
    for (title, plot_obj) in analysis_plots
        println("\n生成图形: $title")
        display(plot_obj)
    end
    
    println("\n" * "="^60)
    println("实际工程应用总结:")
    println("✓ 高压侧OLTC是主要的电压控制手段")
    println("✓ 中压侧分接头通常为无载调节，用于季节性或计划性调整")
    println("✓ 低压侧通常为固定变比，无分接头开关")
    println("✓ OLTC控制需考虑延时和死区，避免频繁动作")
    println("✓ 实际工程中优先使用高压侧OLTC进行电压调节")
    println("✓ 注意：励磁调节是发电机控制概念，与变压器分接头调节不同")
    println("="^60)
    
    return analysis_plots
end
