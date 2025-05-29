using Plots, LaTeXStrings
plotlyjs()

"""
分布参数传输线模型 - 实现理论文档中的完整模型
"""
struct TransmissionLine
    R::Float64      # 电阻 Ω/km
    L::Float64      # 电感 H/km  
    C::Float64      # 电容 F/km
    G::Float64      # 电导 S/km
    length::Float64 # 长度 km
    freq::Float64   # 频率 Hz
end

"""
计算传输线特征参数
"""
function calculate_line_characteristics(line::TransmissionLine)
    ω = 2π * line.freq
    
    # 阻抗和导纳
    Z = line.R + 1im * ω * line.L  # Ω/km
    Y = line.G + 1im * ω * line.C  # S/km
    
    # 特征阻抗
    Z₀ = sqrt(Z / Y)
    
    # 传播常数
    γ = sqrt(Z * Y)
    α = real(γ)  # 衰减常数
    β = imag(γ)  # 相位常数
    
    # 波长和波速
    λ = 2π / β  # km
    v = ω / β   # km/s
    
    return (Z₀=Z₀, γ=γ, α=α, β=β, λ=λ, v=v, Z=Z, Y=Y)
end

"""
计算surge impedance loading (自然功率)
"""
function calculate_surge_impedance_loading(voltage_kv::Union{Int64,Float64}, Z₀::Complex)
    V = Float64(voltage_kv) * 1000
    SIL = V^2 / real(Z₀) / 1e6  # MW
    return SIL
end

"""
分布参数模型的功率传输分析
"""
function distributed_parameter_analysis(line::TransmissionLine, V₁::Float64, V₂::Float64, δ::Float64)
    chars = calculate_line_characteristics(line)
    γl = chars.γ * line.length
    
    # ABCD参数
    A = cosh(γl)
    B = chars.Z₀ * sinh(γl)
    C = sinh(γl) / chars.Z₀
    D = cosh(γl)
    
    # 功率计算
    V₁_complex = V₁ * exp(1im * δ)
    V₂_complex = V₂ + 0im
    
    # 发送端电流
    I₁ = (A * V₂_complex - V₁_complex) / (-B)
    
    # 功率
    S₁ = V₁_complex * conj(I₁)
    P₁ = real(S₁) / 1e6  # MW
    Q₁ = imag(S₁) / 1e6  # MVAr
    
    S₂ = V₂_complex * conj((V₁_complex - A * V₂_complex) / B)
    P₂ = real(S₂) / 1e6  # MW  
    Q₂ = imag(S₂) / 1e6  # MVAr
    
    # 损耗
    P_loss = P₁ - P₂
    Q_loss = Q₁ - Q₂
    
    return (P₁=P₁, Q₁=Q₁, P₂=P₂, Q₂=Q₂, P_loss=P_loss, Q_loss=Q_loss, chars=chars)
end

"""
临界长度分析 - 理论文档中的关键概念
"""
function analyze_critical_length()
    # 典型500kV线路参数
    voltages = [110, 220, 500, 750]  # kV
    distances = 0:10:800  # km
    
    results = Dict()
    
    for V in voltages
        # 根据电压等级设置参数
        if V == 110
            line_params = (R=0.1, L=1.2e-3, C=12e-9, G=0.0)
        elseif V == 220  
            line_params = (R=0.05, L=1.0e-3, C=13e-9, G=0.0)
        elseif V == 500
            line_params = (R=0.025, L=0.8e-3, C=14e-9, G=0.0)
        else  # 750kV
            line_params = (R=0.02, L=0.7e-3, C=15e-9, G=0.0)
        end
        
        charging_powers = Float64[]
        SILs = Float64[]
        
        for d in distances
            line = TransmissionLine(line_params.R, line_params.L, line_params.C, 
                                  line_params.G, d, 50.0)
            chars = calculate_line_characteristics(line)
            
            # 充电功率
            V_actual = Float64(V) * 1000
            Q_charging = V_actual^2 * imag(chars.Y) * d / 1e6
            push!(charging_powers, Q_charging)
            
            # 自然功率 - 使用Float64转换
            SIL = calculate_surge_impedance_loading(Float64(V), chars.Z₀)
            push!(SILs, SIL)
        end
        
        # 找到临界长度 (充电功率 = 自然功率)
        critical_idx = findfirst(i -> charging_powers[i] >= SILs[1], 1:length(distances))
        critical_length = critical_idx !== nothing ? distances[critical_idx] : NaN
        
        results[V] = (distances=distances, charging_powers=charging_powers, 
                     SIL=SILs[1], critical_length=critical_length)
    end
    
    return results
end

"""
绘制临界长度分析图
"""
function plot_critical_length_analysis()
    results = analyze_critical_length()
    
    p = plot(title="传输线临界长度分析", 
             xlabel="距离 (km)", 
             ylabel="充电功率 (MVAr)",
             legend=:topright)
    
    colors = [:red, :blue, :green, :orange]
    
    for (i, (V, data)) in enumerate(results)
        plot!(p, data.distances, data.charging_powers,
              label="$(V)kV 充电功率",
              linewidth=2, 
              color=colors[i])
        
        # 添加自然功率水平线
        hline!(p, [data.SIL], 
               label="$(V)kV SIL ($(round(data.SIL, digits=1))MW)",
               color=colors[i],
               linestyle=:dash)
        
        # 标记临界点
        if !isnan(data.critical_length)
            vline!(p, [data.critical_length],
                   label="$(V)kV 临界长度: $(data.critical_length)km",
                   color=colors[i],
                   linestyle=:dot)
        end
    end
    
    return p
end

"""
电压稳定性分析 - P-V和Q-V特性曲线
"""
function voltage_stability_analysis()
    # 系统参数
    E = 1.0      # 发电机电势 (标幺值)
    X_sys = 0.1  # 系统电抗 (标幺值)
    X_line_values = [0.2, 0.5, 0.8, 1.0]  # 不同线路电抗
    
    # 有功功率负荷水平
    P_loads = [0.3, 0.5, 0.7, 0.9]  # 标幺值
    
    p1 = plot(title="电压稳定性分析 - P-V特性",
              xlabel="接收端电压 (标幺值)",
              ylabel="有功功率 (标幺值)")
    
    p2 = plot(title="电压稳定性分析 - Q-V特性", 
              xlabel="接收端电压 (标幺值)",
              ylabel="无功功率 (标幺值)")
    
    colors = [:red, :blue, :green, :orange]
    
    for (i, P_load) in enumerate(P_loads)
        X_total = X_sys + X_line_values[min(i, length(X_line_values))]
        
        # 电压范围 - 从低压到高压
        voltages = 0.3:0.001:1.1
        P_values = Float64[]
        Q_values = Float64[]
        valid_voltages = Float64[]
        
        for V in voltages
            # 对于给定的P_load，计算对应的功角和无功功率
            # P = (E*V/X)*sin(δ)，求解δ
            sin_delta = P_load * X_total / (E * V)
            
            # 检查是否有解（sin_delta <= 1）
            if abs(sin_delta) <= 1.0
                delta = asin(sin_delta)
                cos_delta = cos(delta)
                
                # 计算无功功率：Q = (E²/X) - (E*V/X)*cos(δ)
                Q = (E^2 / X_total) - (E * V / X_total) * cos_delta
                
                push!(P_values, P_load)
                push!(Q_values, Q)
                push!(valid_voltages, V)
            end
        end
        
        if !isempty(valid_voltages)
            # P-V曲线（水平线，因为P是恒定的）
            plot!(p1, valid_voltages, P_values,
                  label="P = $(P_load) p.u.",
                  linewidth=2,
                  color=colors[i])
                  
            # Q-V曲线（鼻子形状）
            plot!(p2, valid_voltages, Q_values,
                  label="P = $(P_load) p.u.", 
                  linewidth=2,
                  color=colors[i])
            
            # 标记电压崩溃点（Q的最小值点）
            min_Q_idx = argmin(Q_values)
            critical_V = valid_voltages[min_Q_idx]
            critical_Q = Q_values[min_Q_idx]
            
            scatter!(p2, [critical_V], [critical_Q],
                    markersize=4,
                    color=colors[i],
                    label="")
        end
    end
    
    # 标记正常运行区域
    vline!(p1, [0.95], label="额定电压", color=:black, linestyle=:dash)
    vline!(p2, [0.95], label="额定电压", color=:black, linestyle=:dash)
    
    # 标记电压稳定极限区域
    vspan!(p1, [0.7, 0.9], alpha=0.2, color=:red, label="电压不稳定区域")
    vspan!(p2, [0.7, 0.9], alpha=0.2, color=:red, label="电压不稳定区域")
    
    return plot(p1, p2, layout=(1,2), size=(1000,400))
end

"""
ZIP负荷模型 - 无功功率特性
"""
struct ZIPLoad
    Q0::Float64     # 基准无功功率 (标幺值)
    V0::Float64     # 基准电压 (标幺值)
    α_z::Float64    # 恒阻抗比例
    α_i::Float64    # 恒电流比例  
    α_p::Float64    # 恒功率比例
end

"""
计算ZIP负荷的无功功率
Q(V) = Q0 * [α_z*(V/V0)² + α_i*(V/V0) + α_p]
"""
function calculate_zip_reactive_power(load::ZIPLoad, V::Float64)
    V_ratio = V / load.V0
    Q = load.Q0 * (load.α_z * V_ratio^2 + load.α_i * V_ratio + load.α_p)
    return Q
end

"""
电压稳定性分析 - Q-V特性曲线（Q为横坐标，V为纵坐标）
包含ZIP负荷特征
"""
function voltage_stability_analysis()
    # 系统参数
    E = 1.05     # 发电机电势 (标幺值)
    X_values = [0.3, 0.5, 0.8, 1.0]  # 不同系统电抗
    
    # 有功功率水平
    P_load = 0.6  # 恒定有功负荷 (标幺值)
    
    p1 = plot(title="电压稳定性分析 - Q-V特性曲线",
              xlabel="无功功率 Q (标幺值)",
              ylabel="电压 V (标幺值)",
              legend=:bottomright)
    
    colors = [:blue, :red, :green, :orange]
    
    # 绘制系统供电能力曲线
    for (i, X) in enumerate(X_values)
        Q_values = Float64[]
        V_values = Float64[]
        
        # 电压范围
        for V in 0.3:0.001:1.1
            # 对于给定的V和P，计算系统能提供的Q
            # P = (E*V/X)*sin(δ)，求解δ
            sin_delta = P_load * X / (E * V)
            
            if abs(sin_delta) <= 1.0
                delta = asin(sin_delta)
                cos_delta = cos(delta)
                
                # Q = (E²/X) - (E*V/X)*cos(δ)
                Q = (E^2 / X) - (E * V / X) * cos_delta
                
                push!(Q_values, Q)
                push!(V_values, V)
            end
        end
        
        if !isempty(Q_values)
            plot!(p1, Q_values, V_values,
                  label="系统供电能力 X=$(X) p.u.",
                  linewidth=2,
                  color=colors[i])
        end
    end
    
    # 绘制不同ZIP负荷特性
    zip_loads = [
        ZIPLoad(0.3, 1.0, 1.0, 0.0, 0.0),  # 纯恒阻抗负荷
        ZIPLoad(0.3, 1.0, 0.0, 1.0, 0.0),  # 纯恒电流负荷
        ZIPLoad(0.3, 1.0, 0.0, 0.0, 1.0),  # 纯恒功率负荷
        ZIPLoad(0.3, 1.0, 0.5, 0.3, 0.2)   # 混合负荷
    ]
    
    zip_labels = ["恒阻抗负荷 (Z)", "恒电流负荷 (I)", "恒功率负荷 (P)", "混合负荷 (ZIP)"]
    zip_styles = [:dash, :dot, :dashdot, :solid]
    
    for (i, (load, label, style)) in enumerate(zip(zip_loads, zip_labels, zip_styles))
        Q_load_values = Float64[]
        V_load_values = Float64[]
        
        for V in 0.5:0.01:1.1
            Q_load = calculate_zip_reactive_power(load, V)
            push!(Q_load_values, Q_load)
            push!(V_load_values, V)
        end
        
        plot!(p1, Q_load_values, V_load_values,
              label=label,
              linewidth=2,
              color=:black,
              linestyle=style)
    end
    
    # 标记额定运行点
    scatter!(p1, [0.3], [1.0], 
             markersize=6, 
             color=:red, 
             markershape=:star,
             label="额定运行点")
    
    # 添加参考线
    hline!(p1, [0.95, 1.05], 
           label="电压运行范围", 
           color=:gray, 
           linestyle=:dot)
    
    return p1
end

"""
详细的Q-V鼻子曲线分析 - 改为Q为横坐标
"""
function detailed_qv_nose_analysis()
    # 系统参数
    E = 1.05     # 发电机电势 (标幺值)
    X = 0.6      # 系统电抗 (标幺值)
    
    # 不同有功功率水平
    P_loads = [0.2, 0.4, 0.6, 0.8]
    
    p = plot(title="Q-V鼻子曲线详细分析（含ZIP负荷）",
             xlabel="无功功率 Q (标幺值)",
             ylabel="电压 V (标幺值)",
             legend=:bottomright)
    
    colors = [:blue, :red, :green, :orange]
    
    # 绘制系统Q-V特性
    for (i, P) in enumerate(P_loads)
        Q_values = Float64[]
        V_values = Float64[]
        
        # 扫描电压范围
        for V in 0.2:0.001:1.2
            sin_delta = P * X / (E * V)
            
            if abs(sin_delta) <= 1.0
                delta = asin(sin_delta)
                cos_delta = cos(delta)
                Q = (E^2 / X) - (E * V / X) * cos_delta
                
                push!(Q_values, Q)
                push!(V_values, V)
            end
        end
        
        if !isempty(Q_values)
            plot!(p, Q_values, V_values,
                  label="系统特性 P=$(P) p.u.",
                  linewidth=2,
                  color=colors[i])
            
            # 标记鼻子点（V的最小值点）
            if length(V_values) > 10
                min_idx = argmin(V_values)
                nose_Q = Q_values[min_idx]
                nose_V = V_values[min_idx]
                
                scatter!(p, [nose_Q], [nose_V],
                        markersize=5,
                        color=colors[i],
                        markershape=:diamond,
                        label="")
            end
        end
    end
    
    # 添加ZIP负荷特性
    zip_load_mixed = ZIPLoad(0.4, 1.0, 0.4, 0.4, 0.2)  # 典型混合负荷
    
    Q_load_values = Float64[]
    V_load_values = Float64[]
    
    for V in 0.3:0.01:1.1
        Q_load = calculate_zip_reactive_power(zip_load_mixed, V)
        push!(Q_load_values, Q_load)
        push!(V_load_values, V)
    end
    
    plot!(p, Q_load_values, V_load_values,
          label="ZIP负荷特性",
          linewidth=3,
          color=:black,
          linestyle=:dash)
    
    # 计算交点（运行点）
    # 简化计算：找最接近的点
    intersections_Q = Float64[]
    intersections_V = Float64[]
    
    for (i, P) in enumerate(P_loads)
        # 重新计算该P水平下的系统特性
        sys_Q = Float64[]
        sys_V = Float64[]
        
        for V in 0.5:0.01:1.1
            sin_delta = P * X / (E * V)
            if abs(sin_delta) <= 1.0
                delta = asin(sin_delta)
                cos_delta = cos(delta)
                Q = (E^2 / X) - (E * V / X) * cos_delta
                push!(sys_Q, Q)
                push!(sys_V, V)
            end
        end
        
        # 找交点
        if !isempty(sys_Q)
            for j in 1:length(sys_V)
                load_Q = calculate_zip_reactive_power(zip_load_mixed, sys_V[j])
                if abs(sys_Q[j] - load_Q) < 0.05  # 容差
                    push!(intersections_Q, sys_Q[j])
                    push!(intersections_V, sys_V[j])
                    break
                end
            end
        end
    end
    
    # 标记交点
    if !isempty(intersections_Q)
        scatter!(p, intersections_Q, intersections_V,
                markersize=6,
                color=:red,
                markershape=:circle,
                label="运行点")
    end
    
    # 添加参考线
    hline!(p, [0.95, 1.05], 
           label="电压限制", 
           color=:gray, 
           linestyle=:dot)
    
    return p
end

"""
ZIP负荷对电压稳定性的影响分析
"""
function zip_load_impact_analysis()
    # 系统参数
    E = 1.05
    X = 0.6
    P_load = 0.6
    
    # 不同ZIP参数组合
    zip_cases = [
        (ZIPLoad(0.4, 1.0, 1.0, 0.0, 0.0), "100% Z"),
        (ZIPLoad(0.4, 1.0, 0.5, 0.5, 0.0), "50% Z + 50% I"),
        (ZIPLoad(0.4, 1.0, 0.0, 1.0, 0.0), "100% I"),
        (ZIPLoad(0.4, 1.0, 0.0, 0.5, 0.5), "50% I + 50% P"),
        (ZIPLoad(0.4, 1.0, 0.0, 0.0, 1.0), "100% P"),
    ]
    
    p = plot(title="ZIP负荷特性对电压稳定性的影响",
             xlabel="无功功率 Q (标幺值)",
             ylabel="电压 V (标幺值)",
             legend=:bottomright)
    
    colors = [:blue, :green, :orange, :purple, :red]
    
    # 绘制系统Q-V特性
    Q_sys = Float64[]
    V_sys = Float64[]
    
    for V in 0.3:0.001:1.1
        sin_delta = P_load * X / (E * V)
        if abs(sin_delta) <= 1.0
            delta = asin(sin_delta)
            cos_delta = cos(delta)
            Q = (E^2 / X) - (E * V / X) * cos_delta
            push!(Q_sys, Q)
            push!(V_sys, V)
        end
    end
    
    plot!(p, Q_sys, V_sys,
          label="系统供电能力",
          linewidth=3,
          color=:black)
    
    # 绘制不同ZIP负荷特性
    for (i, (load, label)) in enumerate(zip_cases)
        Q_load = Float64[]
        V_load = Float64[]
        
        for V in 0.3:0.01:1.1
            Q = calculate_zip_reactive_power(load, V)
            push!(Q_load, Q)
            push!(V_load, V)
        end
        
        plot!(p, Q_load, V_load,
              label=label,
              linewidth=2,
              color=colors[i],
              linestyle=:dash)
    end
    
    # 标记电压稳定边界
    if !isempty(V_sys)
        min_V_idx = argmin(V_sys)
        critical_Q = Q_sys[min_V_idx]
        critical_V = V_sys[min_V_idx]
        
        scatter!(p, [critical_Q], [critical_V],
                markersize=8,
                color=:red,
                markershape=:star,
                label="电压崩溃点")
        
        vline!(p, [critical_Q], 
               label="临界无功功率",
               color=:red,
               linestyle=:dot)
    end
    
    return p
end

"""
电压稳定裕度分析
"""
function voltage_stability_margin_analysis()
    # 系统参数
    E = 1.05
    X_values = [0.4, 0.6, 0.8, 1.0]
    P_load = 0.6
    
    p1 = plot(title="电压稳定裕度 vs 系统电抗",
              xlabel="系统电抗 X (标幺值)",
              ylabel="临界无功功率 Q_critical (标幺值)")
    
    p2 = plot(title="临界电压 vs 系统电抗",
              xlabel="系统电抗 X (标幺值)",
              ylabel="临界电压 V_critical (标幺值)")
    
    Q_criticals = Float64[]
    V_criticals = Float64[]
    
    for X in X_values
        # 计算临界点
        Q_sys = Float64[]
        V_sys = Float64[]
        
        for V in 0.2:0.001:1.1
            sin_delta = P_load * X / (E * V)
            if abs(sin_delta) <= 1.0
                delta = asin(sin_delta)
                cos_delta = cos(delta)
                Q = (E^2 / X) - (E * V / X) * cos_delta
                push!(Q_sys, Q)
                push!(V_sys, V)
            end
        end
        
        if !isempty(V_sys)
            min_idx = argmin(V_sys)
            push!(Q_criticals, Q_sys[min_idx])
            push!(V_criticals, V_sys[min_idx])
        end
    end
    
    plot!(p1, X_values, Q_criticals,
          linewidth=2,
          marker=:circle,
          color=:blue,
          label="Q_critical")
    
    plot!(p2, X_values, V_criticals,
          linewidth=2,
          marker=:circle,
          color=:red,
          label="V_critical")
    
    # 添加安全运行区域
    hline!(p2, [0.95], 
           label="电压下限",
           color=:orange,
           linestyle=:dash)
    
    return plot(p1, p2, layout=(1,2), size=(1000,400))
end

"""
更新主要的电压稳定性分析函数
"""
function analyze_reactive_power_stability()
    # 使用新的Q-V分析函数
    p1 = voltage_stability_analysis()
    p2 = detailed_qv_nose_analysis()
    p3 = zip_load_impact_analysis()
    p4 = voltage_stability_margin_analysis()
    
    # 组合显示
    combined1 = plot(p1, p2, layout=(1,2), size=(1200,500))
    combined2 = plot(p3, p4, layout=(2,1), size=(1000,800))
    
    return [combined1, combined2]
end

"""
计算传输线的无功功率损耗和充电功率
"""
function calculate_line_reactive_power(; 
    length_km=100, 
    voltage_kv=220, 
    current_a=1000,
    x_ohm_per_km=0.4, 
    b_s_per_km=3e-6,
    power_mw=100)
    
    # 线路参数
    X_total = x_ohm_per_km * length_km  # 总电抗
    B_total = b_s_per_km * length_km    # 总电纳
    
    # 电压和电流
    V = voltage_kv * 1000  # 转换为伏特
    I = current_a
    
    # 无功功率损耗 (感性)
    Q_loss = I^2 * X_total / 1e6  # MVAr
    
    # 充电功率 (容性)
    Q_charging = V^2 * B_total / 1e6  # MVAr
    
    # 净无功功率需求
    Q_net = Q_loss - Q_charging
    
    return (
        Q_loss = Q_loss,
        Q_charging = Q_charging, 
        Q_net = Q_net,
        X_total = X_total,
        B_total = B_total
    )
end

"""
绘制无功功率随距离变化的曲线
"""
function plot_reactive_power_vs_distance()
    distances = 0:10:500  # 0-500km
    Q_losses = Float64[]
    Q_chargings = Float64[]
    Q_nets = Float64[]
    
    for d in distances
        result = calculate_line_reactive_power(length_km=d)
        push!(Q_losses, result.Q_loss)
        push!(Q_chargings, result.Q_charging)
        push!(Q_nets, result.Q_net)
    end
    
    p = plot(distances, Q_losses, 
             label="无功损耗 (感性)", 
             linewidth=2, 
             color=:red,
             title="无功功率随传输距离的变化",
             xlabel="距离 (km)",
             ylabel="无功功率 (MVAr)")
    
    plot!(p, distances, Q_chargings, 
          label="充电功率 (容性)", 
          linewidth=2, 
          color=:blue)
    
    plot!(p, distances, Q_nets, 
          label="净无功需求", 
          linewidth=2, 
          color=:green,
          linestyle=:dash)
    
    # 添加临界点
    critical_idx = findfirst(x -> x <= 0, Q_nets)
    if critical_idx !== nothing
        critical_distance = distances[critical_idx]
        vline!(p, [critical_distance], 
               label="临界距离: $(critical_distance)km", 
               color=:black, 
               linestyle=:dot)
    end
    
    return p
end

"""
电压降随无功功率传输的变化
"""
function plot_voltage_drop_vs_reactive_power()
    Q_values = 0:5:200  # 0-200 MVAr
    voltage_drops = Float64[]
    
    # 系统参数
    V_base = 220e3  # 基准电压 220kV
    R_line = 0.05   # 线路电阻 Ω/km
    X_line = 0.4    # 线路电抗 Ω/km
    length = 200    # 线路长度 km
    P_constant = 100e6  # 恒定有功功率 100MW
    
    R_total = R_line * length
    X_total = X_line * length
    
    for Q in Q_values
        Q_watts = Q * 1e6  # 转换为VAr
        
        # 计算电流幅值
        S = sqrt(P_constant^2 + Q_watts^2)
        I = S / (sqrt(3) * V_base)
        
        # 正确的电压降计算公式: ΔV = (P*R + Q*X)/V
        ΔV_actual = (P_constant * R_total + Q_watts * X_total) / (sqrt(3) * V_base)
        ΔV_percent = (ΔV_actual / V_base) * 100
        
        push!(voltage_drops, ΔV_percent)
    end
    @show voltage_drops  # 显示电压降数据
    @show Q_values  # 显示无功功率数据
    p = plot(Q_values, voltage_drops,
             title="电压降随无功功率传输的变化 (200km线路)",
             xlabel="传输无功功率 (MVAr)",
             ylabel="电压降 (%)",
             linewidth=2,
             color=:purple,
             label="电压降 = (P·R + Q·X)/V")
    
    # 添加5%电压降限制线
    hline!(p, [5.0], 
           label="5%电压降限制", 
           color=:red, 
           linestyle=:dash)
    
    # 添加10%电压降限制线
    hline!(p, [10.0], 
           label="10%电压降限制", 
           color=:orange, 
           linestyle=:dash)
    
    return p
end

"""
比较不同电压等级的无功传输能力
"""
function compare_voltage_levels()
    voltage_levels = [110, 220, 500, 750]  # kV
    distances = 0:20:400
    
    # 标准化参数 - 每个电压等级的典型线路参数
    line_params = Dict(
        110 => (x_ohm_km=0.4, b_s_km=2.8e-6, p_mw=50),    # 110kV线路
        220 => (x_ohm_km=0.4, b_s_km=3.0e-6, p_mw=200),   # 220kV线路  
        500 => (x_ohm_km=0.3, b_s_km=3.5e-6, p_mw=1000),  # 500kV线路
        750 => (x_ohm_km=0.25, b_s_km=4.0e-6, p_mw=2000)  # 750kV线路
    )
    
    p = plot(title="不同电压等级的无功传输特性",
             xlabel="距离 (km)",
             ylabel="无功功率损耗 (MVAr)")
    
    colors = [:red, :blue, :green, :orange]
    
    for (i, V) in enumerate(voltage_levels)
        Q_losses = Float64[]
        params = line_params[V]
        
        for d in distances
            # 根据电压等级计算实际电流
            P_watts = params.p_mw * 1e6
            I_actual = P_watts / (sqrt(3) * V * 1e3 * 0.9)  # 假设功率因数0.9
            
            result = calculate_line_reactive_power(
                length_km=d, 
                voltage_kv=V, 
                current_a=I_actual,
                x_ohm_per_km=params.x_ohm_km,
                b_s_per_km=params.b_s_km,
                power_mw=params.p_mw
            )
            push!(Q_losses, result.Q_loss)
        end
        
        plot!(p, distances, Q_losses,
              label="$(V) kV ($(line_params[V].p_mw) MW)",
              linewidth=2,
              color=colors[i])
    end
    
    return p
end

"""
无功功率的频率响应特性
"""
function plot_frequency_response()
    frequencies = 45:0.1:55  # 45-55 Hz
    f0 = 50  # 基准频率
    
    X_ratios = frequencies / f0      # 电抗比例
    B_ratios = frequencies / f0      # 电纳比例
    
    p = plot(frequencies, X_ratios,
             label="电抗 X ∝ f",
             linewidth=2,
             color=:red,
             title="无功特性的频率响应",
             xlabel="频率 (Hz)",
             ylabel="标幺值")
    
    plot!(p, frequencies, B_ratios,
          label="电纳 B ∝ f", 
          linewidth=2,
          color=:blue)
    
    # 标记额定频率
    vline!(p, [50], 
           label="额定频率 50Hz", 
           color=:black, 
           linestyle=:dot)
    
    return p
end

"""
负荷功率因数对无功需求的影响
"""
function plot_power_factor_impact()
    power_factors = 0.7:0.01:1.0
    P_load = 100  # MW
    
    Q_demands = Float64[]
    
    for pf in power_factors
        Q = P_load * tan(acos(pf))  # 无功功率需求
        push!(Q_demands, Q)
    end
    
    p = plot(power_factors, Q_demands,
             title="功率因数对无功需求的影响",
             xlabel="功率因数",
             ylabel="无功功率需求 (MVAr)",
             linewidth=2,
             color=:green,
             label="Q = P·tan(φ)")
    
    # 标记典型值
    vline!(p, [0.9], 
           label="典型工业负荷 (cos φ = 0.9)", 
           color=:red, 
           linestyle=:dash)
    
    return p
end

"""
演示无功补偿的效果
"""
function demonstrate_reactive_compensation()
    distances = 0:25:300
    
    # 无补偿情况
    Q_losses_uncomp = Float64[]
    voltage_drops_uncomp = Float64[]
    
    # 有补偿情况 (中点补偿)
    Q_losses_comp = Float64[]
    voltage_drops_comp = Float64[]
    
    # 线路参数
    R_ohm_km = 0.05
    X_ohm_km = 0.4
    P_mw = 200
    V_kv = 220
    
    for d in distances
        # 计算电流
        P_watts = P_mw * 1e6
        I = P_watts / (sqrt(3) * V_kv * 1e3 * 0.9)
        
        # 无补偿
        result_uncomp = calculate_line_reactive_power(
            length_km=d, 
            voltage_kv=V_kv,
            current_a=I,
            x_ohm_per_km=X_ohm_km
        )
        push!(Q_losses_uncomp, result_uncomp.Q_loss)
        
        # 电压降计算
        R_total = R_ohm_km * d
        X_total = X_ohm_km * d
        Q_net = result_uncomp.Q_net * 1e6  # 转换为VAr
        V_drop_uncomp = (P_watts * R_total + abs(Q_net) * X_total) / (sqrt(3) * V_kv * 1e3)
        V_drop_percent_uncomp = (V_drop_uncomp / (V_kv * 1e3)) * 100
        push!(voltage_drops_uncomp, V_drop_percent_uncomp)
        
        # 有补偿 (中点50%补偿)
        Q_comp_loss = result_uncomp.Q_loss * 0.25  # 补偿后损耗减少75%
        push!(Q_losses_comp, Q_comp_loss)
        
        Q_comp_net = Q_comp_loss * 1e6
        V_drop_comp = (P_watts * R_total + Q_comp_net * X_total) / (sqrt(3) * V_kv * 1e3)
        V_drop_percent_comp = (V_drop_comp / (V_kv * 1e3)) * 100
        push!(voltage_drops_comp, V_drop_percent_comp)
    end
    
    p1 = plot(distances, Q_losses_uncomp,
              label="无补偿",
              linewidth=2,
              color=:red,
              title="无功补偿效果对比",
              xlabel="距离 (km)",
              ylabel="无功损耗 (MVAr)")
    
    plot!(p1, distances, Q_losses_comp,
          label="中点补偿",
          linewidth=2,
          color=:blue)
    
    p2 = plot(distances, voltage_drops_uncomp,
              label="无补偿",
              linewidth=2,
              color=:red,
              title="电压降对比",
              xlabel="距离 (km)",
              ylabel="电压降 (%)")
    
    plot!(p2, distances, voltage_drops_comp,
          label="中点补偿",
          linewidth=2,
          color=:blue)
    
    # 添加电压降限制线
    hline!(p2, [5.0], 
           label="5%限制", 
           color=:orange, 
           linestyle=:dash)
    
    return plot(p1, p2, layout=(2,1), size=(800,600))
end

"""
添加新函数：分析无功功率对系统稳定性的影响
"""
function analyze_reactive_power_stability()
    # 使用新的详细分析函数
    p1 = voltage_stability_analysis()
    p2 = detailed_qv_nose_analysis()
    p3 = zip_load_impact_analysis()
    p4 = voltage_stability_margin_analysis()
    
    # 组合显示
    combined1 = plot(p1, p2, layout=(1,2), size=(1200,500))
    combined2 = plot(p3, p4, layout=(2,1), size=(1000,800))
    
    return [combined1, combined2]
end

"""
运行所有无功功率传输示例
"""
function run_reactive_power_examples()
    println("=== 无功功率远距离传输限制演示 ===\n")
    
    # 1. 基本计算示例
    println("1. 基本计算示例:")
    result = calculate_line_reactive_power(length_km=200)
    println("   线路长度: 200 km")
    println("   无功损耗: $(round(result.Q_loss, digits=2)) MVAr")
    println("   充电功率: $(round(result.Q_charging, digits=2)) MVAr") 
    println("   净无功需求: $(round(result.Q_net, digits=2)) MVAr\n")
    
    # 2. 绘制图形
    println("2. 生成分析图形...")
    
    plots_to_show = [
        ("无功功率vs距离", plot_reactive_power_vs_distance()),
        ("电压降vs无功传输", plot_voltage_drop_vs_reactive_power()),
        ("不同电压等级对比", compare_voltage_levels()),
        ("频率响应特性", plot_frequency_response()),
        ("功率因数影响", plot_power_factor_impact()),
        ("无功补偿效果", demonstrate_reactive_compensation()),
        ("系统稳定性分析", analyze_reactive_power_stability())
    ]
    
    return plots_to_show
end
