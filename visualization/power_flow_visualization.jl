using LinearAlgebra, SparseArrays, Printf
using Plots, ColorSchemes
using Statistics
include("../data/case9.jl")

"""
    构建节点导纳矩阵
"""
function build_ybus(jpc)
    # 获取基础数据
    baseMVA = jpc["baseMVA"]
    bus = jpc["bus"]
    branch = jpc["branch"]
    
    # 获取节点和支路数量
    nb = size(bus, 1)    # 节点数
    nl = size(branch, 1) # 支路数
    
    # 初始化Ybus为稀疏矩阵
    Ybus = spzeros(Complex{Float64}, nb, nb)
    
    # 处理支路数据
    for i in 1:nl
        # 获取支路起始和终止节点
        from_bus = Int(branch[i, 1])
        to_bus = Int(branch[i, 2])
        
        # 获取支路参数
        r = branch[i, 3]
        x = branch[i, 4]
        b = branch[i, 5]
        
        # 计算支路导纳
        y = 1.0 / complex(r, x)
        
        # 计算分路导纳
        b_sh = complex(0, b)
        
        # 获取变压器变比
        tap = branch[i, 10]
        if tap == 0
            tap = 1.0
        end
        
        # 获取相移角度（弧度）
        shift = branch[i, 11] * (π/180)
        
        # 计算变压器复变比
        tap_complex = tap * exp(im * shift)
        
        # 填充Ybus矩阵
        if tap == 1.0 && shift == 0  # 普通线路
            Ybus[from_bus, from_bus] += y + b_sh/2
            Ybus[to_bus, to_bus] += y + b_sh/2
            Ybus[from_bus, to_bus] -= y
            Ybus[to_bus, from_bus] -= y
        else  # 变压器
            Ybus[from_bus, from_bus] += y / (tap_complex * conj(tap_complex))
            Ybus[to_bus, to_bus] += y
            Ybus[from_bus, to_bus] -= y / conj(tap_complex)
            Ybus[to_bus, from_bus] -= y / tap_complex
        end
    end
    
    # 添加节点分路导纳
    for i in 1:nb
        # 获取节点导纳
        gs = bus[i, 5]  # 电导
        bs = bus[i, 6]  # 电纳
        
        if gs != 0 || bs != 0
            Ybus[i, i] += complex(gs, bs)
        end
    end
    
    return Ybus
end

"""
    初始化电压向量
"""
function initialize_voltage(jpc)
    bus = jpc["bus"]
    gen = jpc["gen"]
    nb = size(bus, 1)
    
    # 初始化电压向量为标称值
    V = ones(Complex{Float64}, nb)
    
    # 设置发电机节点的电压幅值
    for i in axes(gen, 1)
        gen_bus = Int(gen[i, 1])
        V[gen_bus] = gen[i, 6] * exp(im * 0.0)  # 使用发电机设定电压幅值
    end
    
    return V
end

"""
    确定节点类型
"""
function identify_bus_types(jpc)
    bus = jpc["bus"]
    nb = size(bus, 1)
    
    # 初始化节点类型数组
    pv = Int[]  # PV节点
    pq = Int[]  # PQ节点
    ref = Int[] # 参考节点
    
    for i in 1:nb
        bus_type = Int(bus[i, 2])
        if bus_type == 3  # 参考节点
            push!(ref, i)
        elseif bus_type == 2  # PV节点
            push!(pv, i)
        elseif bus_type == 1  # PQ节点
            push!(pq, i)
        end
    end
    
    return ref, pv, pq
end

"""
    计算节点功率注入
"""
function calculate_power_injection(jpc)
    bus = jpc["bus"]
    gen = jpc["gen"]
    baseMVA = jpc["baseMVA"]
    nb = size(bus, 1)
    
    # 初始化功率注入向量
    Sbus = zeros(Complex{Float64}, nb)
    
    # 添加负荷功率
    for i in 1:nb
        pd = bus[i, 3] / baseMVA  # 有功负荷（标幺值）
        qd = bus[i, 4] / baseMVA  # 无功负荷（标幺值）
        Sbus[i] = -complex(pd, qd)  # 负荷为负注入
    end
    
    # 添加发电机功率
    for i in 1:size(gen, 1)
        gen_bus = Int(gen[i, 1])
        pg = gen[i, 2] / baseMVA  # 有功发电（标幺值）
        qg = gen[i, 3] / baseMVA  # 无功发电（标幺值）
        Sbus[gen_bus] += complex(pg, qg)  # 发电为正注入
    end
    
    return Sbus
end

"""
    牛顿-拉夫逊潮流计算（带可视化）
"""
function newton_raphson_power_flow_with_visualization(Ybus, Sbus, V0, pv, pq, ref; 
                                  max_iter=30, tol=1e-6, verbose=true)
    # 初始化
    n = length(V0)           # 节点总数
    V = copy(V0)             # 电压向量
    Va = angle.(V)           # 相角
    Vm = abs.(V)             # 电压幅值
    iteration = 0            # 迭代计数
    converged = false        # 收敛标志
    
    # 设置参考节点的电压和相角
    Va[ref] .= angle.(V0[ref])
    Vm[ref] .= abs.(V0[ref])
    
    # 设置PV节点的电压幅值
    Vm[pv] .= abs.(V0[pv])
    
    # 需要更新的变量索引
    pvpq = [pv; pq]          # PV和PQ节点的相角
    
    # 初始化用于可视化的数组
    iter_history = Int[]
    mismatches_history = Float64[]
    jacobian_plots = []
    voltage_magnitude_plots = []
    voltage_angle_plots = []
    
    if verbose
        println("牛顿-拉夫逊潮流计算开始")
        println("最大迭代次数: $max_iter")
        println("收敛容差: $tol")
        println("节点总数: $n")
        println("PV节点数: $(length(pv))")
        println("PQ节点数: $(length(pq))")
        println("参考节点: $ref")
    end
    
    # 开始迭代
    while !converged && iteration < max_iter
        iteration += 1
        push!(iter_history, iteration)
        
        # 更新电压向量
        V = Vm .* exp.(1im .* Va)
        
        # 计算节点功率注入
        Ibus = Ybus * V
        S = V .* conj(Ibus)
        P = real.(S)
        Q = imag.(S)
        
        # 计算功率不平衡
        dP = real.(Sbus) - P
        dQ = imag.(Sbus) - Q
        
        # 检查收敛性
        F = [dP[pvpq]; dQ[pq]]
        norm_F = norm(F, Inf)
        push!(mismatches_history, norm_F)
        
        if verbose
            @printf("迭代 %3d: 最大功率不平衡 = %.10f\n", iteration, norm_F)
        end
        
        # 保存当前迭代的电压幅值和相角
        vm_plot = bar(1:n, Vm, 
                      title="Voltage Magnitude ($iteration)", 
                      xlabel="Bus", 
                      ylabel="Voltage Magnitude (p.u.)",
                      legend=false,
                      ylims=(0.9, 1.1),
                      size=(600, 400))
        push!(voltage_magnitude_plots, vm_plot)
        
        va_plot = bar(1:n, rad2deg.(Va), 
                      title="Voltage Angle ( $iteration)", 
                      xlabel="Bus", 
                      ylabel="Voltage Angle (degrees)",
                      legend=false,
                      size=(600, 400))
        push!(voltage_angle_plots, va_plot)
        
        if norm_F < tol
            converged = true
            if verbose
                println("牛顿-拉夫逊潮流计算已收敛!")
            end
            break
        end
        
        # 构造雅可比矩阵
        J = create_jacobian(Ybus, V, pv, pq)
        
        # 绘制雅可比矩阵热力图
        j_plot = plot_jacobian_heatmap(J, iteration)
        push!(jacobian_plots, j_plot)
        
        # 求解线性方程组
        dx = J \ F
        
        # 更新状态变量
        # 更新相角 (PV和PQ节点)
        Va[pvpq] = Va[pvpq] + dx[1:length(pvpq)]
        
        # 更新电压幅值 (仅PQ节点)
        if length(pq) > 0
            Vm[pq] = Vm[pq] .* (1.0 .+ dx[length(pvpq)+1:end])
        end
    end
    
    if !converged
        @warn "牛顿-拉夫逊潮流计算未在最大迭代次数内收敛!"
    end
    
    # 最终更新电压向量
    V = Vm .* exp.(1im .* Va)
    
    # 计算最终的功率注入
    Ibus = Ybus * V
    S = V .* conj(Ibus)
    
    # 绘制不平衡量变化图
    mismatch_plot = plot(iter_history, mismatches_history, 
                        title="Power Mismatch Convergence", 
                        xlabel="Iterations", 
                        ylabel="maximum mismatch (p.u.)",
                        marker=:circle,
                        yscale=:log10,
                        legend=false,
                        size=(800, 500))
    
    # 创建动画
    jacobian_anim = @animate for p in jacobian_plots
        plot(p)
    end
    
    voltage_mag_anim = @animate for p in voltage_magnitude_plots
        plot(p)
    end
    
    voltage_ang_anim = @animate for p in voltage_angle_plots
        plot(p)
    end
    
    # 保存动画
    gif(jacobian_anim, "jacobian_heatmap.gif", fps = 1)
    gif(voltage_mag_anim, "voltage_magnitude.gif", fps = 1)
    gif(voltage_ang_anim, "voltage_angle.gif", fps = 1)
    
    return V, converged, iteration, S, mismatch_plot
end

"""
    绘制雅可比矩阵热力图
"""
function plot_jacobian_heatmap(J, iteration)
    p = heatmap(J, 
                color=:viridis, 
                aspect_ratio=:equal,
                title="Heat Map of Jacobian Matrix ($iteration)",
                xlabel="Variable Index",
                ylabel="Power Injection Index",
                colorbar_title="Jacobian Value",
                clim=(-maximum(abs.(J)), maximum(abs.(J))),
                size=(600, 500))
    return p
end

"""
    主函数：执行潮流计算（带可视化）
"""
function run_power_flow_with_visualization(case_data)
    # 构建节点导纳矩阵
    Ybus = build_ybus(case_data)
    
    # 初始化电压向量
    V0 = initialize_voltage(case_data)
    
    # 确定节点类型
    ref, pv, pq = identify_bus_types(case_data)
    
    # 计算节点功率注入
    Sbus = calculate_power_injection(case_data)
    
    # 执行牛顿-拉夫逊潮流计算（带可视化）
    V, converged, iterations, S, mismatch_plot = newton_raphson_power_flow_with_visualization(
        Ybus, Sbus, V0, pv, pq, ref, verbose=true
    )
    
    # 输出结果
    println("\n最终结果:")
    println("收敛状态: ", converged ? "已收敛" : "未收敛")
    println("迭代次数: ", iterations)
    
    # 输出节点电压
    println("\n节点电压:")
    println("节点\t幅值(p.u.)\t相角(度)")
    println("------------------------------------")
    for i in eachindex(V)
        @printf("%d\t%.6f\t%.6f\n", i, abs(V[i]), rad2deg(angle(V[i])))
    end
    
    # 计算线路功率流
    Sf, St = calculate_line_flows(case_data, V)
    
    # 输出线路功率流
    println("\n线路功率流:")
    println("支路\t从节点\t到节点\t从节点功率(MVA)\t\t到节点功率(MVA)")
    println("--------------------------------------------------------------")
    for i in eachindex(Sf)
        from_bus = Int(case_data["branch"][i, 1])
        to_bus = Int(case_data["branch"][i, 2])
        @printf("%d\t%d\t%d\t%.2f + j%.2f\t%.2f + j%.2f\n", 
                i, from_bus, to_bus, real(Sf[i]), imag(Sf[i]), real(St[i]), imag(St[i]))
    end
    
    # 计算系统损耗
    losses = sum(Sf + St)
    @printf("\n系统总损耗: %.2f + j%.2f MVA\n", real(losses), imag(losses))
    
    # 保存不平衡量变化图
    savefig(mismatch_plot, "mismatch_plot.png")
    
    # 显示不平衡量变化图
    display(mismatch_plot)
    
    # 创建最终电压幅值和相角图
    vm_final = bar(1:length(V), abs.(V), 
                  title="Final Voltage Magnitude", 
                  xlabel="Bus", 
                  ylabel="Voltage Magnitude (p.u.)",
                  legend=false,
                  size=(600, 400))
    
    va_final = bar(1:length(V), rad2deg.(angle.(V)), 
                  title="Final Voltage Angle", 
                  xlabel="Bus", 
                  ylabel="Angle (degrees)",
                  legend=false,
                  size=(600, 400))
    
    # 保存最终电压图
    savefig(vm_final, "voltage_magnitude_final.png")
    savefig(va_final, "voltage_angle_final.png")
    
    # 显示最终电压图
    display(vm_final)
    display(va_final)
    
    return V, S, Sf, St
end