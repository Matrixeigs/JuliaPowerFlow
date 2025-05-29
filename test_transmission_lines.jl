include("./reactive_power_examples.jl")

"""
主函数：运行所有无功功率传输理论分析
根据 reactive_power_transmission_theory.md 文档实现完整分析
"""
function run_comprehensive_reactive_power_analysis()
    println("="^70)
    println("无功功率远距离传输限制 - 完整理论分析")
    println("基于 reactive_power_transmission_theory.md 文档")
    println("="^70)
    
    # 1. 基础理论验证
    println("\n1. 传输线基础理论验证")
    println("-"^30)
    
    # 创建典型500kV线路
    line_500kv = TransmissionLine(0.025, 0.8e-3, 14e-9, 0.0, 200.0, 50.0)
    chars = calculate_line_characteristics(line_500kv)
    
    println("500kV线路特征参数:")
    println("   特征阻抗 Z₀: $(round(abs(chars.Z₀), digits=1)) Ω")
    println("   传播常数 γ: $(round(abs(chars.γ)*1000, digits=3)) × 10⁻³ /km")
    println("   自然功率 SIL: $(round(calculate_surge_impedance_loading(500.0, chars.Z₀), digits=1)) MW")
    
    # 2. 分布参数分析
    println("\n2. 分布参数功率传输分析")
    println("-"^30)
    
    result = distributed_parameter_analysis(line_500kv, 500e3, 480e3, π/6)
    println("发送端功率: P₁=$(round(result.P₁, digits=1))MW, Q₁=$(round(result.Q₁, digits=1))MVAr")
    println("接收端功率: P₂=$(round(result.P₂, digits=1))MW, Q₂=$(round(result.Q₂, digits=1))MVAr") 
    println("线路损耗: ΔP=$(round(result.P_loss, digits=1))MW, ΔQ=$(round(result.Q_loss, digits=1))MVAr")
    
    # 3. 临界长度分析
    println("\n3. 临界长度分析")
    println("-"^30)
    
    critical_results = analyze_critical_length()
    for (V, data) in critical_results
        if !isnan(data.critical_length)
            println("$(V)kV线路临界长度: $(data.critical_length)km (SIL: $(round(data.SIL, digits=1))MW)")
        else
            println("$(V)kV线路未达到临界长度 (SIL: $(round(data.SIL, digits=1))MW)")
        end
    end
    
    # 4. 生成所有分析图形
    println("\n4. 生成理论分析图形")
    println("-"^30)
    
    analysis_plots = [
        ("基础距离特性", plot_reactive_power_vs_distance()),
        ("电压降分析", plot_voltage_drop_vs_reactive_power()),
        ("临界长度分析", plot_critical_length_analysis()),
        ("电压稳定性分析", voltage_stability_analysis()),
        ("频率影响分析", frequency_impact_analysis()),
        ("无功功率分布", reactive_power_distribution_analysis()),
        ("温度效应分析", temperature_effect_analysis()),
        ("补偿策略对比", compensation_strategy_analysis()),
        ("不同电压等级对比", compare_voltage_levels()),
        ("功率因数影响", plot_power_factor_impact()),
        ("无功补偿效果", demonstrate_reactive_compensation())
    ]
    
    # 添加新的电压稳定性分析图形 - 使用push!而不是append!
    stability_plots = analyze_reactive_power_stability()
    push!(analysis_plots, ("ZIP负荷电压稳定性分析1", stability_plots[1]))
    push!(analysis_plots, ("ZIP负荷电压稳定性分析2", stability_plots[2]))
    
    # 显示所有图形
    for (title, plot_obj) in analysis_plots
        println("\n生成图形: $title")
        display(plot_obj)
    end
    
    # 5. 理论总结
    println("\n" * "="^70)
    println("理论分析总结")
    println("="^70)
    
    println("\n📊 数学模型验证:")
    println("✓ 分布参数传输线模型: γ = √[(R+jωL)(G+jωC)]")
    println("✓ 特征阻抗计算: Z₀ = √[(R+jωL)/(G+jωC)]")
    println("✓ 功率传输方程: P = (V₁V₂/X)sin(δ), Q = (V₁²/X) - (V₁V₂/X)cos(δ)")
    
    println("\n🔬 物理机理解释:")
    println("✓ 无功损耗 ∝ I²X: 电流平方与电抗的乘积")
    println("✓ 充电功率 ∝ V²B: 电压平方与电纳的乘积")
    println("✓ 临界长度: 充电功率 = 自然功率时的距离")
    
    println("\n📈 工程应用指导:")
    println("✓ 高电压等级传输效率更高: P ∝ V², I ∝ 1/V")
    println("✓ 无功就地平衡原则: 避免长距离传输")
    println("✓ 补偿设备布置: 串联补偿+并联补偿组合")
    println("✓ 电压稳定性: 维持足够的无功储备")
    
    println("\n🎯 关键结论:")
    println("1. 无功功率不能远距离传输的根本原因是电抗损耗")
    println("2. 存在物理极限 - 临界长度概念")
    println("3. 频率和温度变化直接影响无功特性")
    println("4. 综合补偿策略是提高传输能力的关键")
    println("5. 电压稳定性与无功功率密切相关")
    
    return analysis_plots
end

# === 主程序执行 ===
# if abspath(PROGRAM_FILE) == @__FILE__
    # 运行完整的无功功率传输分析
    analysis_results = run_comprehensive_reactive_power_analysis()
    
    println("\n" * "="^70)
    println("分析完成！所有图形已生成并显示。")
    println("详细理论推导请参考: reactive_power_transmission_theory.md")
    println("="^70)
# end