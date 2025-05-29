include("./transformer_tap_control.jl")

# 运行变压器分接头控制分析
println("开始三绕组变压器分接头电压控制分析...")
results = run_transformer_tap_analysis()

println("\n分析完成！")
println("详细理论说明请参考: docs/transformer_tap_voltage_control.md")
