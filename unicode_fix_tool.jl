# Unicode乱码检测与修复脚本
# Unicode Character Detection and Fix Script

"""
检测并修复Markdown文件中的Unicode矢量符号乱码问题
Detect and fix Unicode vector symbol encoding issues in Markdown files
"""

using Unicode

function detect_unicode_issues(file_path::String)
    println("="^70)
    println("Unicode乱码检测报告 / Unicode Issue Detection Report")
    println("文件路径 / File Path: $file_path")
    println("="^70)
    
    # 读取文件内容
    content = read(file_path, String)
    lines = split(content, '\n')
    
    # 问题符号列表
    problematic_chars = [
        "⃗",  # 矢量上标箭头
        "→",  # 右箭头
        "⃗",  # 组合箭头
        "⃗",  # 另一种矢量箭头
    ]
    
    # 常见的矢量符号组合
    vector_patterns = [
        r"[EFBHJD]⃗",  # 矢量符号：E⃗, F⃗, B⃗, H⃗, J⃗, D⃗
        r"dl⃗",        # 线元素
        r"dA⃗",        # 面元素
        r"r̂",         # 单位矢量
    ]
    
    issues_found = 0
    
    for (line_num, line) in enumerate(lines)
        # 检测问题字符
        for char in problematic_chars
            if contains(line, char)
                println("Line $line_num: 发现矢量箭头符号 '$char'")
                println("  内容: $(strip(line))")
                issues_found += 1
            end
        end
        
        # 检测矢量模式
        for pattern in vector_patterns
            matches = eachmatch(pattern, line)
            if !isempty(collect(matches))
                for m in matches
                    println("Line $line_num: 发现矢量符号模式 '$(m.match)'")
                    println("  内容: $(strip(line))")
                    issues_found += 1
                end
            end
        end
    end
    
    println("\n总结 / Summary:")
    if issues_found > 0
        println("⚠️  发现 $issues_found 处Unicode矢量符号问题")
        println("   建议进行修复以确保最佳兼容性")
    else
        println("✅ 未发现Unicode乱码问题")
        println("   文档符号兼容性良好")
    end
    
    return issues_found
end

function fix_unicode_vectors(input_file::String, output_file::String)
    println("\n" * "="^70)
    println("Unicode矢量符号修复 / Unicode Vector Symbol Fix")
    println("="^70)
    
    # 读取原文件
    content = read(input_file, String)
    
    # 矢量符号替换映射
    replacements = [
        # 矢量符号修复
        "E⃗" => "E",
        "F⃗" => "F", 
        "B⃗" => "B",
        "H⃗" => "H",
        "J⃗" => "J",
        "D⃗" => "D",
        "S⃗" => "S",
        "g⃗" => "g",
        "l⃗" => "l",
        "dl⃗" => "dl",
        "dA⃗" => "dA",
        "r̂" => "r̂",  # 保留单位矢量符号，但确保兼容性
        
        # 其他可能的矢量符号
        "∇⃗" => "∇",
        "∮⃗" => "∮",
        "∫∫⃗" => "∫∫",
        
        # 特殊情况处理
        "J⃗_D" => "J_D",
        "J⃗_free" => "J_free",
    ]
    
    fixed_content = content
    fixes_applied = 0
    
    for (old_symbol, new_symbol) in replacements
        old_count = length(collect(eachmatch(Regex(old_symbol), fixed_content)))
        if old_count > 0
            fixed_content = replace(fixed_content, old_symbol => new_symbol)
            println("✅ 修复: '$old_symbol' → '$new_symbol' ($old_count 处)")
            fixes_applied += old_count
        end
    end
    
    # 写入修复后的文件
    write(output_file, fixed_content)
    
    println("\n修复完成 / Fix Complete:")
    println("📁 输入文件: $input_file")
    println("📁 输出文件: $output_file")
    println("🔧 总计修复: $fixes_applied 处矢量符号")
    
    return fixes_applied
end

function generate_compatibility_report()
    println("\n" * "="^70)
    println("兼容性建议报告 / Compatibility Recommendation Report")
    println("="^70)
    
    println("📋 Unicode矢量符号兼容性分析:")
    println()
    println("1. 问题原因 / Root Causes:")
    println("   • Unicode组合字符在某些环境下显示异常")
    println("   • 不同操作系统和软件的字体支持差异")
    println("   • 终端和编辑器的Unicode渲染能力限制")
    println()
    println("2. 修复策略 / Fix Strategy:")
    println("   • 移除矢量上标箭头符号 (⃗)")
    println("   • 保持数学意义不变")
    println("   • 使用标准ASCII字符替代")
    println("   • 确保跨平台兼容性")
    println()
    println("3. 建议做法 / Recommended Practices:")
    println("   • 使用E而非E⃗表示电场矢量")
    println("   • 在注释中说明符号为矢量")
    println("   • 采用粗体或其他格式强调矢量性质")
    println("   • 优先考虑可读性和兼容性")
    println()
    println("4. 验证方法 / Verification Methods:")
    println("   • 在多种环境下测试显示效果")
    println("   • 检查终端输出和PDF生成")
    println("   • 确保代码和文档的一致性")
end

# 主程序执行
function main()
    # 文件路径
    input_file = "docs/electromagnetic_induction_theory.md"
    output_file = "docs/electromagnetic_induction_theory_fixed.md"
    
    # 1. 检测当前文件的Unicode问题
    issues = detect_unicode_issues(input_file)
    
    # 2. 如果发现问题，进行修复
    if issues > 0
        fixes = fix_unicode_vectors(input_file, output_file)
        
        # 3. 验证修复后的文件
        println("\n验证修复结果 / Verifying Fix Results:")
        remaining_issues = detect_unicode_issues(output_file)
        
        if remaining_issues == 0
            println("🎉 修复成功！文件已完全兼容")
        else
            println("⚠️  仍有 $remaining_issues 处问题需要手动检查")
        end
    end
    
    # 4. 生成兼容性报告
    generate_compatibility_report()
    
    println("\n" * "="^70)
    println("Unicode修复程序执行完成 / Unicode Fix Program Completed")
    println("="^70)
end

# 如果直接运行此脚本
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
