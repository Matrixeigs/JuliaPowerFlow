# transformations.jl - Coordinate Transformations

# ============================================================================
# COORDINATE TRANSFORMATIONS
# ============================================================================

"""
Park Transformation: Convert ABC frame to DQ0 frame

Parameters:
- abc: 3-element vector [a, b, c]
- θ: transformation angle (rad)

Returns:
- dq0: 3-element vector [d, q, 0]
"""
function abc_to_dq0(abc::Vector{Float64}, θ::Float64)
    # Park transformation matrix
    T = (2/3) * [cos(θ)         cos(θ - 2π/3)   cos(θ + 2π/3);
                 -sin(θ)        -sin(θ - 2π/3)  -sin(θ + 2π/3);
                 1/2            1/2             1/2]
    
    return T * abc
end

"""
Inverse Park Transformation: Convert DQ0 frame to ABC frame

Parameters:
- dq0: 3-element vector [d, q, 0]
- θ: transformation angle (rad)

Returns:
- abc: 3-element vector [a, b, c]
"""
function dq0_to_abc(dq0::Vector{Float64}, θ::Float64)
    # Inverse Park transformation matrix
    T_inv = [cos(θ)           -sin(θ)          1;
             cos(θ - 2π/3)    -sin(θ - 2π/3)   1;
             cos(θ + 2π/3)    -sin(θ + 2π/3)   1]
    
    return T_inv * dq0
end

"""
Clarke Transformation: Convert ABC frame to αβ0 frame

Parameters:
- abc: 3-element vector [a, b, c]

Returns:
- αβ0: 3-element vector [α, β, 0]
"""
function abc_to_αβ0(abc::Vector{Float64})
    # Clarke transformation matrix
    T_clarke = (2/3) * [1        -1/2      -1/2;
                        0        √3/2      -√3/2;
                        1/2      1/2       1/2]
    
    return T_clarke * abc
end

"""
Inverse Clarke Transformation: Convert αβ0 frame to ABC frame

Parameters:
- αβ0: 3-element vector [α, β, 0]

Returns:
- abc: 3-element vector [a, b, c]
"""
function αβ0_to_abc(αβ0::Vector{Float64})
    # Inverse Clarke transformation matrix
    T_clarke_inv = [1        0        1;
                    -1/2     √3/2     1;
                    -1/2     -√3/2    1]
    
    return T_clarke_inv * αβ0
end

"""
Sequence Component Transformation: Convert ABC to 012 (positive, negative, zero)

Parameters:
- abc: 3-element vector [a, b, c]

Returns:
- seq: 3-element vector [0, 1, 2] (zero, positive, negative sequence)
"""
function abc_to_sequence(abc::Vector{Float64})
    # Sequence transformation matrix
    a = exp(im * 2π/3)  # Complex operator
    T_seq = (1/3) * [1    1     1;
                     1    a     a^2;
                     1    a^2   a]
    
    # Convert to complex for calculation
    abc_complex = complex.(abc)
    seq_complex = T_seq * abc_complex
    
    # Return real parts (assuming balanced system for simplicity)
    return real.(seq_complex)
end

"""
Per-unit Conversion: Convert actual values to per-unit

Parameters:
- actual: actual value
- base: base value

Returns:
- pu: per-unit value
"""
function to_pu(actual::Float64, base::Float64)
    return actual / base
end

"""
Actual Value Conversion: Convert per-unit to actual values

Parameters:
- pu: per-unit value
- base: base value

Returns:
- actual: actual value
"""
function from_pu(pu::Float64, base::Float64)
    return pu * base
end

"""
Phase-to-Line Voltage Conversion

Parameters:
- v_phase: phase voltage magnitude

Returns:
- v_line: line voltage magnitude
"""
function phase_to_line_voltage(v_phase::Float64)
    return √3 * v_phase
end

"""
Line-to-Phase Voltage Conversion

Parameters:
- v_line: line voltage magnitude

Returns:
- v_phase: phase voltage magnitude
"""
function line_to_phase_voltage(v_line::Float64)
    return v_line / √3
end

"""
Power Calculations in DQ frame

Parameters:
- v_d, v_q: d and q axis voltages
- i_d, i_q: d and q axis currents

Returns:
- P: active power
- Q: reactive power
"""
function calculate_power_dq(v_d::Float64, v_q::Float64, i_d::Float64, i_q::Float64)
    P = 1.5 * (v_d * i_d + v_q * i_q)  # Active power
    Q = 1.5 * (v_q * i_d - v_d * i_q)  # Reactive power
    return P, Q
end

"""
RMS Value Calculation from instantaneous ABC values

Parameters:
- abc: 3-element vector of instantaneous values

Returns:
- rms: RMS value
"""
function calculate_rms(abc::Vector{Float64})
    return sqrt(sum(abc.^2) / 3)
end

"""
Phasor Calculation from time-domain signal

Parameters:
- signal: time-domain signal vector
- t: time vector
- f: fundamental frequency

Returns:
- magnitude: phasor magnitude
- phase: phasor phase (radians)
"""
function calculate_phasor(signal::Vector{Float64}, t::Vector{Float64}, f::Float64)
    ω = 2π * f
    N = length(signal)
    
    # Calculate Fourier coefficients
    a = (2/N) * sum(signal .* cos.(ω .* t))
    b = (2/N) * sum(signal .* sin.(ω .* t))
    
    magnitude = sqrt(a^2 + b^2)
    phase = atan(b, a)
    
    return magnitude, phase
end

"""
Symmetrical Components Analysis

Parameters:
- Va, Vb, Vc: three-phase voltages (complex)

Returns:
- V0: zero sequence component
- V1: positive sequence component
- V2: negative sequence component
"""
function symmetrical_components(Va::ComplexF64, Vb::ComplexF64, Vc::ComplexF64)
    a = exp(im * 2π/3)  # Complex operator
    
    # Transformation matrix
    A = [1  1   1;
         1  a   a^2;
         1  a^2 a]
    
    # Calculate sequence components
    V_abc = [Va; Vb; Vc]
    V_012 = (1/3) * A * V_abc
    
    return V_012[1], V_012[2], V_012[3]  # V0, V1, V2
end

"""
Voltage Unbalance Factor Calculation

Parameters:
- Va, Vb, Vc: three-phase voltage magnitudes

Returns:
- unbalance_factor: voltage unbalance factor (%)
"""
function voltage_unbalance_factor(Va::Float64, Vb::Float64, Vc::Float64)
    # Convert to complex assuming phase angles
    Va_complex = Va + 0im
    Vb_complex = Vb * exp(-im * 2π/3)
    Vc_complex = Vc * exp(im * 2π/3)
    
    # Calculate sequence components
    V0, V1, V2 = symmetrical_components(Va_complex, Vb_complex, Vc_complex)
    
    # Unbalance factor = |V2| / |V1| * 100%
    if abs(V1) > 1e-6
        return abs(V2) / abs(V1) * 100.0
    else
        return 0.0
    end
end

"""
Total Harmonic Distortion (THD) Calculation

Parameters:
- signal: time-domain signal
- fs: sampling frequency
- f_fundamental: fundamental frequency

Returns:
- thd: Total Harmonic Distortion (%)
"""
function calculate_thd(signal::Vector{Float64}, fs::Float64, f_fundamental::Float64)
    # Perform FFT
    N = length(signal)
    fft_result = fft(signal)
    freqs = fftfreq(N, fs)
    
    # Find fundamental frequency bin
    fund_bin = argmin(abs.(freqs .- f_fundamental))
    fund_magnitude = abs(fft_result[fund_bin])
    
    # Calculate harmonic magnitudes (up to 50th harmonic)
    harmonic_sum = 0.0
    for h in 2:50
        harm_freq = h * f_fundamental
        if harm_freq < fs/2  # Nyquist limit
            harm_bin = argmin(abs.(freqs .- harm_freq))
            harmonic_sum += abs(fft_result[harm_bin])^2
        end
    end
    
    # THD calculation
    if fund_magnitude > 1e-6
        thd = sqrt(harmonic_sum) / fund_magnitude * 100.0
    else
        thd = 0.0
    end
    
    return thd
end