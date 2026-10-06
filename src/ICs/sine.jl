export Sine

"""
    Sine{D, M} <: InitialCondition

Evaluates a smooth, multi-dimensional periodic sinusoidal wave.

# Fields
- `a::State{M}`: The amplitude state of the wave.
- `period::Space{D}`: The spatial periodicity along each axis.
- `c_offset::State{M}`: The constant background offset state.

The physical state is calculated mathematically as `a * sin(2.0 * pi * sum(pos ./ period)) + c_offset`.
"""
struct Sine{D, M} <: InitialCondition
    a::State{M}
    period::Space{D}
    c_offset::State{M}
end

(ic::Sine)(pos::Space{D}) where {D} = ic.a * sin(2.0 * pi * sum(pos ./ ic.period)) + ic.c_offset

function build_ic(::Val{:sine}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:T]::DataType, ctx[:D]::Int, ctx[:M]::Int
    
    a = State{M, T}(Tuple(T.(conf[:a])))
    period = Space{D, T}(Tuple(T.(conf[:period])))
    c_offset = State{M, T}(Tuple(T.(conf[:c_offset])))
    
    return Sine(a, period, c_offset)
end