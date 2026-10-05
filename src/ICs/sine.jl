export Sine

struct Sine{D, M} <: SmoothInitialCondition
    a::State{M}
    period::Space{D}
    c_offset::State{M}
end

(ic::Sine)(pos::Space{D}) where {D} = ic.a * sin(2.0 * pi * sum(pos ./ ic.period)) + ic.c_offset

function build_ic(::Val{:sine}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:Type]::DataType, ctx[:D]::Int, ctx[:M]::Int
    
    a = State{M, T}(Tuple(T.(conf[:a])))
    period = Space{D, T}(Tuple(T.(conf[:period])))
    c_offset = State{M, T}(Tuple(T.(conf[:c_offset])))
    
    return Sine(a, period, c_offset)
end