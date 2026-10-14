export Sine

"""
    Sine{D, M, T} <: InitialCondition

Evaluates a smooth, multi-dimensional periodic sinusoidal wave.

# Fields
- `a::State{M, T}`: The amplitude state of the wave.
- `period::Space{D, T}`: The spatial periodicity along each axis.
- `c_offset::State{M, T}`: The constant background offset state.

The physical state is `a * sin(2π * sum(pos ./ period)) + c_offset`. It is evaluated
with `sinpi`, which keeps the input precision (no promotion to `Float64` via `pi`)
and is exact at integer multiples of the period.
"""
struct Sine{D, M, T} <: InitialCondition
    a::State{M, T}
    period::Space{D, T}
    c_offset::State{M, T}
end

function Sine(a::State{M}, period::Space{D}, c_offset::State{M}) where {D, M}
    T = promote_type(eltype(a), eltype(period), eltype(c_offset))
    return Sine{D, M, T}(a, period, c_offset)
end

(ic::Sine)(pos::Space{D}) where {D} = ic.a * sinpi(2 * sum(pos ./ ic.period)) + ic.c_offset

function build_ic(::Val{:sine}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:T]::DataType, ctx[:D]::Int, ctx[:M]::Int

    a = State{M, T}(Tuple(T.(conf[:a])))
    period = Space{D, T}(Tuple(T.(conf[:period])))
    c_offset = State{M, T}(Tuple(T.(conf[:c_offset])))

    return Sine(a, period, c_offset)
end
