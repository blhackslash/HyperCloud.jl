export Sine

"""
    Sine{D, M} <: InitialCondition

Evaluates a smooth, multi-dimensional periodic sinusoidal wave.

# Fields
- `a::State{M}`: The amplitude state of the wave.
- `period::Space{D}`: The spatial periodicity along each axis.
- `c_offset::State{M}`: The constant background offset state.

The physical state is `a * sin(2π * sum(pos ./ period)) + c_offset`. It is evaluated
with `sinpi`, which keeps the input precision (no promotion to `Float64` via `pi`)
and is exact at integer multiples of the period.
"""
struct Sine{D, M} <: InitialCondition
    a::State{M}
    period::Space{D}
    c_offset::State{M}
end

(ic::Sine)(pos::Space{D}) where {D} = ic.a * sinpi(2 * sum(pos ./ ic.period)) + ic.c_offset

function build_ic(::Val{:sine}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:T]::DataType, ctx[:D]::Int, ctx[:M]::Int

    a = State{M, T}(Tuple(T.(conf[:a])))
    period = Space{D, T}(Tuple(T.(conf[:period])))
    c_offset = State{M, T}(Tuple(T.(conf[:c_offset])))

    return Sine(a, period, c_offset)
end
