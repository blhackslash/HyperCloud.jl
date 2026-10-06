export build_tableau, build_timestepper

# --- Tableau Builder ---
"""
    build_tableau(time_conf::Dict, context::Dict)

Parses the user-specified `:stepper` string and builds the corresponding Butcher tableau.
Supports explicit Runge-Kutta tableaus (e.g., `:Euler`, `:RK2`, `:RK3`, `:RK4`) and implicit-explicit (IMEX) tableaus (e.g., `:ARS222`, `:ARS233`, `:PRSSP3`, `:SSP332`, `:IMEXEuler`).
"""
function build_tableau(time_conf::Dict, context::Dict)
    name = time_conf[:stepper]::Symbol
    return build_tableau(Val(name), time_conf, context)
end

"""
    build_tableau(name::Val, conf::Dict, ctx::Dict)

Generic fallback for Butcher tableaus.
Users can extend this by defining `build_tableau(::Val{:my_stepper}, ...)`.
Supports explicit Runge-Kutta tableaus (e.g., `:Euler`, `:RK2`) and implicit-explicit (IMEX) tableaus (e.g., `:ARS222`, `:IMEXEuler`).
"""
build_tableau(name::Val, conf::Dict, ctx::Dict) = error("Unknown Time Stepper Tableau: $(typeof(name))")

# Explicit RK Tableaus
build_tableau(::Val{:Euler}, conf::Dict, ctx::Dict) = RK1_Euler_Tableau(ctx[:T])
build_tableau(::Val{:RK2}, conf::Dict, ctx::Dict)   = RK2_Ralston_Tableau(ctx[:T])
build_tableau(::Val{:RK3}, conf::Dict, ctx::Dict)   = RK3_SSP_Tableau(ctx[:T])
build_tableau(::Val{:RK4}, conf::Dict, ctx::Dict)   = RK4_Classical_Tableau(ctx[:T])

# IMEX Tableaus
build_tableau(::Val{:ARS233}, conf::Dict, ctx::Dict) = IMEX_ARS233_Tableau(ctx[:T])
build_tableau(::Val{:PRSSP3}, conf::Dict, ctx::Dict) = IMEX_PRSSP3_Tableau(ctx[:T])
build_tableau(::Val{:ARS222}, conf::Dict, ctx::Dict) = IMEX_ARS222_Tableau(ctx[:T])
build_tableau(::Val{:SSP332}, conf::Dict, ctx::Dict) = IMEX_SSP2332_Tableau(ctx[:T])
build_tableau(::Val{:IMEXEuler}, conf::Dict, ctx::Dict) = IMEX_Euler_Tableau(ctx[:T])

# --- Main Time Stepper Builder ---
"""
    build_timestepper(context::Dict)

Deduces and builds the overarching time integrator based on the constructed tableau.
Automatically returns a `GeneralRKTimeStepper` if an `RKButcherTableau` is present in the context, or a `GeneralIMEXTimeStepper` if an `IMEXButcherTableau` is detected.
"""
function build_timestepper(context::Dict)
    tableau = context[:Tableau]
    return _build_timestepper(tableau, context)
end

# Deduce Explicit RK Stepper from the RKButcherTableau
function _build_timestepper(tableau::RKButcherTableau, context::Dict)
    return GeneralRKTimeStepper(
        context[:Equation], 
        context[:Scheme], 
        context[:MOOD],
        context[:ExplicitSources], 
        tableau
    )
end

# Deduce IMEX Stepper from the IMEXButcherTableau
function _build_timestepper(tableau::IMEXButcherTableau, context::Dict)
    return GeneralIMEXTimeStepper(
        context[:Equation], 
        context[:Scheme], 
        context[:MOOD],
        context[:ExplicitSources], 
        context[:ImplicitSources], 
        tableau
    )
end