# plotting utilities
function piticklabel(x::Rational, ::Val{:latex})
    iszero(x) && return L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return L"%$S%$N\pi"
    L"%$S\frac{%$N\pi}{%$d}"
end

function pitick(start, stop, denom; mode=:text)
    # Compute lower index
    a = Int(floor(start / (π / denom)))
    # Compute upper index
    b = Int(ceil(stop / (π / denom)))
    tick = range(a * π / denom, b * π / denom; step=π / denom)
    ticklabel = piticklabel.((a:b) .// denom, Val(mode))
    tick, ticklabel
end

# Hsingle(c, μ) = ∫(∂ψ̂' * ∂ψ̂ - μ * ψ̂' * ψ̂ + c * (ψ̂')^2 * ψ̂^2, (-Inf, +Inf));
Hsingle(c, μ) = ∫(2 * ∂ψ̂' * ∂ψ̂ - 2 * μ * ψ̂' * ψ̂ + 4 * c * (ψ̂')^2 * ψ̂^2, (-Inf, +Inf));
Hcoupled(c, μ) = ∫(
    (∂ψ̂₁' * ∂ψ̂₁ - μ * ψ̂₁' * ψ̂₁ + c * (ψ̂₁')^2 * ψ̂₁^2 +
     ∂ψ̂₂' * ∂ψ̂₂ - μ * ψ̂₂' * ψ̂₂ + c * (ψ̂₂')^2 * ψ̂₂^2 +
     2 * c * (ψ̂₁') * (ψ̂₂') * ψ̂₂ * ψ̂₁), (-Inf, +Inf));

function _find_groundstate(D, H, cT=YangGaudinCMPS; gradtol=1e-10, optalg=LBFGS(80; verbosity=0, maxiter=3000, gradtol=gradtol), state=nothing)
    if isnothing(state)
        T = Float64 # scalar type
        Q₀ = randn(T, (D, D))
        R₀ = randn(T, (D, D))
        state = cT(Constant(Q₀), Constant(R₀))
    end

    state, ρL, ρR, E, e, normgrad, numfg, history = CMPSKit.groundstate_unconstrained(H, state; optalg=optalg, verbosity=1, gradtol=gradtol)
    println("D = $D | $(typeof(state))")
    return state
end

function enlarge_state(state, Df; k=30, ϵ=0.0, T=YangGaudinCMPS)
    Qi, Ri = state.Q[], state.Rs[1][]
    Di = size(Qi, 1)
    dD = Df - Di
    Qf = zeros(eltype(Qi), Df, Df)
    Rf = zeros(eltype(Ri), Df, Df)

    # Q → Q ⊕ -kI
    Qf[1:Di, 1:Di] .= Qi
    Qf[Di+1:end, Di+1:end] .= -k * I(dD)

    # R → R ⊕ 0
    Rf[1:Di, 1:Di] .= Ri

    Qf .+= ϵ * rand(size(Qf))
    Rf .+= ϵ * rand(size(Rf))

    return T(Constant(Qf), Constant(Rf))
end

function find_groundstate(Ds, H, T=YangGaudinCMPS; state=nothing, k=50.0, kwargs...)
    if !isnothing(state)
        Ds = [size(state.Q[], 1)]
    end

    for D in Ds
        println("Optimizing D=$D")
        @time state = _find_groundstate(D, H, T; state=isnothing(state) ? state : enlarge_state(state, D, k=k, T=T), kwargs...)
        println("---------------")
    end

    return state
end