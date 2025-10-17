module CMPSExcitations

using Reexport
@reexport using KrylovKit, OptimKit, LinearAlgebra, CMPSKit
@reexport using Plots, LaTeXStrings
@reexport using CMPSKit: YangGaudinCMPS

theme(:wong)
default(fontfamily="Computer Modern", label=nothing, dpi=100, framestyle=:box)

include("utils.jl")

export pitick, piticklabel, _find_groundstate, find_groundstate, enlarge_state

end