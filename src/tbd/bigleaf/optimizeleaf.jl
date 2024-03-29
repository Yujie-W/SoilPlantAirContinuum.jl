#=
###############################################################################
#
# Optimize leaf invstment Vcmax+Jmax and Leaf area
#
###############################################################################
"""
    optimize_leaf!(
                node::SPACSimple{FT},
                weather::Array{FT,2},
                printing::Bool,
                Δt::FT = FT(1)
    ) where {FT<:AbstractFloat}

Optimize leaf area (LAI within 0-20) and photosynthetic capacity (within
    5-200), given
- `node` [`SPACSimple`] type struct
- `weather` Weather profile in a growing season
- `printing` Optional. If true, printing progress
- `Δt` Time period in `[h]`
"""
function optimize_leaf!(
            node::SPACSimple{FT},
            weather::Array{FT,2},
            printing::Bool = false,
            Δt::FT = FT(1)
) where {FT<:AbstractFloat}
    # 1. use the opt_laba and opt_vmax to initialize
    @inline f(x) = (
        tmp_node = deepcopy(node);
        leaf_allocation!(tmp_node, x[1], x[2]);
        tmp_prof = annual_profit(tmp_node, weather, Δt);
        if printing
            @info "\tLABA=$(x[1])\tVMAX=$(x[2])\tPROF=$(tmp_prof)";
        end;
        return tmp_prof
    );

    ms = ReduceStepMethodND{FT}(
                x_mins=FT[1,5],
                x_maxs=FT[node.gaba*20, 200],
                x_inis=FT[node.opt_laba, node.opt_vmax],
                Δ_inis=FT[100,10]);
    st = SolutionToleranceND{FT}(FT[0.9, 0.09], 50);
    lv = find_peak(f, ms, st);

    # 2. update the opt_laba and opt_vmax to node
    leaf_allocation!(node, lv[1], lv[2]);
    node.opt_laba = lv[1];
    node.opt_vmax = lv[2];

    return nothing
end
=#
