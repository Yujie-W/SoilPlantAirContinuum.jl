"""

    update_gsw!(spac::SPACMono{FT}, sm::EmpiricalStomatalModel{FT}, ind::Int, Δt::FT; β::FT = FT(1), τ::FT = FT(600)) where {FT<:AbstractFloat}

Update stomatal conductance, given
- `spac` Soil plant air continuum struct
- `sm` Empirical stomatal model
- `ind` Canopy layer index
- `Δt` Time step in `[s]`
- `β` Tuning factor on Anet (not on Vcmax)
- `τ` Time constant for prognostic stomtal conductance in `[s]`

"""
function update_gsw!(spac::SPACMono{FT}, sm::EmpiricalStomatalModel{FT}, ind::Int, Δt::FT; β::FT = FT(1), τ::FT = FT(600)) where {FT<:AbstractFloat}
    # calculate steady state values
    for _iLF in 1:spac.plant_ps[ind].n_leaf
        _gsw_ss = max(0, stomatal_conductance(sm, spac.plant_ps[ind], spac.envirs[ind], β, _iLF));
        spac.plant_ps[ind].g_sw[_iLF] += (_gsw_ss - spac.plant_ps[ind].g_sw[_iLF]) / τ * Δt;

        # bug fix: update g_lw and others as well
        spac.plant_ps[ind].g_lw[_iLF] = 1 / ( 1/spac.plant_ps[ind].g_sw[_iLF]  + 1/spac.plant_ps[ind].g_bw[_iLF] );
        spac.plant_ps[ind].g_sc[_iLF] = spac.plant_ps[ind].g_sw[_iLF] / FT(1.6);
        spac.plant_ps[ind].g_lc[_iLF] = 1 / ( 1/spac.plant_ps[ind].g_sc[_iLF] + 1/spac.plant_ps[ind].g_m[_iLF] + 1/spac.plant_ps[ind].g_bc[_iLF] );
    end;

    return nothing
end
