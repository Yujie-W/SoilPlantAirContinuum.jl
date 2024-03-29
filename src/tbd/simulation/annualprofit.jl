#=
###############################################################################
#
# Calculate annual profit
#
###############################################################################
"""
    annual_profit(
                node::SPACSimple{FT},
                weather::Array{FT,2},
                Δt::FT = FT(1)
    ) where {FT<:AbstractFloat}

Calculate the profit in the growing season so as to optimize leaf investment,
    given
- `node` [`SPACSimple`] type struct
- `weather` Weather profile in a growing season
- `Δt` Time period in `[h]`
"""
function annual_profit(
            node::SPACSimple{FT},
            weather::Array{FT,2},
            Δt::FT = FT(1)
) where {FT<:AbstractFloat}
    # 0. unpack required values
    #    make sure construction cost must be postive
    @unpack c_cons, c_vmax, elevation, gaba, laba, latitude, vtoj = node;
    cons = laba * (c_cons + node.ps.PSM.v_cmax25 * c_vmax);

    # if cons > 0
    if cons > 0
        # 1. update the environmental constants based on node geographycal info
        ratio            = atmospheric_pressure_ratio(elevation);
        node.envir.P_AIR = ratio * P_ATM(FT);
        node.envir.P_O₂  = node.envir.P_AIR * FT(0.209);

        # 2. calculate the growing season canopy profit
        day  ::FT = FT(0);
        hour ::FT = FT(0);
        _tair::FT = FT(0);
        _dair::FT = FT(0);
        p_co2::FT = FT(0);
        r_all::FT = FT(0);
        wind ::FT = FT(0);
        rain ::FT = FT(0);
        gscp ::FT = FT(0);
        anet ::FT = FT(0);

        N = size(weather)[1]
        for i in 1:N
            # 2.1 read and update the hourly data
            day   = weather[i,2 ]
            hour  = weather[i,3 ]
            _tair = weather[i,7 ]
            _dair = weather[i,9 ]
            p_co2 = weather[i,10] * ratio
            r_all = weather[i,4 ]
            wind  = weather[i,6 ]
            rain  = weather[i,5 ]
            update_environment!(node.envir; p_CO₂ = p_co2, t = _tair + T_0(FT), vpd = FT(_dair) * 1000, wind = wind);

            # 2.2 if day time
            zenith = zenith_angle(latitude, day, hour, FT(0));
            if (r_all>0) & (zenith<=85)
                # 2.2.1 update the leaf partitioning
                big_leaf_partition!(node, zenith, r_all);
                @unpack frac_sh, frac_sl = node.container2L;

                # 2.2.2 optimize flows in each layer
                optimize_flows!(node);
                leaf_gas_exchange!(node, node.opt_f_sl, node.opt_f_sh);
                flow = node.opt_f_sl + node.opt_f_sh;
                anet = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;

                # 2.2.3 update drought history
                pressure_profile!(node.hs, node.p_soil, node.opt_f_sl, node.opt_f_sh, frac_sl);

            # 2.3 if night time
            else
                node.ps.t = max(200, leaf_temperature(node, r_all, FT(0)));
                photosystem_temperature_dependence!(node.ps.PSM, node.envir, node.ps.t);
                node.ps.p_H₂O_sat = saturation_vapor_pressure(node.ps.t);
                flow = FT(0);
                anet = -node.ps.PSM.r_d;
            end

            # 2.4 update soil moisture by converting flow to Kg per hour
            flow /= KG_H_2_MOL_S(FT);
            rain *= gaba * ρ_H₂O(FT) / 1000;
            soil_moisture!(node, flow-rain);

            # 2.5 add up the gscp
            gscp += anet
        end

        # 3. scale gscp to mol per year
        gscp *= laba * FT(1e-6) * 3600 * Δt;

        # 4. substract the construction costs
        gscp -= cons;

    else
        gscp = FT(-Inf);
    end

    return gscp
end
=#
