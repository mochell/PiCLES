using ...Architectures: AbstractGrid, AbstractGridStatistics

function SphericalPropagationCorrection(phi; R=6.3710E+6)

    """
    Function to correct the change of angle due to the spherical propagation of the wave particle
    function GreatCircle_correction(cg_bar)
        phi Latitude in degrees
        R = 6.3710E+6 # earth Radius iun meters
        return (cg_bar * tan(phi)) / R
    end
    """
    coefficient = (sign.(phi) * min(sign.(phi) * tand.(phi), 60)) / R
    #coefficient = tand.(phi) / R

    function GreatCircle_correction(cg_x_bar)
        return cg_x_bar * coefficient
    end

    return GreatCircle_correction
end


function SphericalPropagationCorrection(phi, cg_x_bar; R=6.3710E+6)

    """
    function GreatCircle_correction(cg_x_bar)
        phi Latitude in degrees
        R = 6.3710E+6 # earth Radius iun meters
        return (cg_bar * tan(phi)) / R
    end
    """
    coefficient = (sign.(phi) * min(sign.(phi) * tand.(phi), 60)) / R
    #coefficient = tand.(phi) / R

    function GreatCircle_correction(cg_x_bar)
        return cg_x_bar * coefficient
    end

    return GreatCircle_correction(cg_x_bar)
end


"""
function SphericalPropagationCorrection(ij_mesh, gridstats)
    SphericalPropagationCorrection(gridd.data.y)
end
"""
function SphericalPropagationCorrection(ij_mesh::NamedTuple, gridstats::AbstractGridStatistics)
    SphericalPropagationCorrection(ij_mesh.y)
end


function SphericalPropagationCorrection_dummy(ij_mesh::NamedTuple, gridstats::AbstractGridStatistics)
    SphericalPropagationCorrection_dummy = x -> 0.0
    return SphericalPropagationCorrection_dummy
end
