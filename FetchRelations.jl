module FetchRelations
# fetch relations following JONSWAP

fetch_grouth_parameter= 3.5#3.622#3.5

"""
function X̃_j(X::Float64, U::Float64, g::Float=9.81)
    returns non-dimensional fetch using JONSWAP given X fetch (m), U (m/s) and gravity (m)
"""
function X̃_j(U::Number, X::Number, g::Number=9.81)
    g * X / U^2
end

"""
function X̃_j_U_freq(X::Float64, U::Float64, g::Float=9.81)
    returns non-dimensional fetch using JONSWAP given peak frequency f_max (Hz), U (m/s) and gravity (m)
"""
function X̃_j_U_freq(U10::Number, f_max::Number, g::Number=9.81)
    fetch_grouth_parameter^3.0815 * g^3 /(U10^3 *f_max^3)
end

"""
function X_j_U_freq(X::Float64, U::Float64, g::Float=9.81)
    returns fetch (meters) using JONSWAP given peak frequency f_max (Hz), U (m/s) and gravity (m)
"""
function X_j_U_freq(U10::Number, f_max::Number,  g::Number=9.81)
    fetch_grouth_parameter^3.0815 * g^2 /(U10 *f_max^3)
end

"""
function X̃_j_U_freq(X::Float64, U::Float64, g::Float=9.81)
    returns fetch (meters) using JONSWAP given peak frequency f_max (Hz), U (m/s) and gravity (m)
"""
function X̃_j_U_tau(U10::Number, tau::Number,  g::Number=9.81)
    (tau * g / (14 * pi * U10) )^(3/2)
end

"""
function α_j(f_max::Float64, U::Float64,, g::Float=9.81 )
    returns fetch amplitude using JONSWAP given peak frequency f_max (Hz), U (m/s) and gravity (m)
"""
function α_j(U::Number,  f_max::Number, g::Number=9.81 )
    0.033*(f_max*U/g)^0.67
end

"""
function τ_j(U::Number,  X::Number, g::Number=9.81 )
    returns equivalent fetch time using JONSWAP given U (m/s) and fetch (meters)and gravity (m)
"""
function τ_j(U::Number,  X::Number, g::Number=9.81 )
    14 * pi * (U/g) * X̃_j(U, X)^(2/3.0)
end

# tau1 = τ_j(U, X)
#
# U = 15.0
# X = 500e3
# X_tilde =X̃_j(U, X)
#
# X_tilde
#
# f_max = fₘ(U, X)
# f_max
#
# X/X_j_U_freq(U, f_max)
#
# X_j_U_freq(U, f_max)/X
#
# X_tilde/X̃_j_U_freq( U, f_max)
#
# fₘ(U, X )/ fₘ_given_U_tau(U, tau1)
#
# X̃_j_U_tau(U, tau1)
#
#



"""
function fₘ(U10::Float64, X::Float64)
    returns peak frequency using JONSWAP given U10 (m/s) and fetch (m)
"""
function fₘ(U10::Number, X::Number, g::Number=9.81)
    fetch_grouth_parameter* (g / U10) *  X̃_j(U10, X)^(-0.33)
end

"""
function fₘ_given_U_tau(U10::Float64, tau::Float64)
    returns peak frequency using JONSWAP given U10 (m/s) and time (sec)
"""
function fₘ_given_U_tau(U10::Float64, tau::Float64)
    X_tilde = (9.81 * tau/ ( 14 * pi * U10) )^(3.0/2.0)
    f_max   =  fetch_grouth_parameter * (9.81/ U10) * X_tilde^(-1/3.0)
    f_max * 1.035 # empirical adjustment factor
end

"""
function c_g_U_tau(U10::Float64, tau::Float64)
    returns peak frequency using JONSWAP given U10 (m/s) and time (sec)
"""
function c_g_U_tau(U10::Float64, tau::Float64)
    g= 9.81
    g / (4 * pi * fₘ_given_U_tau(U10, tau) )
end


"""
function Eⱼ(U10::Float64, tau::Float64)
    returns JONSWAP wave energy given U10 (m/s) and time (sec)
"""
function Eⱼ(U10::Float64, tau::Float64)
    f_max = fₘ_given_U_tau(U10::Float64, tau::Float64)
    alpha_j = α_j(U10,  f_max )
    0.31 * 9.81^2 * alpha_j * (f_max * 2 *pi)^(-4)
end



"""
JONSWAP_omega(U10::Float64, ωₚ , ω::Vector, g::Number=9.81)
    return the Pierson–Moskowitz spectrum parameters f_peak, Hs, Energy
    following (Bouws, 1998) given U10 (m/s) as the 10 meters wind speed.
"""
function JONSWAP_omega(U10::Float64, ωₚ , ω::Vector, g::Number=9.81)
    S = (2 * pi * α_j(U, ωₚ) * g^2 ) ./ ω.^5  .* exp.(-(5/4.0) .* (ωₚ./ω).^4 )
    sigma= (ω .> ωₚ ) .* 0.09 + (ω .<= ωₚ ) .* 0.07
    Gamma_j = exp.( - ( ω .- ωₚ).^2 ./ (2 .* sigma.^2 .* ωₚ^2) )
    S .*= 3.3.^Gamma_j
    S
end

# ω = collect(range(1/100, 1, step=1/100))
# ω.^3
# ωₚ = 2 * pi *f_max
#
# using Plots
# JONSWAP_omega(U, ωₚ , ω)
# plot( ω, JONSWAP_omega(U, ωₚ/2 , ω) )

"""
JONSWAP_frequency(U10::Float64, fₚ , freq::Vector, g::Number=9.81)
    return the Pierson–Moskowitz spectrum parameters f_peak, Hs, Energy
    following (Bouws, 1998) given U10 (m/s) as the 10 meters wind speed.
"""
function JONSWAP_frequency(U10::Float64, fₚ , freq::Vector, g::Number=9.81)
    JONSWAP_omega(U10, fₚ * 2 * pi , ω * 2 * pi)
end


#plot( ω/2 /pi, JONSWAP_frequency(U, ωₚ/2/pi/2 , ω/2/pi) )


"""
PMSpectrum(f::Vector, U10::Float64)
    returns Pierson Moskoviz spectum given f as the frequency in (Hz) and U10 (m/s) as the 10 meters wind speed.
    --- never tested!!
"""
function PMSpectrum(U10::Float64, f::Vector)
    """
    see Ocean Surface waves - S. R. Massel eq.3.79 and eq.3.80
    """
    g   = 9.81
    wp  = 0.879*g / U10
    w   = 2*np.pi*f
    sigma=0.04 *g / wp^2.0
    alpha=5.0* (wp^2.0 *sigma / g)^2.0
    alpha * w^(-5.0) * g^2.0 *np.exp(-5/4.0 * (w/wp)^-4)#
end

"""
PMParameters(f::Vector, U10::Float64)
    return the Pierson–Moskowitz spectrum parameters f_peak, Hs, Energy
    following (Bouws, 1998) given U10 (m/s) as the 10 meters wind speed.
    ---- never tested!!
"""
function PMParameters(U10::Float64)
    f_peak = 0.816  * 9.81 / (2 * np.pi * U10)
    Hs     = 0.0246 * U10^2
    E      = (Hs/4) ^4
    f_peak, Hs, E
end



end
