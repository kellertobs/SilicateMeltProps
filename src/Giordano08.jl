"""
    Giordano08(oxd_wt::Vector{_T}, T::_T) where {_T<:Number}

Calculate silicate melt viscosity `eta` after Iacovino & Till (2019)

Input: 
- `oxd_wt`  Melt composition as 9-element vector containing concentrations 
            in [wt%] of the following oxides ordered in the exact sequence 
            [SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O]
              1    2    3     4   5   6   7    8   9  
- `T`       Melt emperature in [degC]

Output:
- `eta`     Melt viscosity in [Pas]

Reference:
- Giordano D, Russell JK, & Dingwell DB (2008). Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science Letters, 271, 123-134. [doi:10.1016/j.epsl.2008.03.038](https://dx.doi.org/10.1016/j.epsl.2008.03.038)

Translated to Julia by Tobias Keller: [kellertobs @ GitHub](https://github.com/kellertobs), 03/2023
Original Matlab code by JK Russel: [FTP @ EOAS.ubc.ca](https://www.eoas.ubc.ca/~krussell/ftp/jkrftp.html), 03/2023

Example: calculate melt viscosity of hydrous N-MORB at 1200 C and 1.5 kbar
===
```julia
julia> oxd_wt = [50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.140, 1.0]  # [wt%]
julia> T      = 1200.0  # [deg C]
julia> eta    = Giordano08(oxd_wt,T)  # [Pas]
```
"""
function Giordano08(oxd_wt::Vector{_T}, T::_T) where {_T<:Number}

    # Fitting parameter values
    AT = -4.55
    bb = (159.56, -173.34, 72.13, 75.69, -38.98, -84.08, 141.54, -2.43, -0.91, 17.62)
    cc = (2.75, 15.72, 8.32, 10.20, -12.29, -99.54, 0.3)

    # Molar weights
    MW = (60.0855, 79.88, 101.96, 71.85, 40.3, 56.08, 61.98, 94.2, 18.02)

    # Convert from wt % to mol%
    # Note all data are normalized on anhydrous components first
    oxd_wt  = [oxd_wt[1:8].*(100.0 .- oxd_wt[9])./sum(oxd_wt[1:8]); oxd_wt[9]]
    oxd_mol = oxd_wt ./ MW
    oxd_mol = oxd_mol ./ sum(oxd_mol) .* 100

    # Load composition-basis matrix for multiplication against model-coefficients
    siti = oxd_mol[1] .+ oxd_mol[2]
    tial = oxd_mol[2] .+ oxd_mol[3]
    fmg  = oxd_mol[4] .+ oxd_mol[5]
    nak  = oxd_mol[7] .+ oxd_mol[8]
    b1   = siti
    b2   = oxd_mol[3]
    b3   = oxd_mol[4]
    b4   = oxd_mol[5]
    b5   = oxd_mol[6]
    b6   = oxd_mol[7] .+ oxd_mol[9]
    b7   = oxd_mol[9] .+ log.(1.0 .+ oxd_mol[9])
    b12  =  siti .* fmg
    b13  = (siti .+ oxd_mol[3]) .* (nak .+ oxd_mol[9])
    b14  = oxd_mol[3] .* nak

    c1   = oxd_mol[1]
    c2   = tial
    c3   = fmg
    c4   = oxd_mol[6]
    c5   = nak
    c6   = log.(1 .+ oxd_mol[9])
    c11  = (oxd_mol[3] .+ fmg .+ oxd_mol[6]) .* (nak .+ oxd_mol[9])
    bcf  = (b1, b2, b3, b4, b5, b6, b7, b12, b13, b14)
    ccf  = (c1, c2, c3, c4, c5, c6, c11)

    BT   = sum(bb .* bcf)
    CT   = sum(cc .* ccf)

    # Convert temperature to [K]
    TK   = T .+ 273.15

    # Calculate melt viscosity in [Pas]
    eta  = 10 .^ (min.(12, max.(-6, AT .+ BT ./(TK .- CT))))

    return eta

end