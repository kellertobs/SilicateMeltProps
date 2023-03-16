using SilicateMeltProps
using Test


@testset "DensityX" begin

# calculate melt density of anhydrous N-MORB at 1200 C and 1.5 kbar
#         SiO2   TiO2  Al2Oe  FeO(t) MgO  CaO    Na2O  K2O    H2O
oxd_wt = [50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.140, 0.0]  # [wt%]
T = 1200.0  # [deg C]
P = 1.5     # [kbar]

rho = DensityX(oxd_wt,T,P)  # [kg/m3]

# compare to expected value (calculated with DensityX.xls)
rho_exp = 2718.187297500990
@test rho ≈ rho_exp


# calculate melt density of hydrous N-MORB (1 wt% H2O) at 1200 C and 1.5 kbar
#         SiO2   TiO2  Al2Oe  FeO(t) MgO  CaO    Na2O  K2O    H2O
oxd_wt = [50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.140, 1.0]  # [wt%]
T = 1200.0  # [deg C]
P = 1.5     # [kbar]

rho = DensityX(oxd_wt,T,P)  # [kg/m3]

# compare to expected value (calculated with DensityX.xls)
rho_exp = 2647.515936086210
@test rho ≈ rho_exp #atol=1.0e-12

end


@testset "Giordano08" begin

# calculate melt viscosity of anhydrous N-MORB at 1200 C and 1.5 kbar
#         SiO2   TiO2  Al2Oe  FeO(t) MgO  CaO    Na2O  K2O    H2O
oxd_wt = [50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.140, 0.0]  # [wt%]
T = 1200.0  # [deg C]

eta = Giordano08(oxd_wt,T)  # [Pas]

# compare to expected value (calculated with grdViscosity.xls)
eta_exp = 10^1.757695047576570
@test eta ≈ eta_exp


# calculate melt viscosity of hydrous N-MORB (1 wt% H2O) at 1200 C and 1.5 kbar
#         SiO2   TiO2  Al2Oe  FeO(t) MgO  CaO    Na2O  K2O    H2O
oxd_wt = [50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.140, 1.0]  # [wt%]
T = 1200.0  # [deg C]

eta = Giordano08(oxd_wt,T)  # [Pas]

# compare to expected value (calculated with grdViscosity.xls)
eta_exp = 10^1.05530134571660
@test eta ≈ eta_exp atol=1.0e-12

end
