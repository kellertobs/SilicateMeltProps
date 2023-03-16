using SilicateMeltProps
using Test


println("Testing package SilicateMeltProps.")

println("Testing function DensityX.")

# calculate melt density of anhydrous N-MORB at 1200 C and 1.5 kbar
#         SiO2   TiO2  Al2Oe  FeO(t) MgO  CaO    Na2O  K2O    H2O
oxd_wt = [50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.140, 0.0]  # [wt%]
T = 1200.0  # [deg C]
P = 1.5     # [kbar]

println("Calculate density of anhydrous N-Morb at 1200 C, 1.5 kbar.")
rho = DensityX(oxd_wt,T,P)  # [kg/m3]

# compare to expected value (calculated with DensityX.xls)
println("Calculated density is: ",rho)
rho_exp = 2718.187297500990
println("Expected density is: ",rho_exp)
err = rho - rho_exp

# report test result
if err<1e-9
    println("Calculated density matches expected value, test passed.")
else
    error("Calculated density does not match expected value, test failed.")
end

# calculate melt density of hydrous N-MORB (1 wt% H2O) at 1200 C and 1.5 kbar
#         SiO2   TiO2  Al2Oe  FeO(t) MgO  CaO    Na2O  K2O    H2O
oxd_wt = [50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.140, 1.0]  # [wt%]
T = 1200.0  # [deg C]
P = 1.5     # [kbar]

println("Calculate density of anhydrous N-Morb at 1200 C, 1.5 kbar.")
rho = DensityX(oxd_wt,T,P)  # [kg/m3]

# compare to expected value (calculated with DensityX.xls)
println("Calculated density is: ",rho)
rho_exp = 2647.515936086210
println("Expected density is: ",rho_exp)
err = rho - rho_exp

# report test result
if err<1e-9
    println("Calculated density matches expected value, test passed.")
else
    error("Calculated density does not match expected value, test failed.")
end


println("Testing function Giordano08.")

# calculate melt viscosity of anhydrous N-MORB at 1200 C and 1.5 kbar
#         SiO2   TiO2  Al2Oe  FeO(t) MgO  CaO    Na2O  K2O    H2O
oxd_wt = [50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.140, 0.0]  # [wt%]
T = 1200.0  # [deg C]

println("Calculate viscosity of anhydrous N-Morb at 1200 C, 1.5 kbar.")
eta = Giordano08(oxd_wt,T)  # [Pas]

# compare to expected value (calculated with grdViscosity.xls)
eta_exp = 10^1.757695047576570
err = eta - eta_exp

# report test result
if err<1e-12
    println("Calculated viscosity matches expected value, test passed.")
else
    error("Calculated viscosity does not match expected value, test failed.")
end

# calculate melt viscosity of hydrous N-MORB (1 wt% H2O) at 1200 C and 1.5 kbar
#         SiO2   TiO2  Al2Oe  FeO(t) MgO  CaO    Na2O  K2O    H2O
oxd_wt = [50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.140, 1.0]  # [wt%]
T = 1200.0  # [deg C]

println("Calculate viscosity of anhydrous N-Morb at 1200 C, 1.5 kbar.")
eta = Giordano08(oxd_wt,T)  # [Pas]

# compare to expected value (calculated with grdViscosity.xls)
eta_exp = 10^1.05530134571660
err = eta - eta_exp

# report test result
if err<1e-12
    println("Calculated viscosity matches expected value, test passed.")
else
    error("Calculated viscosity does not match expected value, test failed.")
end