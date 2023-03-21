
# CLINICAL RANGE
CLINICAL_RANGE_PERC = 0.8

#CONSTANT FOR DEPTH-DOSE PROFILE FITTING

CONVERSION_SIGMA = 2.35482

BRAGG_GAMMA_FUNC = 1.565 # no unit - gamma function value (Bortfeld model)
#BRAGG_ALPHA = 0.0022 #cm*MeV^{-p} in water
#BRAGG_P = 1.77 #no unit, in water
BRAGG_ALPHA_AL = 0.00146137480 #cm*MeV^{-p} - alpha parameter in the analytical relationship Energy - Range for aluminum
BRAGG_P_AL = 1.70525347 #no unit - p parameter in the analytical relationship Energy - Range for aluminum

BRAGG_ALPHA_CU = 0.00185679 #cm*MeV^{-p} - alpha parameter in the analytical relationship Energy - Range for copper
BRAGG_P_CU = 1.68973459 #no unit - p parameter in the analytical relationship Energy - Range for copper

BRAGG_ALPHA = (BRAGG_ALPHA_AL+BRAGG_ALPHA_CU)/2
BRAGG_P = (BRAGG_P_AL+BRAGG_P_CU)/2

#BRAGG_ALPHA = 0.001288253135486655 #cm*MeV^{-p} - alpha parameter in the analytical relationship Energy - Range
#BRAGG_P = 1.7293451613367292 #no unit - p parameter in the analytical relationship Energy - Range
BRAGG_BETA = 0.012 #cm^{-1} - beta parameter value (Bortfeld model)
BRAGG_GAMMA = 0.6 #no unit - gamma parameter value (Bortfeld model)
BRAGG_SIGMA_MONO_FACT = 0.012 # - (Bortfeld model)
BRAGG_SIGMA_MONO_EXP = 0.935 #no unit - (Bortfeld model)
BRAGG_SIGMA_E0_FACTOR = 0.01 #no unit - (Bortfeld model)
BRAGG_EPSILON_MIN = 0.0
BRAGG_EPSILON_MAX = 0.2

#QUBE CONSTANTS:
TO_WE = 0.234 # cm we per channel - each ionization chamber module corresponds to xx cm WEPL
TO_AL = 0.118 #cm aluminum equivalent per channel - each ionization chamber module corresponds to 0.118 cm AluminumEPL
