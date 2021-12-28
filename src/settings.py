# This is a setting file containing some parameters used in the hepParser.

# The section of mz data simulator

# the min and max dp of theory hp for fitting the isotopic peaks ratio.
min_dp = 2
max_dp = 20

# the number of isotopic peaks in one isotopic cluster
iso_num = 10

# the min and max mass that can be simulate. (not m/z range)
min_mass = 400
max_mass = 6000

# the number of hp that be simulated (not the number of all simulated data)
hp_num = 100

# the max row range that one m/z can shift in simulator
max_row_num = 20

# the max scans number that one distribution can cover
max_columns_num = 400

# the max bias of m/z (ppm)
ppm = 20

# the max charge that one hp can get in hepParser
max_charge = 10

# the max mean num of noise mz between two adjust isotopic peaks (charge=1)
max_noise_mz_num = 20

# the std of noise mz num which follows the gauss distribution
noise_mz_num_std = 0.2

# the ratio that noise/(matched peak)
# which mean that the most intensive noise is no more than 0.01 times of max intensity peak
noise_int_ratio = 0.01

# the ratio that one position intensity was set to 0
p_set0 = 0.1

# the max number of peaks follow normal distribution in noise simulator
max_num_norm = 4
