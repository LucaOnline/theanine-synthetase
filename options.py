"""The `options` module has several hardcoded values used to configure the program."""

# Number of clusters to use in mismatch clustering variance analysis.
# Each number here represents a separate set of trials.
CLUSTER_COUNTS = [5, 10, 15, 20, 30]

# Number of simulations to perform.
SIMULATION_COUNT = 1000

# Window sizes (in base pairs) to use during sliding-window dN/dS analysis.
# If the window size is too small, there can be situations in which dS is 0.
DNDS_WINDOW_SIZES = [180, 360, 540]
