# BiocManager::install("rhdf5")

library(rhdf5)

test=rhdf5::h5read("sim.test.hdf5","labels")

any(test == "")
sum(test == 8)

