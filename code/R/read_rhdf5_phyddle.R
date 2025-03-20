# BiocManager::install("rhdf5")

library(rhdf5)

test=rhdf5::h5read("sim.test.hdf5","labels")

emp = rhdf5::h5read("sim.empirical.hdf5","label_names")

any(test == "")
sum(test == 8)

