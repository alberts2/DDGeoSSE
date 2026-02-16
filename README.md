# Diversity-dependent GeoSSE (DDGeoSSE)

DDGeoSSE is a generative, geographical model for species diversification that accounts for state-dependence and diversity-dependence. As a generative model, DDGeoSSE can simulate phylogenetic trees with the following events: 

- **Within-region speciation** describes local speciation event in a particular region $i$.
- **Extinction describes** extinction event in a particular region $i$.
- **Dispersal** describes a dispersal event of a species from region $i$ to region $j$.
- **Between-region speciation** describes a between-region speciation where an ancestral range of a parent species splits into two daughter species with range $k$ and $\ell$.

As an inference model, DDGeoSSE is designed to work with [phyddle](https://phyddle.org/), which is a software for exploring phylogenetic models using deep learning. This allows us to perform two main inference tasks under our model: 

- **Parameter estimation**
- **Model selection**

## Diversity-dependent dynamics 

DDGeoSSE assigns a unique diversity-dependent effect parameter on each of the four events described above. These scalars take real values and can have positive, negative, or neutral correlation with local species richness. When the correlation is positive, the absolute rate of an event will increase with increasing species richness. When the correlation is negative, the absolute rate of an event will decreasing with increasing species richness. Lastly, when there is no correlation, the model is reduced to diversity-independent scenario where species richness has no influence on event rates. The model is flexible to test different hypotheses on the influence of species richness on diversification history. These diversiy-dependent dynamics then may induce local equilibrium diversity. This quantity is a result from the data-generating process, instead of being estimated separately from model parameters as commonly employed in other diversity-dependent models. Thus, avoiding the risk of statistical biases. 

## Local equilibrium diversity 

By definition, local equilibrium diversity represents a stationary state of the process, in which the number of local species diversity remains unchanged over time under a given model of diversification. Currently, to be able to derive this quantity, the model assumes that each region share the equal equilibrium diversity. However, this assumption can be relaxed if biologists have a good estimate on what these numbers for all the other regions in the system. 

## Simulation 

DDGeoSSE is written in [julia](https://julialang.org/) to simulate phylogenetic trees with tip states. To run the simulation script, enter the `jl` directory:
```shell
cd ~/DDGeoSSE/code/jl
```

Then run the `sim_one_func.jl` script using the following command: 

```shell
julia sim_one_func.jl simulate sim 0 1
```

This will simulate a complete phylogenetic tree with index 0  and store the output in the `simulate` directory. Currently, users can simulate trees with the following stopping conditions:

- **numtaxa** : number of extant taxa
- **time** : simulation run time
- **both**  : number of extant taxa and run time

Users can specify the number of regions by changing the value in `num_regions` in the simulation script. For each replicate, the script will return the following outputs: 

- `sim.*.tre`: a complete phylogenetic tree written as standard Newick string.
- `sim.*.labels.csv`: the true generating parameters of the replicate.
- `sim.*.extant.dat.csv`: tip states denoted by presence/absence of a species in regions.
- `sim.*.dat.csv`: tip states of species (including extinct species).
- `anagenetic_changes_tree_*.csv` : a history of anagenetic events on the tree.
- `cladogenetic event history_*.csv`: a history of cladogenetic events on the tree.
- `species_count_range_tree_*.csv` : a list containing number of extant species across geographical ranges at each sampled event time.
- `species_count_region_tree_*.csv` : a list containing number of extant species across different regions at each sampled event time.
- `cts_species_count_range_tree_*.csv`: a list containing number of extant species across geographical ranges in continuous time steps.
- `cts_species_count_region_tree_*.csv`: a list containing number of extant species across different regions in continuous time steps.
- `full_history_*.csv`: a list containing a complete history of the tree. 

## Phyddle analysis

Users can also run the simulation script directly from $\texttt{phyddle}$ by running the following script:
```shell
phyddle -s S -c config.py
```
Users can specify the size of the training dataset by modifying the `config.py` file (details on https://phyddle.org/overview.html#simulate)

To run the format step, users only require to provide `sim.*.tre`, `sim.*.labels.csv`, and `sim.*.dat.csv` in their training dataset as an input, and run the following command:
```shell
phyddle -s F -c config.py
```

Finally, users can run the rest of the pipeline, as follows:
```shell
# train network with tensor data
phyddle -s T -c config.py

# make prediction 
phyddle -s E -c config.py

# generate figures and store in plot
phyddle -s P -c config.py
```

To make prediction on parameter values, users need to specify the following parameters in `param_est` field of the `config.py` file, as follows 

- 'log_base_w'          : 'num' # This is the base rate for within-region speciation in log scale
- 'log_base_e'          : 'num', # This is the base rate for extinction in log scale
- 'log_base_d'          : 'num', # This is the base rate for dispersal in log scale 
- 'log_base_b'          : 'num', # This is the base rate for between-region speciation in log scale
- 'exp_scalar_w'        : 'num', # This is the diversity-dependent scalar for within-region speciation in exponential scale
- 'exp_scalar_e'        : 'num', # This is the diversity-dependent scalar for extinction in exponential scale
- 'exp_scalar_d_j'      : 'num' # This is the diversity-dependent scalar for incoming dispersal into region $j$ in exponential scale

For model selection, users need to specify the following parameters in `param_est` field of the `config.py` file, as follows 

- 'speciation_status'          : 'cat' # 0 = diversity-independent within-region speciation, 1 = diversity-dependent within-region speciation 
- 'extinction_status'          : 'cat' # 0 = diversity-independent extinction, 1 = diversity-dependent extinction
- 'dispersal_status'           : 'cat' # 0 = diversity-independent incoming dispersal, 1 = diversity-dependent incoming dispersal

## Important R scripts 

- `solve_balance.R`       : this script solves for the local equilibrium diversity numerically given parameters under a fully unconstrained DDGeoSSE model.
- `plot_reg_occupancy.R`  : this script visualizes region occupancy through time from a DDGeoSSE simulation.
- `plot_range_occupancy.R` : this script visualizes number of species in each range through time from a DDGeoSSE simulation.

## Important directories 
- [empiric_anolis](https://github.com/alberts2/DDGeoSSE/tree/main/data/empiric_anolis/use_this) contains Anolis phylogeny and biogeographic range data.
- [empiric_viburnum_rescaled](https://github.com/alberts2/DDGeoSSE/tree/main/data/empiric_viburnum_rescaled) contains Viburnum phylogeny and biogeographic range data.
- [theory_validation](https://github.com/alberts2/DDGeoSSE/tree/main/data/theory_validation) contains simulated dataset used to validate theories related to local equilibrium diversity (Figure 3).
- [simulator_validation](https://github.com/alberts2/DDGeoSSE/tree/main/data/simulator_validation/sim_geosse) contains simulated dataset used to validate the simulator under a GeoSSE model (Supplemental Figures 23-24).
- [conceptual_fig](https://github.com/alberts2/DDGeoSSE/tree/main/data/conceptual_fig) contains example to visualize DDGeoSSE (Figure 1)
- [deep_learning](https://github.com/alberts2/DDGeoSSE/tree/main/data/deep_learning/server/DD-GeoSSE/5_reg/FINAL/mixture_dataset) contains trained neural networks used for analyses on the manuscript. 

## Miscellaneous

- For model selection, we recommmend users to train a network to make a prediction on one diversity-dependent feature at a time, based on the finding we found on the manuscript.
- The `rel_extant_age_tol` in the `config.py` can cause incorrectly pruned tree of only extant species, if the number specified is too small. We recommend assigning large value when dealing with simulated data, and small value when dealing with empirical data.
- If you have any questions, please e-mail me at [soewongsono@wustl.edu](mailto:soewongsono@wustl.edu).
                          
                           
                            
                           
                            

