#==============================================================================#
# Config:       phyddle config file DDGeoSSE                                    #
# Authors:     Albert Soewongsono and Michael Landis                          #
# Date:         250627                                                         #
# Description:  DDGeoSSE for testing density-dependent within-region speciation  #
#==============================================================================#


args = {

    #-------------------------------#
    # Project organization          #
    #-------------------------------#
    'step'    : 'SFTEP',                   # Step(s) to run
    'verbose' : 'T',                       # print verbose phyddle output?
    'prefix'  : 'sim',                     # Prefix for output for all setps
    'dir'     : './',
    #'output_precision'   : 12,             # Number of digits (precision) for numbers in output files

    #-------------------------------#
    # Multiprocessing               #
    #-------------------------------#
    'use_parallel'   : 'T',                 # use multiprocessing to speed up jobs?
    'use_cuda'       : 'T',
    'num_proc'       : 90,                  # how many CPUs to use (-2 means all but 2)

    #-------------------------------#
    # Simulate Step settings        #
    #-------------------------------#
    'sim_command'       : f'julia sim_one_func.jl', # exact command string, argument is output file prefix
    'sim_logging'       : 'verbose',        # verbose, compressed, or clean
    'start_idx'         : 1,                # first simulation replicate index
    'end_idx'           : 400000,             # last simulation replicate index
    'sim_batch_size'    : 1000,

    #-------------------------------#
    # Format Step settings          #
    #-------------------------------#
    'encode_all_sim'    : 'T',
    'num_char'          : 5,                # number of evolutionary characters
    'num_states'        : 2,                # number of states per character
    'min_num_taxa'      : 10,               # min number of taxa for valid sim
    'max_num_taxa'      : 60000,              # max number of taxa for valid sim
    'rel_extant_age_tol': 1E-4,
    'tree_width'        : 300,              # tree width category used to train network
    'tree_encode'       : 'extant',         # use model with serial or extant tree
    'brlen_encode'      : 'height_brlen',   # how to encode phylo brlen? height_only or height_brlen
    'char_encode'       : 'integer',        # how to encode discrete states? one_hot or integer
    'param_est'         : {                 # model parameters to predict (labels)
                            'speciation_status'          : 'cat'
                            },
    'param_data'        : {},                 # model parameters that are known (aux. data)
    'tensor_format'     : 'hdf5',           # save as compressed HDF5 or raw csv
    'char_format'       : 'csv',
    'save_phyenc_csv'   : 'F',              # save intermediate phylo-state vectors to file

    #-------------------------------#
    # Train Step settings           #
    #-------------------------------#
    'num_epochs'        : 200,              # number of training intervals (epochs)
    'prop_test'         : 0.05,             # proportion of sims in test dataset
    'prop_val'          : 0.05,             # proportion of sims in validation dataset
    'prop_cal'          : 0.20,             # proportion of sims in CPI calibration dataset
    'cpi_coverage'      : 0.80,             # coverage level for CPIs
    'cpi_asymmetric'    : 'T',              # upper/lower ('T') or symmetric ('F') CPI adjustments
    'trn_batch_size'    : 4096,             # number of samples in each training batch
    'loss_numerical'    : 'mae',            # loss function for learning (real-valued) targets
    'optimizer'         : 'adam',           # optimizer for network weight/bias parameters
    'phy_channel_plain'  : [32, 64],        # Output channel sizes for plain convolutional layers for phylogenetic state input
    'phy_channel_stride' : [32, 64],             # Output channel sizes for stride convolutional layers for phylogenetic state input
    'phy_channel_dilate' : [32, 64],             # Output channel sizes for dilate convolutional layers for phylogenetic state input
    'aux_channel'        : [64, 32],        # Output channel sizes for dense layers for auxiliary data input
    'lbl_channel'        : [64, 32, 16, 8],        # Output channel sizes for dense layers for label outputs
    'phy_kernel_plain'   : [3, 5],            # Kernel sizes for plain convolutional layers for phylogenetic state input
    'phy_kernel_stride'  : [7, 9],               # Kernel sizes for stride convolutional layers for phylogenetic state input
    'phy_kernel_dilate'  : [3, 5],               # Kernel sizes for dilate convolutional layers for phylogenetic state input
    'phy_stride_stride'  : [3, 6],               # Stride sizes for stride convolutional layers for phylogenetic state input
    'phy_dilate_dilate'  : [3, 5],
    'early_stop_rule'   : 'total',
    'num_early_stop'    : 20,
    'loss_aggregation'  :'sum',
    'learning_rate'     : '0.0010',
    #-------------------------------#
    # Estimate Step settings        #
    #-------------------------------#

    #-------------------------------#
    # Plot Step settings            #
    #-------------------------------#
    'plot_train_color'      : 'blue',       # plot color for training data
    'plot_test_color'       : 'purple',     # plot color for test data
    'plot_val_color'        : 'red',        # plot color for validation data
    'plot_aux_color'        : 'green',      # plot color for input auxiliary data
    'plot_label_color'      : 'orange',     # plot color for labels (params)
    'plot_emp_color'        : 'black'      # plot color for estimated data/values
    #'plot_pca_noise'        : 0.01          # Add noise to PCA plot to dampen contrast from point values
 }
