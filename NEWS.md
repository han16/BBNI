# BBNI 0.2.1

* Major performance optimization: ~14x speedup via vectorization in Error_LLH and implementing Floyd-Warshall algorithm in update_ancestor_matrix, keeping strict numerical equivalence with v0.1.1
* Vignette expanded and successfully compiled to demonstrate new independent (non-timeseries) mode and visualization features
* Real-world yeast dataset application realized in the vignette
* Minor code reformatting for readability

# BBNI 0.2.0

* Added new visualization functions: plot_bbni(), plot_trace(), and plot_network()
* Enhanced plot_bbni() to compare inferred networks against true networks and fixed a reversed edge direction bug
* Implemented independent (non-timeseries) mode across core algorithm and data generation functions
* Upgraded run_bbni() with a progress bar, MCMC summary, burn-in parameters, and posterior edge probabilities
* Optimized MCMC mixing with logic fixes to ProposalConstruction
* Added default parameters for key user-facing functions
* Significantly expanded documentation and examples across all primary functions
* Included public yeast dataset from original paper for user testing and for vignette

# BBNI 0.1.1

* Rewrote documentation, vignette, and README for clarity
* Reformatted code for readability
* Removed unused/dead code/comments
* Fixed spelling and minor typos

# BBNI 0.1.0

* Initial development version.
* Refactored legacy Bayesian Boolean Network Inference code into a modular, documented R package.
* Added `run_bbni()` as the primary user-facing function.
* Added a vignette demonstrating network recovery from simulated data.
* Added unit tests for core network-validity and likelihood functions.
