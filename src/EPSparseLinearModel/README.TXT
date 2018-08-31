------------------------------------------------------------------- EXPERIMENTAL DESIGN FOR EFFICIENT IDENTIFICATION OF GENE REGULATORY
 NETWORKS USING SPARSE BAYESIAN MODELS

 Florian Steinke, Matthias Seeger, Koji Tsuda
 
 Max-Planck Institute for Biological Cybernetics
 Spemannstrasse 38
 P.O. Box 21 69
 72012 Tuebingen, Germany
 steinke,seeger,tsuda@tuebingen.mpg.de
-------------------------------------------------------------------

Contents and Usage
==================

This package includes two parts:

1) In directory ToyExample, we present a very simple toy system that highlights the use and the advantages of the sparse linear model and its Bayesian treatment. We demonstrate experimental design on simple linear regression from R to R. The example also show the usage of the iAUC evaluation metric.

To run the example type the following in Matlab

	cd ToyExample
	run_example
	
More details about the example model formulation are found in run_example.m.	
	
2) In directory BMCPaper you will find the code for reproducing the essential results of the paper. For better readability of the code, we have removed some options and extensions in the code that were used for specific evaluations and comparisons. We hope this way, the method is more easy to apply for others.

To reproduce the central figure 3 in the paper, type the following in Matlab

	cd BMCPaper
	runmex
	main_generate_networks
	main('prepare','design')
	main('run','design')
	main('eval','design')
	main('plot','design')
	
To reproduce the figure 4 left in the paper, type the following in Matlab

	cd BMCPaper
	runmex
	main_generate_networks
	main('prepare','inputshape')
	main('run','inputshape')
	main('eval','inputshape')
	main('plotclass','inputshape',2)
	
Note that these calls may take a very long time, since many repetitions are performed for statistical soundness. Note also that a large part of the time is spent for sampling from the stochastic differential equations and for evaluation with the ROC measure. This would not be necessary in a real experiment.


Prerequisites
=============

This implementation is based on Matlab with some core routines included into Matlab as MEX files written in C++. The package requires Matlab Version 7 or higher, and MEX file binaries are included in this package for 32-bit Windows and Linux.

One may be able to build these mex files also for other architectures, but the source code has to be downloaded separately from

	http://www.kyb.tuebingen.mpg.de/sparselinearmodel

No support is provided for building such binaries.

Usage
=====

The code is published under the GNU GPL licence. If you want to use it in a scientific publication, please cite

@article{SteSeeTsu07,
	author	={Steinke, Florian. And Seeger, Matthias. and Tsuda, Koji},
	title	= {Experimental Design for Efficient Identification of Gene Regulatory Networks using Sparse Bayesian Models},
	journal	= {BMC Systems biology},
	note 	= {(under review)},
	year 	= {2007},
}



