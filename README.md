# LGG Reachability algorithm
Reachability analysis consists in determining future (or past) trajectories of systems of differential equations. The systems are assumed to be uncertain, which usually involves the following assumptions: (i) uncertainty in the initial states; (ii) uncertainty in the coefficients of the system; (iii) uncertainties in the inputs of the system (e.g. noise). These assumptions lead to consider reachability flow-pipes (which encompass an infinite number of trajectories obtained by the numeric solution of usual ODE-solvers). A convenient computational representation of reachability flow-pipes is given by [convex polytopes](https://en.wikipedia.org/wiki/Convex_polytope).  

## What this tutorial is about 

This is a tutorial on the LGG algorithm. LGG, that stands for the author's names Colas Le Guernic and Antoine Girard (see [linear](http://www.sciencedirect.com/science/article/pii/S1751570X09000387) and [hybrid](http://link.springer.com/chapter/10.1007/978-3-642-02658-4_40) papers), is a set-based reachability algorithm, that can be applied to uncertain hybrid systems with piecewise-linear dynamics. It is one of the algorithms included in the [SpaceEx](http://spaceex.imag.fr/) state space explorer tool, it can be downloaded for free at http://spaceex.imag.fr/. 

## What are the methods that we use 

Since we assume no familiarity with Sage, several explanations are included. There are also some exercises with solution to go further in the computational aspects. Some required notions in convex geometry are defined and computed, such as support functions and Hausdorff distance computations. We make frequent use of the following libraries:
* Polyhedra library 
* Linear programming
* Polynomial optimization, with the tool [ncpol2sdpa](http://peterwittek.github.io/ncpol2sdpa/)
* Numerical integration
* Random number generators
* Visualization (2d and 3d plots, explicit and implicit)

## How to read the project

The tutorial was designed in the collaborative mathematics platform SageMathCloud (SMC). It can also be run locally with SageMath + Jupyter notebook. Any of the notebooks (files with extension .ipynb) can be read with the [nbviewer](https://nbviewer.jupyter.org/), and searching for [this project](https://nbviewer.jupyter.org/github/mforets/LGG-Reachability-algorithm/tree/master/). 

## Who may be interested 

* If you would like to understand how does one of the algorithms in SpaceEx works.
* If you would like to start using Sage in your own work, and would like to see working examples of some Sage libraries (such as Polyhedra, LP), or to perform numerical computations with Numpy (such as performing numerical integrals or using random number generators).
* If you are interested in set-based methods in numerical analysis, reachability, and Hausdorff distance computations.
