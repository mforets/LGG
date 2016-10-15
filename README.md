# What is reachability?  

Reachability analysis consists, roughly speaking, in determining future (or past) trajectories of dynamical systems. We will work on the continuous-time setting, and with dynamical systems that are mathematically described by a set of ordinary differential equations. The reachability problem gets more interesting when we consider *uncertainties*, which usually involve one or more of these assumptions:

(i) uncertainty in the initial states; 

(ii) uncertainty in the coefficients of the system; 

(iii) uncertainties in the inputs of the system (e.g. noise). 

These assumptions lead to consider reachability flow-pipes (which encompass an infinite number of trajectories obtained by the numeric solution of usual ODE-solvers). A convenient computational representation of reachability flow-pipes is given by [convex polytopes](https://en.wikipedia.org/wiki/Convex_polytope).  

## What is LGG?

Reachability analysis consists in determining future (or past) trajectories of systems of differential equations. The systems are assumed to be uncertain, which usually involves the following assumptions: (i) uncertainty in the initial states; (ii) uncertainty in the coefficients of the system; (iii) uncertainties in the inputs of the system (e.g. noise). These assumptions lead to consider reachability flow-pipes (which encompass an infinite number of trajectories obtained by the numeric solution of usual ODE-solvers). A convenient computational representation of reachability flow-pipes is given by [convex polytopes](https://en.wikipedia.org/wiki/Convex_polytope).  

## What is LGG?

LGG stands for the name of the authors of the algorithm, Colas Le Guernic and Antoine Girard, see [linear](http://www.sciencedirect.com/science/article/pii/S1751570X09000387), [hybrid](http://link.springer.com/chapter/10.1007/978-3-642-02658-4_40), and references therein. It is a set-based reachability algorithm, that can be applied to uncertain hybrid systems with piecewise-linear dynamics. Moreover, it is one of the algorithms included in the [SpaceEx](http://spaceex.imag.fr/) state space explorer tool, it can be downloaded for free at http://spaceex.imag.fr/. 

## What can I learn in the tutorial?

In the folder ```/tutorial```, you will find a set of Jupyter notebooks with a tutorial. It contains further explanations of the implementation in SageMath, some worked examples, and exercises. Some but not all mathematical developments are included; I recommend that you refer to the original papers for further details. 

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
