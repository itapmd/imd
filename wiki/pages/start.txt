~~NOTOC~~
====== Welcome to the IMD wiki ======

//IMD// is a software package for classical molecular dynamics simulations. Several types of interactions are supported, such as central pair potentials, EAM potentials for metals, Stillinger-Weber and Tersoff potentials for covalent systems, and Gay-Berne potentials for liquid crystals. A rich choice of simulation options is available: different integrators for the simulation of the various thermodynamic ensembles, options that allow to shear and deform the sample during the simulation, and many more. There is no restriction on the number of particle types.

The main design goals were to create a flexible and modular software reaching high performance on contemporary computer architectures, while being as portable as reasonably possible. //IMD// runs efficiently on both single processor workstations and massively parallel supercomputers, but is not very well suited for vector processors. The parallelization uses the standard [[http://www.mcs.anl.gov/research/projects/mpi/|MPI]] library for message passing. As //IMD// is written in ANSI C, it should be easily portable to any Unix-like environment.

//IMD// even received a [[http://parallel.rz.uni-mannheim.de/sc/suparcup/win97/|prize]], and holds a [[http://imd.itap.physik.uni-stuttgart.de/worldrecord.html|world record]].

===== Collaborations and new developments =====
  * //IMD// @ [[http://www.icm.edu.pl/~rudnicki|GPU]]
  * //IMD// with [[https://openkim.org|OpenKIM]]
  * //IMD// and [[http://www.scafacos.de|ScaFaCos]]: Interface programmed by [[http://www.math.uni-bielefeld.de/~gaehler/|Franz Gähler]]
  * Load balancing //IMD//: Implementation by [[http://www.icams.de/content/people/staff-members/index.html?view=details&member=684|Christoph Begau]]
  * //IMD// combined with Mesoscale Simulations: [[http://www.hlrs.de/organization/people/brinkmann/|Steffen Brinkmann]]


===== Download =====

//IMD// is available as open source software released under the [[license|GPL]].

The latest version can be downloaded [[download|here]].


**We are currently reconstructing //IMD//. In the cource of this process we will also update the user guide.**