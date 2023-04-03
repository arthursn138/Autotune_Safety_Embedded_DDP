# Autotune Barrier Methods for Safe Trajectory Optimization
Aim to autotune the weights of barrier states and penalty for control barrier functions such that an optimal trajectory can be calculated with quick convergence, i.e. we avoid the conflict of the optimal path with too conservative barrier penalties.


Arthur Scaquetti do Nascimento (nascimento@gatech.edu), Hassan Almubarak (halmubarak@gatech.edu) - ACDS Lab @ Georgia Tech

Last Update May 2022

#
Follow those steps to run any example
1. Call the system's dynamics
2. Generate the obstacle course and define the safe set function (h)
3. Call DBaS_dyn to generate the DBaS dynamics
4. Call Safety_Embedding_dynamics to augment the DBaS to the system's dynamics
5. Define DDP and optimization paramters and run vanilla DDP

NOTE: to avoid pentrating the obstacles in some cases due to the discrete formulation, use disc_ddp_alg_penalty which penalizes the interior of the unsafe regions as done in the penalty methods. disc_ddp_alg_penalty takes h as an extra input to penalize the interiors of the obstacles.
This gives advantages of good planning from safety embedded ddp and the advantage of penalty methods. 


**Safety Embedded Differential Dynamic Programming using Discrete Barrier States (DBaS)**
https://arxiv.org/pdf/2105.14608.pdf

	@article{DBLP:journals/corr/abs-2105-14608,
 	 author    = {Hassan Almubarak and
               Kyle Stachowicz and
               Nader Sadegh and
               Evangelos A. Theodorou},
	  title     = {Safety Embedded Differential Dynamic Programming using Discrete Barrier
               States},
  	journal   = {CoRR},
	  volume    = {abs/2105.14608},
	  year      = {2021},
	  url       = {https://arxiv.org/abs/2105.14608},
	  eprinttype = {arXiv},
 	 eprint    = {2105.14608},
	  timestamp = {Wed, 02 Jun 2021 11:46:42 +0200},
	  biburl    = {https://dblp.org/rec/journals/corr/abs-2105-14608.bib},
	  bibsource = {dblp computer science bibliography, https://dblp.org}
	}