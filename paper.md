# Introduction to Monte Carlo methods for an Ising Model of a Ferromagnet

'If God has made the world a perfect mechanism, He has at least conceded so much to our imperfect intellects that in order to predict little parts of it, we need not solve innumerable differential equations, but can use dice with fair success.'

Max Born

Jacques Kotze

# Contents

# 1 Theory 2

1.1 Introduction . 2   
1.2 Background 2   
1.3 Model 2   
1.4 Computational Problems 4   
1.5 Sampling and Averaging 5   
1.6 Monte Carlo Method . 8   
1.7 Calculation of observables 8   
1.8 Metropolis Algorithm 9

# 2 Results 11

2.1 Energy Results . . . 11   
2.2 Magnetization Results 12   
2.3 Exact Results . 18   
2.4 Finite size scaling . 19   
2.5 Conclusion 23

A Monte Carlo source code

26

# Acknowledgements

Thank you to Prof B.A. Bassett for encouraging me to submit this report and my office colleagues in Room 402 for their assistance in reviewing it. A special thanks must be extended to Dr A.J. van Biljon who was always available for invaluable discussion, advice and assistance. Dr L. Anton and Dr K. MüllerNederbok also helped in sharing their expertise with regards this report. Lastly I would like to thank Prof H.B. Geyer for his patience and constructive criticism.

# 1 Theory

# 1.1 Introduction

This discussion serves as an introduction to the use of Monte Carlo simulations as a useful way to evaluate the observables of a ferromagnet. Key background is given about the relevance and effectiveness of this stochastic approach and in particular the applicability of the Metropolis-Hastings algorithm. Importantly the potentially devastating effects of spontaneous magnetization are highlighted and a means to avert this is examined.

An Ising model is introduced and used to investigate the properties of a two dimensional ferromagnet with respect to its magnetization and energy at varying temperatures. The observables are calculated and a phase transition at a critical temperature is also illustrated and evaluated. Lastly a finite size scaling analysis is underatken to determine the critical exponents and the Curie temperature is calculated using a ratio of cumulants with differing lattice sizes. The results obtained from the simulation are compared to exact calculations to endorse the validity of this numerical process. A copy of the code used, written in C++, is enclosed and is freely available for use and modification under the General Public License.

# 1.2 Background

In most ordinary materials the associated magnetic dipoles of the atoms have a random orientation. In effect this non-specific distribution results in no overall macroscopic magnetic moment. However in certain cases, such as iron, a magnetic moment is produced as a result of a preferred alignment of the atomic spins.

This phenomenon is based on two fundamental principles, namely energy minimization and entropy maximization. These are competing principles and are important in moderating the overall effect. Temperature is the mediator between these opposing elements and ultimately determines which will be more dominant.

The relative importance of the energy minimization and entropy maximization is governed in nature by a specific probability

$$
P ( \alpha ) = \exp { \frac { - E ( \alpha ) } { k T } } .
$$

which is illustrated in Figure 1 and is known as the Gibbs distribution.

# 1.3 Model

A key ingredient in the understanding of this theory is the spin and associated magnetic moment. Knowing that spin is a quantum mechanical phenomenon it is easy to anticipate that a complete and thorough exposition of the problem would require quantum rules of spin and angular momentum. These factors prove to be unnecessary complications.

We thus introduce a model to attain useful results. The idea central to a model is to simplify the complexity of the problem to such a degree that it is mathematically tractable to deal with while retaining the essential physics of the system. The Ising Model does this very effectively and even allows for a good conceptual understanding.

![](images/f70f2c5cf764acb693401e79faf5d329bf42927140479a038c812d01cbabec6e.jpg)  
Figure 1: This figure shows the Boltzmann probability distribution as a landscape for varying Energy ( $E$ ) and Temperature $( T )$ .

$$
\begin{array} { c c c c c } { { } } & { { } } & { { 1 } } & { { 1 } } & { { 1 } } \\ { { } } & { { } } & { { } } & { { } } & { { } } \\ { { \uparrow } } & { { \downarrow } } & { { \uparrow } } & { { \uparrow } } & { { } } \\ { { } } & { { \downarrow } } & { { \downarrow } } & { { \downarrow } } & { { } } \\ { { } } & { { \uparrow } } & { { \downarrow } } & { { \downarrow } } & { { } } \\ { { \uparrow } } & { { \downarrow } } & { { \uparrow } } & { { \downarrow } } & { { } } \end{array}
$$

Figure 2: Two Dimensional lattice illustration of an Ising Model. The up and down arrows represent a postive and a negative spin respectively.

The Ising Model considers the problem in two dimensions1 and places dipole spins at regular lattice points while restricting their spin axis to be either up ( $+$ y) or down (-y). The lattice configuration is square with dimensions $L$ and the total number of spins equal to $N = L \times L$ . In its simplest form the interaction range amongst the dipoles is restricted to immediately adjacent sites (nearest neighbours). This produces a Hamiltonian for a specific spin site, $i$ , of the form

$$
H _ { i } = - J \sum _ { j _ { n n } } s _ { i } s _ { j }
$$

where the sum $j _ { n n }$ runs over the nearest neighbours of $i$ . The coupling constant between nearest neighbours is represented by $J$ while the $s _ { i }$ and $s _ { j }$ are the respective nearest neighbour spins. The nature of the interaction in the model is all contained in the sign of the interaction coupling constant $J$ . If $J$ is positive it would mean that the material has a ferromagnetic nature (parallel alignment) while a negative sign would imply that the material is antiferromagentic (favours anti-parallel alignment). $J$ will be taken to be $+ 1$ in our discussion and the values for spins will be $+ 1$ for spin up and -1 for spin down. A further simplification is made in that $J / k _ { b }$ is taken to be unity. The relative positioning of nearest neighbours of spins is shown in Figure 3 with the darker dot being interacted on by its surrounding neighbours.

![](images/8d5aa9436305af1f2493984fd46e2419ac78613a5edac031a8c33ebe4e472485.jpg)  
Figure 3: Nearest neighbour coupling. The dark dot, at position (x,y), is being interacted upon by its nearest neighbours which are one lattice spacing away from it.

To maximize the interaction of the spins at the edges of the lattice they are made to interact with the spins at the geometric opposite edges of the lattice. This is referred to as periodic boundary condition (pbc) and can be visualized better if we consider the 2d lattice being folded into a 3d torus with spins being on the surface of this topological structure.

![](images/ac01c37663b41480006171dc4753662b9f42d186c7bc2e00d420ee20cc3ed8e6.jpg)  
Figure 4: An illustration of a three dimensional torus which is repreentative of a two dimensional lattice with periodic boundary conditions.

# 1.4 Computational Problems

With the help of an Ising Model we can proceed with the anticipation of achieving solutions for observables. If the energy of each possible state of the system is specified, then the Boltzmann distribution function, equation (1), gives the probability for the system to be in each possible state (at a given temperature) and therefore macroscopic quantities of interest can be calculated by doing probability summing. This can be illustrated by using magnetization and energy as examples. For any fixed state, $\alpha$ , the magnetization is proportional to the 'excess' number of spins pointing up or down while the energy is given by the Hamiltonian (2).

$$
M ( \alpha ) = N _ { u p } ( \alpha ) - N _ { d o w n } ( \alpha )
$$

The expected value for $M$ is given by

$$
\left. M \right. = \sum _ { \alpha } M ( \alpha ) P ( \alpha ) ,
$$

and the expected value for $E$ is given by

$$
\langle E \rangle = \sum _ { \alpha } E ( \alpha ) P ( \alpha ) .
$$

These calculation pose a drastic problem from a practical perspective. Considering we have two spin orientations (up $\&$ down) and there are $N$ spins which implies that there are $2 ^ { N }$ different states. As $N$ becomes large it is evident that computations become a daunting task if calculated in this manner.

It may seem a natural suggestion to use a computer simulation to do these calculations but by examining equations (4) and (5) more closely it becomes apparent that using this procedure would waste as much computing effort on calculating an improbable result as it does on a very probable result. Thus a better numerical alternative would be to use a simulation to generate data over the 'representative states'. These representative states constitute the appropriate proportions of different states 2. This is a form of biased sampling which essentially boils down to satisfying the following condition

# GENERATED FREQUENCY=ACTUAL PROBABILITY. (computer) (theory)

We now examine, in a more formal setting, how to accomplish this objective.

# 1.5 Sampling and Averaging

The thermal average for an observable $A ( x )$ is defined in the canonical ensemble

$$
\langle A ( x ) \rangle _ { T } = \frac { 1 } { Z } \int e ^ { - \beta H ( x ) } A ( x ) d x
$$

where $x$ is a vector in the phase space and $\beta = 1 / k _ { b } T$ . The Partition function, $Z$ , is given by

$$
Z = \int e ^ { - \beta H ( x ) } d x
$$

while the normalized Boltzmann factor is

$$
P ( x ) = { \frac { 1 } { Z } } e ^ { - \beta H ( x ) } .
$$

This probability gives the actual statistical weight with which the configuration $x$ occurs in the thermal equilibrium. We now want to consider the discrete case of the formal definitions above. If we are to consider a finite portion of the phase space it would produces an average of the form

$$
\langle A ( x ) \rangle = \frac { \sum _ { l = 1 } ^ { M } e ^ { - \beta H ( x _ { l } ) } A ( x _ { l } ) } { \sum _ { l = 1 } ^ { M } e ^ { - \beta H ( x _ { l } ) } } .
$$

If we were to take $M \to \infty$ in equation (8) it would reduce to equation (6). The problem with taking a simple sample of this nature in the phase space is that it would not guarantee that the probability distribution is peaked in the region under consideration (not representative). Figure 5 illustrates this problem.

![](images/fb1725b56d6cbdf22ab46b74aabed0f0419fcf613785e0c1c2b1d45c08b1c727.jpg)  
Figure 5: Example of a simple sampling producing a Gaussian distribution, centered around zero, while the crucial data is peaked outside the sampling region.

It thus makes sense to attempt a smarter sampling technique to include the important areas in the phase space. We want a process that selects points, $x _ { l }$ , with an associated probability, $P ( x _ { l } )$ in the phase space. Estimating the thermal average now for a chosen set, $x _ { l }$ , reduces equation (8) to

$$
\langle A ( x ) \rangle = \frac { \sum _ { l = 1 } ^ { M } e ^ { - H ( x _ { l } ) \beta } A ( x _ { l } ) / P ( x _ { l } ) } { \sum _ { l = 1 } ^ { M } e ^ { - H ( x _ { l } ) \beta } / P ( x _ { l } ) } .
$$

The most sensible choice for $P ( x _ { l } )$ is $P ( x _ { l } ) \propto e ^ { - H ( x _ { l } ) \beta }$ .This construct produces the simple arithmetic average for (9) by canceling out the Boltzmann factors, thus

$$
\langle A ( x ) \rangle = \frac { 1 } { M } \sum _ { l = 1 } ^ { M } A ( x _ { l } ) .
$$

If we stop to reflect on what we are trying to achieve at this point, we discover that we are attempting to reduce a probability distribution at equilibrium of the infinite phase space to a representative distribution with a finite set of points from the phase space, $x _ { l }$ . The question is how do we generate this distribution.

Metropolis et al. advanced conventional thinking at the time by introducing the idea of using a Markov process of successive states $x _ { l }$ where each state $x _ { l + 1 }$ is constructed from a previous state $x _ { l }$ via a suitable transition probability $W ( x _ { l } \to x _ { l + 1 } )$ . In order to implement this idea successfully a detailed balance equation has to be imposed,

$$
P _ { e q } ( x _ { l } ) W ( x _ { l } \to x _ { l ^ { \prime } } ) = P _ { e q } ( x _ { l ^ { \prime } } ) W ( x _ { l ^ { \prime } } \to x _ { l } ) .
$$

This equation simply prescribes that at equilibrium there is an equal probability for $x _ { l } \to x _ { l ^ { \prime } }$ and $x _ { l ^ { \prime } } \to x _ { l }$ . If we now take the ratio of the transition probabilities it becomes evident that a move $x _ { l } \to x _ { l ^ { \prime } }$ and the inverse move $x _ { l ^ { \prime } } \to x _ { l }$ is only dependent on the energy change $\delta H = H ( x _ { l ^ { \prime } } ) - H ( x _ { l } )$ .

$$
\frac { W ( x _ { l }  x _ { l ^ { \prime } } ) } { W ( x _ { l ^ { \prime } }  x _ { l } ) } = e ^ { - \delta H \beta }
$$

This doesn't however help to specify $W ( x _ { l } \to x _ { l ^ { \prime } } )$ uniquely. In order to do this we introduce

$$
W ( x _ { l } \to x _ { l ^ { \prime } } ) = \left\{ \begin{array} { l l } { { e ^ { - \delta H \beta } } } & { { \mathrm { i f } \delta H < 0 , } } \\ { { } } & { { } } \\ { { 1 } } & { { \mathrm { o t h e r w i s e } , } } \end{array} \right.
$$

It can be shown that using this transition probability $W ( x _ { l } \to x _ { l + 1 } )$ the distribution $P ( x _ { l } )$ of states generated by the Markov process tend to the equilibrium distribution as $M \to \infty$ . Thus the construct holds and approximates the theory with an increasing degree of accuracy as we consider a larger set of points, $\{ x _ { l } \}$ , in the phase space.

The changes in the probability distribution over time are governed by the Markovian Master equation

$$
\frac { d P ( x , t ) } { d t } = - \sum _ { x ^ { \prime } } W ( x _ { l }  x _ { l ^ { \prime } } ) + \sum _ { x ^ { \prime } } W ( x _ { l ^ { \prime } }  x _ { l } ) .
$$

In the thermal limit where $P ( x _ { l } ) = P _ { e q } ( x _ { l } )$ the detailed balance equation, (11), which was imposed, comes into play resulting in ${ d P ( x , t ) } / { d t } = 0$ as we would expect. Furthermore since we are considering a finite system it is logical to conclude that the system is ergodic3.

$$
\left. A ( t ) \right. = \frac { 1 } { t _ { M } } \int A ( t ) d t
$$

This relation reduces in a similar fashion to the arithmetic average previously discussed if we consider the number of Monte Carlo steps (mcs) as a units of 'time'. We thus are confronted by the question of whether the system is

ergodic in order for the time averaging to be equivalent to the canonical ensemble average. This condition is thus forced upon the system if we consider the mcs as a measure of time.

$$
\langle A ( t ) \rangle = \frac { 1 } { M } \sum _ { t = 1 } ^ { M } A ( x ( t ) )
$$

Thus the Metropolis sampling can be interpreted as time averaging along a stochastic trajectory in phase space, controlled by the master equation (14) of the system.

# 1.6 Monte Carlo Method

The fact that we can make a successful stochastic interpretation of the sampling procedure proves to be pivotal in allowing us to introduce the Monte Carlo method as a technique to solve for the observables of the Ising Model.

A Monte Carlo calculation is defined, fundamentally, as explicitly using random variates and a stochastic process is characterized by a random process that develops with time. The Monte Carlo method thus lends itself very naturally to simulating systems in which stochastic processes occur. From what we have established with regards to the sampling being analogous to a time averaging along a stochastic trajectory in the phase space it is possible to simulated this process by using the Monte Carlo method. The algorithm used is design around the principle of the Metropolis sampling.

From an implementation point of view the issue of the random number lies at the heart of this process and its success depends on the fact that the generated number is truly random. 'Numerical Recipes in $\mathrm { C }$ , [9], deals with this issue extensively and the random number generator 'Ran1.c' from this text was used in the proceeding simulation. This association with randomness is also the origin of the name of the simulation since the glamorous location of Monte Carlo is synonymous with luck and chance.

# 1.7 Calculation of observables

The observables of particular interest are $\langle E \rangle$ , $\langle E ^ { 2 } \rangle , \langle M \rangle , \langle | M | \rangle$ and $\langle M ^ { 2 } \rangle$ These are calculated in the following way,

$$
\left. M \right. = \frac { 1 } { N } \sum _ { \alpha } ^ { N } M ( \alpha )
$$

similarly $\langle | M | \rangle$ and $\langle M ^ { 2 } \rangle$ are calculated using the above equation.

To calculate energy we use the Hamiltonian given in equation (2).

$$
\langle E \rangle = { \frac { 1 } { 2 } } \langle \sum _ { i } ^ { N } H _ { i } \ \rangle = { \frac { 1 } { 2 } } \langle \ - J \sum _ { i } ^ { N } \sum _ { j _ { n n } } s _ { i } s _ { j } \ \rangle
$$

the factor of a half is introduced in order to account for the spins being counted twice. Equation (18) is used in a similar way to determine $\left. E ^ { 2 } \right.$ .

At the Curie temperature we expect a marked fuctuation in these quantities. A good candidate to illustrate this fluctuation would be the variance

$( \Delta A ) ^ { 2 } = \langle A ^ { 2 } \rangle - \langle A \rangle ^ { 2 }$ . This leads us to the logical conclusion of calculating the heat capacity, $C$ , and the susceptibility, $\chi$ .

$$
C = { \frac { \partial E } { \partial T } } = { \frac { ( \Delta E ) ^ { 2 } } { k _ { b } T } } = { \frac { \langle E ^ { 2 } \rangle - \langle E \rangle ^ { 2 } } { k _ { b } T ^ { 2 } } }
$$

$$
\chi = { \frac { \partial M } { \partial T } } = { \frac { ( \Delta M ) ^ { 2 } } { k _ { b } T } } = { \frac { \langle M ^ { 2 } \rangle - \langle M \rangle ^ { 2 } } { k _ { b } T } }
$$

A cumulant is also calculated. This will be used to ultimately determine the Curie temperature.

$$
U _ { L } = 1 - \frac { \langle M ^ { 4 } \rangle _ { L } } { 3 \langle M ^ { 2 } \rangle _ { L } }
$$

# 1.8 Metropolis Algorithm

The algorithm implemented for this simulation is the Metropolis Algorithm [8]. The steps executed in the program are best summarized in a fowchart. From the flowchart, Figure 6, it is possible to attain a better conceptual feel for what the algorithm is attempting to achieve.

![](images/e759680c0524ada1042a8fd44049d65d4b1ba1e7b91c3e82e8d1d13b8a98641a.jpg)  
Figure 6: Metropolis Flowchart

•In the first step the lattice is INITIALIZED to a starting configuration. This may either be a homogeneous or random configuration. A random configuration has a benefit in that it uses less computing time to reach a equilibrated configuration with the associated heat bath.

• In the following PROCESS the random number generator is used to select a position on the lattice by producing a uniformly generated number between 1 and $N$ .

•A DECISION is then made whether the change in energy of flipping that particular spin selected is lower than zero. This is in accordance with the principle of energy minimization.

If the change in energy is lower than zero then a PROCESS is invoked to flip the spin at the selected site and the associated change in the observables that want to be monitored are stored.

If the change in energy is higher than zero then a DECISION has to be used to establish if the spin is going to be fipped, regardless of the higher energy consideration. A random number is generated between 0 and 1 and then weighed against the Boltzmann Probability factor. If the random number is less than the associated probability, $e ^ { - \delta \beta H }$ , then the spin is flipped (This would allow for the spin to be flipped as a result of energy absorbed from the heat bath, as in keeping with the principle of entropy maximization) else it is left unchanged in its original configuration.

• The above steps are repeated N times and checked at this point in a DECISION to determine if the loop is completed. The steps referred to here do not include the initialization which is only required once in the beginning of the algorithm.

• Once the N steps are completed a PROCESS is used to add all the progressive changes in the lattice configuration together with the original configuration in order to produce a new lattice configuration.

• All these steps are, in turn, contained within a Monte Carlo loop. A DECISION is used to see if these steps are completed.

•Once the Monte Carlo loop is completed the program is left with, what amounts to, the sum of all the generated lattices within the N loops. A PROCESS is thus employed to average the accumulated change in observables over the number of spins and the number of Monte Carlo steps.

Lastly this data can be OUTPUT to a file or plot.

This run through the algorithm produces a set of observables for a specific temperature. Considering that we are interested in seeing a phase transition with respect to temperature we need to contain this procedure within a temperature loop in order to produce these observables for a range of temperatures.

The program that was implemented for this discussion started at a temperature of $T = 5$ and progressively stepped down in temperature to $T = ~ 0 . 5$ with intervals of $\delta T = 0 . 1$ . The different lattice sizes considered where $2 \times 2 , ~ 4 \times 4$ ,

$8 \times 8$ and $1 6 \times 1 6$ . A million Monte Carlo steps (mcs) per spin where used in order to ensure a large sample of data to average over. The physical variables calculated were $\langle E \rangle$ , $\left. E ^ { 2 } \right.$ , $\langle \left| M \right| \rangle$ and $\langle M ^ { 2 } \rangle$ .

A problem that occurs after the initialization is that the configuration will, more than likely, not be in equilibrium with the heat bath and it will take a few Monte Carlo steps to reach a equilibrated state. The results produced during this period are called transient and aren't of interest. We thus have to make provision for the program to disregard them. This is achieved by doing a run of a thousand Monte Carlo steps preceding the data collection of any statistics for a given temperature in order to ensure that it has reached a equilibrated state. This realization is only significant for the initial configuration and the problem is avoided by the program in the future by using small temperature steps and the configuration of the lattice at the previous temperature. A very small number of mcs are thus required for the system to stabilize its configuration to the new temperature.

# 2 Results

# 2.1 Energy Results

![](images/7177b3835bfaf43233e69e2464688d21e3f31f33baaad4b1cf92dfc6ec167ec7.jpg)  
Figure 7: This plot shows the differing results of the Energy for varying lattice sizes, $L \times L$ .

In Figure 7 the energy per spin as a function of temperature can be seen. The curve of the graph becomes more pronounced as the lattice size increases but there isn't a marked difference between the $L = 8$ and $L = 1 6$ lattices. The steep gradient in the larger lattices points towards a possible phase transition but isn't clearly illustrated. The energy per spin for higher temperatures is relatively high which is in keeping with our expectation of having a random configuration while it stabilizes to a $E / N = - 2 J = - 2$ at low temperatures. This indicates that the spins are all aligned in parallel.

![](images/8be43f57cba6d3011514095a31d08de21849febf3973f6e928f9130a3dec07bf.jpg)  
Figure 8: This plot shows the differing results of Specific Heat Capacity for varying lattice sizes, $L \times L$ .

In Figure 8 the specific heat capacity per spin is shown as a function of temperature. We concluded previously that a divergence would occur at a phase transition and thus should be looking for such a divergence on the graph. It is however clear that there is no such divergence but merely a progressive steepening of the peak as the lattice size increases. The point at which the plot is peaked should be noted as a possible point of divergence. The reason for not explicitly finding a divergence will be discussed in Section 2.4.

# 2.2 Magnetization Results

Figure 9 of the magnetization results shows very beautifully that the shape of the gradient becomes more distinct as the lattice size is increased. Furthermore, as opposed to Figure 7, there is a far more apparent difference that the larger lattices produce in the curves and this illustrates a more apparent continuous phase transition. The behaviour of the magnetization at high and low temperature are as the theory prescribes (random to stable parallel aligned configuration).

At this juncture it is prudent to point out that the susceptibility cannot be calculated using the ordinary technique in equation (20) given in the discussion on the calculation of observables. The reason is focused around a subtle fact that has drastic implications. To comprehend the problem at work we have to consider one of the constraints of our model, namely the finite nature of our lattice. This manifests in the fact that spontaneous magnetization can occur for a finite sized lattice. In this instance the effect is of particular interest below the critical temperature.

This can be illustrated by considering the following example of collected data in Figure 10. This data is taken at a temperature that is considerably less than the Curie temperature and we would thus expect it to have a stable nature and yet it clearly displays a fluctuation that is uncharacteristic, resulting in a complete flip of the magnetization. It has already been highlighted that because we are dealing with a limited lattice size there is finite probability for this kind of behaviour to take place. This probability is directly proportional to the number of mcs used and inversely proportional to the lattice size, this is compounded by the preiodic boundary conditions used.

![](images/96ab93290aadd51ec373b6b0793856449dd9f83d5d41d4600ab064a3df6436fb.jpg)  
Figure 9: This plot shows the differing results of the Magnetization for varying lattice sizes, $L \times L$ .

Figure 11 schematically depicts this fact. The valley shown linking the two peaks of the probability will thus be dropped for lower temperatures and bigger lattice configurations. It should be noted that even though the probability may be less it does always exist and this has to be accounted for in data collection or it may corrupt the results. We expect the configuration to be relatively stable at the peaks but if its magnetization has slipped down (fluctuated) to the center of the valley then it has an equal probability of climbing up either side of the peaks, this is the crux of the spontaneous flipping. This aspect of symmetry proves to also be the seed for a possible solution to this problem.

As an example of what has just been mentioned we note from Figure 10 where a fluctuation occurs just before 5000 mcs and the magnetization peaks at 0 from $^ { - 1 }$ . The configuration is now in the middle of the valley and happens to go back to its previous state. The same phenomenon occurs just after 5000 mcs but in this instance chooses to fip to an opposite, but equally probable, magnetization, from -1 to 1.

If we now were to think of the implications of this spontaneous flipping we come to the realization that it would cause an averaging out of the mean magnetization, $\langle M \rangle$ . This of course has a detrimental effect on calculating the variance of the magnetization and thus the susceptibility.

This can be illustrated in the Figure 12 where the plot shows that $\langle M \rangle ^ { 2 }$ remains zero for an extended period at low temperatures. This would cause the variance to peak at lower temperatures. As the lattice size increases the spontaneous magnetization is less likely to occur and the critical point moves progressively to higher temperatures, this implies that the peak for the susceptibility would approach the Curie temperature from the left (lower temperatures). This is inconsistent with what the theory prescribes [3]. We can also conclude that the greater the number of mcs we use the more likely we are to introduce spontaneous flipping and thus averaging out of the mean magnetization which would move the peak of the susceptibility more to the left.

![](images/c5ebf307a986bc6f9ac7700341761cc811d71db3d3a197dc9ad8238b955f7a67.jpg)  
Figure 10: This plot shows a spontaneous flip in magnetization for a $2 \times 2$ lattice at $T = 1$ .

The solution to this problem lies in a characteristic that is at the heart of the problem, namely that there is an equal probability for the magnetization to change to an opposite configuration or go back to its previous configuration. Thus if we were to use the absolute magnetization we would effectively be considering only one probability peak (positive one) and the complications of averaging out the mean magnetization would be overcome. Figure 11 thus changes to a distribution shown in Figure 13 if we where to use the absolute magnetization.

This modification doesn't come without a cost. In this instance it can be seen that we have averted the problems at low temperatures but end up with weaker statistics at high temperatures. This is apparent from the fact that previously we had frequent fluctuations, in Figure 11, between positive and negative magnetization at high temperatures resulting in a zero mean magnetization. However we have only positive values to consider for the averaging of the mean magnetization producing a non zero average. This effect is reduced slightly since the valley is raised at high temperatures resulting in the magnetization having a high probability of being close to zero. Fortunately this nonzero average for the magnetization at higher temperatures is inconsequential since it doesn't influence the Curie temperature and appears only in the region above it.

At lower temperatures the shape of the distribution changes, as indicated by the dotted line in Figure 13. Thus the magnetization remains relatively stable in a homogeneous configuration. This is exactly the behaviour we expect and produces good results.

![](images/67fc46b6e6cf144c436f0e56bf0b2c627c4c12403bdbd4cced55b9dcf37781dd.jpg)  
Figure 11: A schematic illustration of the probability density of the magnetization and how the representative spin distribution would populate a square latice of size $L$ . The darker and lighter regions depict negative and positive spin repsectively.

The susceptibility that is produced in our data is thus not exactly equivalent to the theoretical susceptibility, $\chi$ , and we will be distinguished as $\chi ^ { \prime }$ . The scaling characteristic of this susceptibility is, however, equivalent to the theoretical value and only varies by a constant factor above the Curie temperature.

$$
\chi ^ { \prime } = \frac { \langle M ^ { 2 } \rangle - \langle | M | \rangle ^ { 2 } } { k _ { b } T }
$$

A comparison can be made between $\chi ^ { \prime }$ and $\chi$ in Figures 14 and 15 respectively. It is clear that a marked difference in results occurs. This mistake becomes even more evident if you were to use $\chi$ to get the critical exponent using finite size scaling. Only $\chi ^ { \prime }$ produces the correct finite size scaling.

We now evaluate the plots produced using this technique discussed thus far. The more distinctive character of the differing plots of Figure 9 produce more dramatic peaks in Figure 14 of the magnetic susceptibility ( $\chi ^ { \prime }$ ) per spin versus temperature as a result. This, once again, doesn't show an exact divergence but shows a sharp peak for the $L = 1 6$ lattice. This should be strong evidence eluding to a second order phase transition.

![](images/9bf2d5ec8116af7d68f090a0e9c59ffe2f7a30442bdca3e0e29fbe226fd48599.jpg)  
Figure 12: This is an illustration of the differences in the normalized values of $\langle M \rangle ^ { 2 } ; \langle | M | \rangle ^ { 2 }$ and $\langle M ^ { 2 } \rangle$ with respect to temperature. The varying lattice sizes considerd are $2 x 2$ (top left); $4 x 4$ top right; $8 x 8$ (bottom left) and $1 6 x 1 6$ (bottom right).

![](images/8c0737827e599a040491e877e2d7e5d66f2aae6fbf40439ae4f5846e926fd8d2.jpg)  
Figure 13: The solid line shows the revised probability density when using the absolute magnetization as opposed to the dotted line which represents the orginal propability density for magnetization.

![](images/eab34f22bd052bb322aab7267adfb95550c4952fb39f5326eed6557fa7388fa1.jpg)  
Figure 14: This plot shows the differing results of the susceptibility for varying lattice sizes, $L \times L$ .

![](images/250237a61c33e19b5b80176d58d6c726211b398fd840b0385de237fb8e0d9ca7.jpg)  
Figure 15: This plot shows the differing results of the susceptibility for varying lattice sizes, $L \times L$ .

# 2.3 Exact Results

As pointed out in the beginning of this work there are significant challenges to the calculation of exact solutions but we do need to evaluate the dependability of the numerical process used in this discussion. For this purpose we restrict the comparison between simulation and exact results to only a $2 \times 2$ This only serves as a good first order approximation but the corroboartion should imporve as the lattice size is increased.

The $2 ^ { 4 }$ different configurations that can be listed can be reduced by symmetric considerations to only four configurations. Configuration A is the fully aligned spin configuration and configuration B has one spin in the opposite direction to the other three while C and D have two spins up and two spins down. To generate the full sixteen different configurations we need only consider geometric variations of these four configurations.

![](images/6455299ed5728c3f2a458ed7198577616d8abaae35cb28b78b237e7b3b11fbb4.jpg)  
Figure 16: The four significant lattice configurations for a $2 \times 2$ lattice.

Using equations (18) and (17) we can calculate the energy and magnetization for each of the four configurations. Taking into account the degeneracy of each of these configurations allows us to generate the exact equations for calculating the relevant observables. These results are listed in Table 1.

Table 1: Energy and Magnetization for respective configurations illustrated in Figure 16.   

<table><tr><td>Configuration</td><td>Degeneracy</td><td>Energy(E)</td><td>Magnetization(M)</td></tr><tr><td>A</td><td>2</td><td>-8</td><td>+4,-4</td></tr><tr><td>B</td><td>8</td><td>0</td><td>+2,-2</td></tr><tr><td>C</td><td>4</td><td>0</td><td>0</td></tr><tr><td>D</td><td>2</td><td>8</td><td>0</td></tr></table>

Using equations (4) and (5) in conjunction with the results of Table 1 produces the following equations

$$
Z = { 2 \ e ^ { 8 \beta J } + 1 2 + 2 \ e ^ { - 8 \beta } }
$$

$$
\langle E \rangle = - { \frac { 1 } { Z } } [ ~ 2 ( 8 ) ~ e ^ { 8 \beta } + 2 ( - 8 ) ~ e ^ { - 8 \beta } ~ ]
$$

$$
\langle E ^ { 2 } \rangle = \frac { 1 } { Z } [ \ 2 ( 6 4 ) \ e ^ { 8 \beta } + 2 ( 6 4 ) \ e ^ { - 8 \beta } \ ]
$$

$$
\langle | M | \rangle = \frac { 1 } { Z } [ \ 2 ( 4 ) \ e ^ { 8 \beta } + 8 ( 2 ) \ ]
$$

$$
\langle M ^ { 2 } \rangle = \frac { 1 } { Z } [ \ 2 ( 1 6 ) \ e ^ { 8 \beta } + 8 ( 4 ) \ ] .
$$

Applying equations (24); (25); (26) and (27) to the formulas given in (19) and (22) we can obtain the heat capacity and susceptibility respectively. The comparison of these results are shown in Figure 17. The results achieved by the Monte Carlo method match the exact calculation exceptionally well and only deviate very slightly from the prescribed heat capacity at high temperatures.

![](images/a65ff4b22ffe523a50355c39e9bf9f60c2572c3bda340b30be69c22bb0a8042c.jpg)  
Figure 17: This plot shows a favourable comparison between the exact calculations done for the Heat Capacity per spin $( C / N )$ and the Susceptibility per spin $( \chi ^ { \prime } / N )$ with the numerically generated solutions of the Monte Carlo simulation.

# 2.4 Finite size scaling

One of the limitations that the Ising Model confronts us with is the finite size of our lattice. This results in a problem of recognizing the specific point at which the phase transition occurs. This should be at a theoretical point of divergence but we are limited by the size of the lattice under consideration and thus dont see this divergence. This effect is minimized by using periodic boundary conditions but would only be resolved if we where to consider an infinitely sized lattice as with the associated theoretical values for the phase transition. It is thus necessary to use a construct that will allow us to extrapolate the respective theoretical value given the limited resource of a finite sized lattice. The aptly named procedure of finite size scaling is used to do just this.

It becomes useful to define a critical exponent to better understand the nature of the divergence near the critical temperature. The critical exponent, $\lambda$ , is given by $\begin{array} { r } { \lambda = \operatorname* { l i m } _ { t  0 } \ \frac { \ln | F ( t ) | } { \ln | t | } } \end{array}$ or more commonly written as $F ( t ) \sim | t | ^ { \lambda }$ where $t = ( T - T _ { c } )$ . This exponent is important in the sense that it offers a more universal characteristic for differing data collected. This attribute will be taken advantage of and illustrated by showing that in a reduced unit plot the data collapses to a common curve with the same critical exponent.

The critical exponents relevant to the Ising model are as follows:

$$
\begin{array} { r } { \xi ( T ) \sim | T - T _ { c } | ^ { - \nu } } \\ { M ( T ) \sim ( T _ { c } - T ) ^ { \beta } } \\ { C \sim | T - T _ { c } | ^ { - \alpha } } \\ { \chi \sim | T - T _ { c } | ^ { - \gamma } } \end{array}
$$

If we examine the relationship between the lattice size, $L$ , with respect to the temperature relationship, $| T - T _ { c } |$ , we discover that $| T - T _ { c } | < < 1$ as $L  \infty$ . Thus a critical exponent is also applicable for the lattice size. This produces $L \sim | T _ { c } ( L = \infty ) - T _ { c } ( L ) | ^ { - \nu }$ which can in turn be used to reduce the above exponents for the Ising model to a more appropriate form, in terms of lattice size.

$$
\begin{array} { r } { \xi ( T ) \sim | T - T _ { c } | ^ { - \nu }  L } \\ { M ( T ) \sim ( T _ { c } - T ) ^ { \beta }  L ^ { - \beta / \nu } } \\ { C ( T ) \sim | T - T _ { c } | ^ { - \alpha }  L ^ { \alpha / \nu } } \\ { \chi \sim | T - T _ { c } | ^ { - \gamma }  L ^ { \gamma / \nu } } \end{array}
$$

This proportionality is illustrated in a schematic form in Figure 18 of the relationship between the plot of a certain observable and the respective critical exponent. In this instance susceptibility was used as an example.

![](images/40ed5049c075194e871282ba9c0ba1a0acc8d63d02bdad5df8010a5870dc02f4.jpg)  
Figure 18: Theoretical llustrationof Finite Size Scaling relationship for susceptibility $\chi$ .

We can now attempt to use equations (28) - (31) and determine their appropriate exponents. This is simply done by taking the peak values for the collected data of the observables and plotting a ln-ln graph that should yield a straight line with the gradient being equal to the respective critical exponents.

This procedure is made easier since $\nu$ is equal to 1 for a two dimensional lattice. An additional calculation for the observables of a $L = 3 2$ lattice were done, since this would increase the number of data points. This would also offset the poor statistics associated with the $L = 2$ lattice and ultimately allow for a more accurate result.

![](images/6663fef3c726ba1770a452f9686c0a064aead1990f59a93d22197f03806e2bdc.jpg)  
Figure 19: Plots obtained for the calculation of the critical exponents.

In Figure 19 the calculation for the critical exponent is shown. The answers are listed in Table 2. The graph for the heat capacity isn't a straight line and shows a curvature. The reason is that $\alpha$ is zero in the two dimensional Ising model and should rather be interpreted as $C \sim C _ { 0 }$ ln $L$ .

Table 2: Displays the calculated and theoretical critical exponents.   

<table><tr><td>Quantity</td><td>Exponent</td><td>Finite Size Scaling (2d)</td><td>Theoretical (2d)</td></tr><tr><td>Magnetization</td><td>β</td><td>0.128 ± 0.008</td><td>0.125</td></tr><tr><td>Susceptibility</td><td>γ</td><td>1.76 ± 0.01</td><td>1.75</td></tr><tr><td>Heat Capacity</td><td>C0</td><td>0.518 ± 0.02</td><td>0.500</td></tr></table>

We can in a naive sense attempt to determine the Curie temperature by implementing the critical relation for the lattice size. This is listed in Table 3. It doesn't prove to be very accurate and we thus have to implement a different notion in order to ascertain the critical temperature accurately.

The transition point can be determined by using the cumulant

$$
U _ { L } = 1 - \frac { \langle M ^ { 4 } \rangle _ { L } } { 3 \langle M ^ { 2 } \rangle _ { L } } .
$$

![](images/2773b47b9ba7ef3783912b338a9fc84e059e910a05ee5a101f27d9b3275756a4.jpg)  
Figure 20: This graph shows a reduced unit plot for $\chi ^ { \prime }$ with a critical exponent of $\gamma = 1 . 7 5$ .

Table 3: Listed critical temperatures calculated from the critical lattice size relation.   

<table><tr><td>Lattice Size (L)</td><td>Estimated Tc</td><td>Estimated Tc(L → ∞)</td><td>Theoretical Tc(L = ∞)</td></tr><tr><td>2</td><td>3.0</td><td>2.50</td><td>2.269</td></tr><tr><td>4</td><td>2.8</td><td>2.55</td><td>2.269</td></tr><tr><td>8</td><td>2.5</td><td>2.43</td><td>2.269</td></tr><tr><td>16</td><td>2.4</td><td>2.34</td><td>2.269</td></tr><tr><td>32</td><td>2.3</td><td>2.27</td><td>2.269</td></tr></table>

This calculation has to use double precision in order to retain a high level of accuracy in calculating the critical point. To determine the critical point we choose pairs of linear system sizes $( L , L ^ { \prime } )$ . The critical point is fixed at $U _ { L } = U _ { L ^ { \prime } }$ . Thus taking the ratio of different cumulants for different sized lattices, $U _ { L } / U _ { L ^ { \prime } }$ , we will get an intersection at a particular temperature. This is the desired critical temperature. This procedure is not as straightforward as it may seem and requires the cumulants to be collected very near to the transition point. Thus an iterative process may need to be employed in order to narrow down the region of where the critical temperature is located.

This analysis is done until an unique intersection is found for all cumulants. This method is illustrated in Figure 21. The $L = 2$ lattice isn't shown since it doesn't exhibit the level of accuracy that we desire. From the fitted graph of the cumulants it can be seen that intersection is common to all the cumulants except for the $U _ { 4 }$ . This points towards a poor value for this particular statistic and is thus not used, along with $U _ { 2 }$ , in the ratios of the cumulants to determine the Curie temperatures. The final analysis produces a result of a $T _ { c } = 2 . 2 6 6$ which agrees favourably with the theoretical value of $T _ { c } = 2 . 2 6 9$ .

![](images/278e1fb13941cf1d2e9a4585f81c1fd5f44cd82a32ce63385dc700ccc27e3dde.jpg)  
Figure 21: Plots of the cumulants $U _ { L }$ ).

# 2.5 Conclusion

In this review we pointed out why an exact exposition of an Ising model for a ferromagnet is not easily achieved. An argument was proposed and motivated for making the computational problem far more tractable by considering a stochastic process in combination with the Metroplis-Hastings sampling algorithm. The numerical results produced by the Monte Carlo simulation compare favourably with the theoretical results and are a viable and efficient alternative to an exact calculation. The requirements for producing accurate results are to consider large lattice sizes and a large number of Monte Carlo steps. The accuracy is very compelling even for small lattice sizes.

It is important to note and make provision for the potential of spontaneous magnetization to occur in a finite sized lattice. This can have serious consequences on the accuracy of the magnetization and the susceptibility which in turn will lead to incorrect results for the finite size scaling of these observables. The occurrence and severity of spontaneous magnetization is directly proportional to the number of Monte Carlo steps used and inversely proportional to the lattice size considered. A practical means to overcome this complication is to use the absolute magnetization in the variance of the susceptibility instead of just the magnetization. This is an effective solution that produces good statistics only deviating slightly at high temperatures from the theoretical values.

In conclusion a finite size scaling analysis was undertaken to determine the critical exponents for the observables. These where in good agreement with the theoretical exponents. The critical temperature was also calculated using a ratio of cumulants with differing lattice sizes and generated results which were in good agreement with the theoretical values.

# References

[1] Statistical Mechanics of Phase Transitions. Clarendon Press, Oxford, 1992.   
[2] K. Binder and D.W. Heermann. Monte Carlo Simulation in Statistical Physics. Springer, Berlin, 1997.   
[3] John Cardy. Scaling and Renormalization in Statistical Physics. Cambridge University Press, London, 1999. [4] D.P.Landau. Finite-size behavior of the ising square lattice. Physical Review B, 13(7):29973011, 1 April 1976.   
[5] N.J. Giordano. Computation Physics. Prentice Hall, New Jersey, 1997.   
[6] K.H. Hoffmann and M. Schreiber. Computation Physics. Springer, Berlin, 1996. [7] W. Kinzel and G. Reents. Physics by Computer. Springer, Berlin, 1998.   
[8] N. Metropolis, A.W. Rosenbluth, M.N. Rosenbluth, A.H. Teller, and E. Teller. Journal of Chemical Physics, 1953.   
[9] W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery. Numerical Recipes in C. Cambridge University Press, Cambridge, 1992.   
[10] L.E. Reichel. A Modern Course in Statistical Physics. University of Texas Press, Austin, 1980.   
[11] Robert H. Swendsen and Jian-Sheng Wang. Nonuniversal critical dynamics in monte carlo simulations. Physical Review Letters, 58(2):8688, 12 January 1987.

<table><tr><td></td><td colspan="3">A Monte Carlo source code</td></tr><tr><td colspan="2">Monte Carlo Simulation of a 2D Ising Model Page 1/2 #include &lt;iostream.h&gt;</td><td>Monte Carlo Simulation of a 2D Ising Model Page 2/2</td><td></td></tr><tr><td colspan="2">#include &lt;math.h&gt; #include &lt;stdlib.h&gt; #include &lt;fstream.h&gt;</td><td>//function for disregarding transient results transient_results() lat_type pos; int de=0;</td><td></td></tr><tr><td colspan="2">//random number generator from &quot;Numerical Recipes in C&quot; (Ranl.c) #include &lt;&quot;random.h&quot;&gt;</td><td>for(int a=1;a&lt;=transient;a++) for(int b=1;b&lt;=n;b++)</td><td></td></tr><tr><td colspan="2">//file to output data into ofstream DATA(&quot;DATA.1.dat&quot;, ios: :out);</td><td>if(test_flip(pos,de)) flip(pos);</td><td></td></tr><tr><td colspan="2">//structure for a 2d lattice with coordinates x and y struct lat_type</td><td></td><td></td></tr><tr><td colspan="2">int y; int x; const int size=2; //lattice size //array size for lattice</td><td>//function for calculating total magnetization of lattice</td><td></td></tr><tr><td colspan="2">const int lsize=size-1; const int n=size*size; //number of spin points on lattice const float minT=0.5; float T=5.0; //minimum temperature //starting point for temperature //size of steps for temperature loop</td><td>int total_magnetization() int m=0; for(int y=size;y&gt;=1;y--)</td><td></td></tr><tr><td colspan="2">float change=0.1; int lat [size+1][size+1]; lotgat signed intimes=10000; //2d lattice for spins //number of Monte Carlo steps int transient=1000; //number of transient steps //normalization for averaging //seed for random number generator</td><td>for(int x=1;x&lt;=size;x++)</td><td></td></tr><tr><td colspan="2">double norm=(1.0/float(mcs*n)); long int seed=436675;</td><td>m+=1at [x] [y]; return m;</td><td></td></tr><tr><td colspan="2">//function for random initialization of lattice</td><td></td><td></td></tr><tr><td colspan="2">initialize(int lat[size+1][size+1]) for(int y=size;y&gt;=1;y--)</td><td>//function for calculating total energy of lattice</td><td></td></tr><tr><td colspan="2">for(int x=1;x&lt;=size;x++)</td><td>int total_energy() lat_type pos; int e=0;</td><td></td></tr><tr><td colspan="2">if (ran1(&amp;seed) &gt;=0.5)</td><td>for(int y=size;y&gt;=1;y--)</td><td></td></tr><tr><td colspan="2">lat [x] [y]=1; else 1at [x] [y]=−1;</td><td>pos.y=y; for(int x=1;x&lt;=size;x++)</td><td></td></tr><tr><td colspan="2"></td><td>pos.x=x; e+=energy_pos(pos);</td><td></td></tr><tr><td colspan="2"></td><td>return e;</td><td></td></tr><tr><td colspan="2">//output of lattice configuration to the screen output(int lat[size+1] [size+1])</td><td></td><td></td></tr><tr><td colspan="2">for(int y=size;y&gt;=1;y--)</td><td>//main program void main()</td><td></td></tr><tr><td colspan="2">for(int x=1;x&lt;=size;x++)</td><td>double E=0,Esq=0, Esq_avg=0, E_avg=0, etot=0, etotsq=0; //declaring variables to be_used in calculating the observables</td><td></td></tr><tr><td colspan="2">else</td><td>doub1e M=0, Msq=0, Msq_avg=0,M_avg=0, mtot=0, mtotsq=0; double Mabs=0,Mabs_avg=0,Mq_avg=0,mabstot=0,mqtot=0; int de=0;</td><td></td></tr><tr><td colspan="2">cout&lt;&lt;&quot;+&quot;; cout&lt;&lt;endl;</td><td>lat_type pos; //initialize lattice to random configuration</td><td></td></tr><tr><td colspan="2"></td><td>initialize(lat); //Temperature loop</td><td></td></tr><tr><td colspan="2">//function for choosing random position on lattice choose_random_pos_lat(iat_type &amp;pos)</td><td>for(;T&gt;=minT;T=T-change)</td><td></td></tr><tr><td colspan="2">pos.x=(int)ceil(ranl(&amp;seed)*(size)); pos.y=(int)ceil(ranl(&amp;seed)*(size));</td><td>//transient function transient_results(); //observables adopt equilibrated lattice configurations values</td><td></td></tr><tr><td colspan="2">if(pos.x&gt;sizel|pos.y&gt;size) cout&lt;&lt;&quot;error in array size allocation for random position on lattice!&quot; ; exit;</td><td>Mabsabs(total_magnetization()); M=total_magnetization(); E=total_energy();</td><td></td></tr><tr><td colspan="2"></td><td>//initialize summation variables at each temperature step etot=0;</td><td></td></tr><tr><td colspan="2">//function for calculating energy at a particular position on lattice</td><td>etotsq=0; mot=0; mtotsq=0;</td><td></td></tr><tr><td colspan="2">int energy_pos(lat_type &amp;pos)</td><td>mabstot=0; mqtot=0;</td><td></td></tr><tr><td colspan="2">//periodic boundary conditions int up,down,left,right,e; if(pos.y==size)</td><td>//Monte Carlo loop. for(inta=l;a&lt;=mcs;a++)</td><td></td></tr><tr><td colspan="2">up=1; else up=pos.y+1;</td><td>//Metropolis loop for(int b=1;b&lt;=n;b++)</td><td></td></tr><tr><td colspan="2">if (pos.y==1) down=size;</td><td>choose_random_pos_lat (pos); if(test_flip(pos,de))</td><td></td></tr><tr><td colspan="2">else down=pos.y-1; if(ps.x==1)</td><td>flip(pos);</td><td></td></tr><tr><td colspan="2">left=size; else</td><td>//adjust observables</td><td></td></tr><tr><td colspan="2">1ef=pos-1; if(pos.x==size)</td><td>E+=2*de; M+=2*1at [pos.x] [pos,y];</td><td></td></tr><tr><td colspan="2">right=1; else</td><td>Mabs+=abs(lat [pos.x] [pos.y]);</td><td></td></tr><tr><td colspan="2">right=pos.x+1;</td><td></td><td></td></tr><tr><td colspan="2">//energy for specific position</td><td>//keep summation of observables etot+=E/2.0;//so as not to count the energy for each spin twice</td><td></td></tr><tr><td colspan="2">e=-1*1at [pos.x] [pos.y] *(lat[left] [pos.y]+lat [right] [pos.y]+lat [pos.x] [up]+1at [pos.x] [down]);</td><td>etotsq+=E/2.0*E/2.0;</td><td></td></tr><tr><td colspan="2">return e;</td><td>mtot+=M; mtotsq+=M*M;</td><td></td></tr><tr><td colspan="2"></td><td>mqtot+=M*M*M*M; mabstot+=(sqrt(M*M));</td><td></td></tr><tr><td colspan="2">//function for testing the validity of flipping a spin at a selected position</td><td></td><td></td></tr><tr><td colspan="2">bool test_flip(lat_type pos, int &amp;de)</td><td>//average observables E_avg=etot*norm;</td><td></td></tr><tr><td colspan="2">de=-2*energy_pos (pos); //change in energy for specific spin</td><td>Esq_avg=etotsq*norm; M_avg=mtot*norm;</td><td></td></tr><tr><td colspan="2">if(de&lt;0) return true; //flip due to lower energy //flip due to heat bath</td><td>Msq_avg=mtotsq*norm;</td><td></td></tr><tr><td colspan="2">else if(ran1(&amp;seed)&lt;exp(-de/T)) //no flip</td><td>Mabs_avgmabstot*norm; Mq_avg=mqtot*norm;</td><td></td></tr><tr><td colspan="2">return true; else</td><td></td><td></td></tr><tr><td colspan="2">return false;</td><td>//output data to file</td><td>//temperature</td></tr><tr><td colspan="2"></td><td>DATA&lt;&lt;T&lt;&lt;</td><td>//&lt;M&gt;;&lt;|M|&gt;;&lt;M^2&gt; per spin</td></tr><tr><td colspan="2"></td><td>&quot;t&quot;&lt;&lt;M_avg&lt;&lt;mt&quot;&lt;&lt;Mabs_avg&lt;&lt;m1tm&lt;&lt;Msq_avg&lt;</td><td>//susceptibility per spin (x)</td></tr><tr><td colspan="2"></td><td>&quot;t&quot;&lt;&lt;(Msq_avg-(M_avg+M_avg*n) )/(T)&lt;&lt;</td><td></td></tr><tr><td colspan="2"></td><td>t&quot;&lt;&lt;(Msq_avg-(Mabs_avg*Mabs_avg*n))/(T)&lt;&lt;//susceptibility per spin (X′)</td><td>//&lt;E&gt;;&lt;E^2&gt; per spin</td></tr><tr><td colspan="2"></td><td></td><td>//heat capacity () per spin</td></tr><tr><td colspan="2">//flip spin at given position</td><td>&quot;t&quot;&lt;&lt;E_avg&lt;&lt;&quot;\t&quot;&lt;&lt;Esq_avg&lt;&lt;</td><td></td></tr><tr><td colspan="2"></td><td>t&quot;&lt;&lt;(Esq_avg-(E_avg*E_avgn))/(T*T) &lt;&lt;</td><td></td></tr><tr><td colspan="2">flip(lat_type pos)</td><td></td><td></td></tr><tr><td colspan="2"></td><td>&quot;t&quot;&lt;&lt;1-((Mq_avg)/(3*Msq_avg)) &lt;&lt;endl;</td><td></td></tr><tr><td colspan="2"></td><td></td><td>//cumulant (U_L)</td></tr><tr><td colspan="2"></td><td></td><td></td></tr><tr><td colspan="2"></td><td></td><td></td></tr><tr><td colspan="2"></td><td></td><td></td></tr><tr><td colspan="2">lat[pos.x] [pos.y]=-lat[pos.x] [pos.y];</td><td></td><td></td></tr><tr><td colspan="2"</table>