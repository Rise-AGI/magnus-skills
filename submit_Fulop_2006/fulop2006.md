# A Numerical Study of Landau Damping with PETSc-PIC

Daniel S. Finn $^ { - 1 }$ , Matthew G. Knepley $^ 2$ , Joseph V. Pusztay $^ 3$ , and Mark F. Adams

$^ { 1 , 2 , 3 }$ University at Buffalo, The State University of New York $^ 4$ Lawrence Berkeley National Laboratory

March 23, 2023

# Abstract

We present a study of the standard plasma physics test, Landau damping, using the Particle-In-Cell (PIC) algorithm. The Landau damping phenomenon consists of the damping of small oscillations in plasmas without collisions. In the PIC method, a hybrid discretization is constructed with a grid of finitely supported basis functions to represent the electric, magnetic and/or gravitational fields, and a distribution of delta functions to represent the particle field. Approximations to the dispersion relation are found to be inadequate in accurately calculating values for the electric field frequency and damping rate when parameters of the physical system, such as the plasma frequency or thermal velocity, are varied. We present a ful derivation and numerical solution for the dispersion relation, and verify the PETSC-PIC numerical solutions to the Vlasov-Poisson for a large range of wave numbers and charge densities.

Key words— Simulation, Particle-In-Cell, PETSc, Plasma, Landau damping

# 1 Introduction

In 1936, Lev Landau first formulated a simple kinetic model, now referred to as the Fokker-Plank equation in Landau form or simply just the Landau equation, for the description of charged particles in a plasma performing Coulomb collisions [24]. Ten years later, Landau furthered this discovery by predicting the damping of non-relativistic, collisionless plasma oscilations, or Langmuir waves, for the first time [25]. The basic concept proposed in that paper, that a conservative phenomenon exhibits irreversible behaviors, has since influenced hundreds of papers and become one of the foundational problems in plasma physics. Thus, the phenomenon is now referred to as Landau damping. In his seminal paper, Landau used the solution to the Cauchy problem for the linearized Vlasov-Poisson equation around a spatially homogeneous Maxwellian equilibrium. Landau solved the equation analytically using Fourier and Laplace transforms and concluded that the electric field damps exponentially and that the decay is a function of the wavenumber, $k$ , of the perturbation. In [5], Bohm and Gross provide a simple explanation for the damping in plasmas. In essence, plasmas exhibit a tendency to remain approximately field free. Therefore, if electric fields are introduced, either by external disturbance or by an incomplete space charge neutralization, the newly introduced fields will be forced out by a reaction from the free charges.

Through the years, numerous others have extensively examined Landau damping in literature [12, 19, 35]. In 2009, a rigorous solution to the nonlinear Vlasov-Poisson equation was given by Villani and Mouhot in [30]. In their paper, the damping phenomenon is reinterpreted in terms of transfer of regularity between kinetic and spatial variables, rather than exchanges of energy, with phase mixing being the driving mechanism.

Developed in parall to the theory behind Landau damping, numerical methods for approximating solutions to the kinetic plasma system were pioneered by Vlasov [36]. The Particle-In-Cell (PIC) method has been a popular choice for numerically simulating plasmas since its inception [17, 18], as it can considerably reduce the complexity of the system in comparison to a direct $N$ -body methods. The PIC method is a hybrid discretization algorithm comprised of two separate sets of bases for evaluation of different aspects of the problem. These bases are the particle basis, where the particle is represented by some (usually radially symmetric) shape function, and the mesh basis, where a mean field approach may be taken to computing different field quantities from external and self consistent forces. Typically, the continuum field solve is handled by employing the finite element method, although other formulations have used splines [10], finite difference methods [4], etc.

In this paper, we present a Particle-In-Cell (PIC) method for solving the Vlasov-Poisson system using the Portable Extensible Toolkit for Scientific Computing (PETSc) [2, 3]. PETSc-PIC uses symplectic integration schemes [1] for particle pushing while conducting feld solves with a finite element method [21, 26]. The gal of PETSc is to provide composable pieces from which optimal simulations can be constructed. PETSc user level APIs allow applications to delay implementation choices, such as solver details, until runtime using dynamic configuration [6]. PETSc-PIC solvers fully conserve the moments, mass, momentum and energy, at each time step while also preserving entropy monotonicity. Recent advances in the PETSc-PIC code [32] alsoinclude conservative projections between the finite element and partice basis, a keysteptowards ybrd FEM-particle algorithms.

# 2 Problem Formulation

Consider the Vlasov-Poisson system, a common variation of the more general Vlasov-Maxwell-Landau system of equations in the non-relativistic case where the magnetic and collisional ffects are neglected. It can be an effective model for strongly non-Maxwellian plasmas. The Vlasov equation,

$$
\frac { \partial f } { \partial t } + \mathbf { v } \cdot \frac { \partial f } { \partial \mathbf { x } } - \frac { q _ { e } } { m } \mathbf { E } \cdot \frac { \partial f } { \partial \mathbf { v } } = 0 ,
$$

describes the evolution of the phase space distribution, $f ( \mathbf { x } , \mathbf { v } , t )$ , defined over the domain $( \mathbf { x } , \mathbf { v } ) \in \mathbb { R } ^ { D } \times \mathbb { R } ^ { D }$ where $D$ is the spatial dimension. The electric field is obtained using Poisson's equation,

$$
\Delta \phi ( { \bf x } , t ) = - \frac { \rho } { \epsilon _ { 0 } } ,
$$

where $\phi$ is the electric potential, $\rho$ is the charge density and $\mathbf { E } = - \nabla \phi$ . The charge density contains a neutralizing background term, $\sigma$ , such that,

$$
\rho ( \mathbf { x } , \mathbf { v } , t ) = \sigma - q _ { e } \int _ { \mathbb { R } ^ { D } } f ( \mathbf { x } , \mathbf { v } , t ) d \mathbf { v } .
$$

This neutralizing background simulates the effect of ions on the electrons in the domain. The use of a stationary, uniform background charge is based on the assumption that the ions are much heavier than the electrons and thus feel little influence from them.

In order to study the linear Landau damping phenomenon, we consider the initial particle distribution,

$$
\displaystyle \boldsymbol { f } ( \boldsymbol { x } , \boldsymbol { v } , t = 0 ) = \frac { 1 } { \sqrt { 2 \pi v _ { t h } ^ { 2 } } } e ^ { - v ^ { 2 } / 2 v _ { t h } ^ { 2 } } ( 1 + \alpha \cos ( k x ) )
$$

where $v _ { t h } = \sqrt { K T _ { e } / m }$ , $\alpha = 0 . 0 1$ , $k = 0 . 5$ , $v _ { m a x } = 1 0$ and the boundaries are periodic. An important piece of PIC methods for the Vlasov-Poisson system is the reduction of noise. A primary source of noise in PIC methods can be traced to the initial discrete distribution of particles in phase space. We mimic a "quiet start" [7, 13, 14] continuum initialization in this work by placing particles at the center of the spatial and velocity cells and weighting them based on the initial distribution function $f ( x , v , t = 0 )$ . This method is also used in [31], where particles are further remapped back to the cell centers every few steps. The remapping step provides enough noise reduction to accurately observe nonlinear effects in damping, however, as we are, for now, concerned only with the linear case of Landau damping, we will ignore the remapping phase.

# 2.1 Linear Landau Damping

We seek to first derive a set of equations to understand the damping of plasma oscillations in our system and to calculate expected values for the damping rate and electric field oscilation frequency. These

expressions are found by first deriving the dispersion relation for a plasma. The derivation shown follows from [9]. Consider a uniform plasma with an initial distribution $f _ { 0 } ( v )$ with zero initial electric and magnetic fields, ${ \bf E } _ { 0 } = { \bf B } _ { 0 } = 0$ . To first order, the perturbation in $f ( x , v , t )$ is denoted by $f _ { 1 } ( x , v , t )$ such that,

$$
f ( x , v , t ) = f _ { 0 } ( v ) + f _ { 1 } ( x , v , t ) .
$$

Plugging (5) in to (1) gives,

$$
\frac { \partial f _ { 1 } } { \partial t } + \mathbf { v } \cdot \frac { \partial f _ { 1 } } { \partial \mathbf { x } } - \frac { q _ { e } } { m _ { e } } \mathbf { E } _ { 1 } \cdot \frac { \partial f _ { 0 } } { \partial \mathbf { v } } = 0 .
$$

Assuming that the ions are massive and fixed and that the waves are one-dimensional plane waves $f _ { 1 } \propto e ^ { i ( k x - \omega t ) }$ , (6) becomes,

$$
f _ { 1 } = \frac { i q _ { e } E _ { x } } { m _ { e } } \frac { \partial f _ { 0 } / \partial v _ { x } } { \omega - k v _ { x } } .
$$

Recall the Poisson equation (2), with the potential $\phi$ replaced by the divergence of the electric field,

$$
\begin{array} { l } { \displaystyle \nabla \cdot { \mathbf E } = \nabla \cdot { \mathbf E } _ { 1 } = - \frac { \rho } { \epsilon _ { 0 } } } \\ { = - \frac { 1 } { \epsilon _ { 0 } } \left( \sigma - q _ { e } \int \left( f _ { 0 } ( v ) + f _ { 1 } ( x , v , t ) \right) d v \right) . } \end{array}
$$

With zero initial electric field, the electric field vector is replaced by the electric perturbation, $\mathbf { E } _ { 1 }$ , which takes the form $\mathbf { E _ { 1 } } = E _ { x } e ^ { i ( k x - \omega t ) } \hat { \mathbf { x } }$ Furthermore, at equilibrium, the neutralizing background is equal to the total weight of the electron distribution, $\sigma = q _ { e } \int f _ { 0 } d v$ , leaving only the perturbation term $f _ { 1 }$ in the Poisson equation. Thus we are left with,

$$
i k \epsilon _ { 0 } E _ { x } = - q _ { e } \int f _ { 1 } d v .
$$

Substituting (7) into (9) and dividing by $i k \epsilon _ { 0 } E _ { x }$ , we have,

$$
1 = \frac { q _ { e } ^ { 2 } } { k m _ { e } \epsilon _ { 0 } } \int \frac { \partial f _ { 0 } / \partial v } { \omega - k v } d v .
$$

The integral in (10) is a three-dimensional integral, however for the Maxwellian or other factorable distribution, integration in the 2nd and 3rd dimension is simple.Evaluating the (10) integral in the 2nd and 3rdimsion, nsubstiutnthe pase $\omega _ { p } = \left( n _ { e } q _ { e } ^ { 2 } / m \epsilon _ { 0 } \right) ^ { 1 / 2 }$ , lp n,

$$
1 = \frac { \omega _ { p } ^ { 2 } } { k ^ { 2 } } \int _ { - \infty } ^ { \infty } \frac { \partial f _ { 0 } / \partial v _ { x } } { v _ { x } - ( \omega / k ) } d v _ { x } .
$$

Landau showed that this problem can be solved rigorously by means of the Laplace transform method. Importantly, it is necessary to go around the singularity in the integrand in (11) in the complex plane. The solution to (11) takes the form,

$$
\omega = \omega _ { r } + i \gamma ,
$$

where $\omega _ { r }$ represents the real oscillations of the plasma and $\gamma$ the imaginary, which Landau showed to be the par o the solution divng the damping f the cillations.FoowigLandau's method [9], aapproiatin for the oscillation and damping terms can be derived, given by,

$$
\begin{array} { l } { \displaystyle \omega _ { r } = 1 + \frac { 3 } { 2 } \hat { k } ^ { 2 } , } \\ { \gamma = - \sqrt { \frac { \pi } { 8 } } \frac { 1 } { \hat { k } ^ { 3 } } \exp \left[ - \frac { 1 } { 2 \hat { k } ^ { 2 } } \right] . } \end{array}
$$

A normalized form of the wavenumber $k$ has been introduced to simplify the equations going forward. The normalized wavenumber, $\hat { k }$ is given by,

$$
\hat { k } = \frac { k v _ { t h } } { \omega _ { p } }
$$

where $v _ { t h } = \sqrt { K T / m }$ is the thermal velocity. For all examples, we non-dimensionalize so that $v _ { t h } = 1$ . The real part of the solution to (11) was similarly derived by Vlasov in [36], however Vlasov did not account for the imaginary damping term.

These approximations are valid for the case where $\hat { k } \ll 1$ but their accuracy degrades considerably as $\hat { k }$ approaches 1 and higher. Even when $\hat { k } = 0 . 5$ , the calculated values for $\omega _ { r }$ and $\gamma$ differ from the numerical results by at least $5 \%$ . In [29], McKinstrie draws similar conclusions, electing to derive more accurate forms of (13) by expanding $\omega _ { r }$ in powers of $\hat { k }$ ,

$$
\begin{array} { l } { \displaystyle \omega _ { r } = 1 + \frac { 3 } { 2 } \hat { k } ^ { 2 } + \frac { 1 5 } { 8 } \hat { k } ^ { 4 } + \frac { 1 4 7 } { 1 6 } \hat { k } ^ { 6 } , } \\ { \displaystyle \gamma = - \sqrt { \frac { \pi } { 8 } } \left( \frac { 1 } { \hat { k } ^ { 3 } } - 6 \hat { k } \right) \exp \left[ - \frac { 1 } { 2 \hat { k } ^ { 2 } } - \frac { 3 } { 2 } - 3 \hat { k } ^ { 2 } - 1 2 \hat { k } ^ { 4 } \right] . } \end{array}
$$

These new expressions are more accurate for $\hat { k }$ up to 0.4 but still diverge from the correct values as $\hat { k }$ increases further. Shalaby et. al. provided further refinements to these equations in [33], using a numerical fitting formula, taking the form,

$$
\begin{array} { r l } & { \omega = 1 + \cfrac { 3 } { 2 } \hat { k } ^ { 2 } + \cfrac { 1 5 } { 8 } \hat { k } ^ { 4 } + \cfrac { 1 4 7 } { 1 6 } \hat { k } ^ { 6 } + 7 3 6 . 4 3 7 \hat { k } ^ { 8 } - 1 4 7 2 9 . 3 \hat { k } ^ { 1 0 } } \\ & { \qquad + 1 0 5 4 2 9 \hat { k } ^ { 1 2 } - 3 7 0 1 5 1 \hat { k } ^ { 1 4 } + 6 4 5 5 3 8 \hat { k } ^ { 1 6 } - 4 4 8 1 9 0 \hat { k } ^ { 1 8 } , } \\ & { \gamma = - \sqrt { \cfrac { \pi } { 8 } } \left( \cfrac { 1 } { \hat { k } ^ { 3 } } - 6 \hat { k } - 4 0 . 7 1 7 3 \hat { k } ^ { 3 } + 3 9 0 0 . 2 3 \hat { k } ^ { 5 } - 2 4 6 2 . 2 5 \hat { k } ^ { 7 } - 2 7 4 . 9 9 \hat { k } ^ { 9 } \right) } \\ & { \qquad \exp \left[ - \cfrac { 1 } { 2 \hat { k } ^ { 2 } } - \cfrac { 3 } { 2 } - 3 \hat { k } ^ { 2 } - 1 2 \hat { k } ^ { 4 } - 5 7 5 . 5 1 6 \hat { k } ^ { 6 } + 3 7 9 0 . 1 6 \hat { k } ^ { 8 } \right. } \\ & { \qquad \left. - 8 8 2 7 . 5 4 \hat { k } ^ { 1 0 } + 7 2 6 6 . 8 7 \hat { k } ^ { 1 2 } \right] . } \end{array}
$$

These equations give good estimates for $\omega _ { r }$ and $\gamma$ in the case where $\hat { k } = 0 . 5$ , which is of particular interest in this paper. In fact, the values obtained from (17) in the case where $\hat { k } = 0 . 5$ and all other parameters $( \omega _ { p } , v _ { t h } , q _ { e } , \mathrm { e t c . } )$ are assumed to be 1.0 match those commonly listed as "analytic solutions" [9, 31, 38]. That being said, the accuracy of the numerical fit still decreases considerably for $\hat { k } > 0 . 6$ .

An alternate, and as we will show, more accurate way to calculate $\omega _ { r }$ and $\gamma$ for given values of $\hat { k }$ is to find them by computing the zeros of (11). This was done by Canosa in [8] for values of $k$ ranging from 0.25 to 2.0 in increments of 0.05 (see Table 1 for a selection of values). A comparison of the approximations by Landau, McKinstrie and Shalaby to the zero finding results from Canosa is shown in Section 4.

<table><tr><td>k</td><td>Wr</td><td>γ</td></tr><tr><td>0.25</td><td>1.1056</td><td>-0.0021693</td></tr><tr><td>0.5</td><td>1.4156</td><td>-0.15336</td></tr><tr><td>0.75</td><td>1.7371</td><td>-0.46192</td></tr><tr><td>1.0</td><td>2.0459</td><td>-0.85134</td></tr><tr><td>1.5</td><td>2.6323</td><td>-1.7757</td></tr><tr><td>2.0</td><td>3.1891</td><td>-2.8272</td></tr></table>

# 3 PETSc-PIC

PETSc, the Portable, Extensible Toolkit for Scientific Computation, is a well-known library for numerical methods. It provides parallel data management, structured and unstructured meshes, linear and nonlinear algebraic solvers and preconditioners, optimization algorithms, time integrators and many more functions. The PETSc-PIC algorithm relies on two modules to handle the particle and mesh solves simultaneously. The first, DMPlex [22, 23, 26], is a PETSc module for generic unstructured mesh creation, manipulation, and I/O [16]. It decouples user applications from the implementation details of common mesh and discretization tasks. The other important module for this work, DMSwarm [28], provides a fully parallel solution for pure particle methods (e.g. DEM, SPH, EFG) and for particle-mesh methods (e.g. PIC, FLIP, MPM, GIMP).

We start with discussion of the particle methods in the PETSc-PIC algorithm. A method must first be chosen to represent the particle space, and for interpolation between the mesh and particle representations. Thereare numerous choices in shape functions for this purpose, however in ur case asimple elta functin representation of particles is chosen.Thus the approximation of the distribution function is defined in the particle space as,

$$
f _ { p } = \sum _ { p } \overrightarrow { \omega _ { p } } \delta \left( \mathbf { x } - \mathbf { x _ { p } } \right) ,
$$

where $\overrightarrow { \omega _ { p } }$ is the vector of weights, $\mathbf { x }$ are the configuration space variables and $\mathbf { x _ { p } }$ represents the particle position and velocity, respectively. The finite element representation, using a function space $\nu$ , is given by the weighted sum of basis functions,

$$
f _ { F E } = \sum _ { i } f _ { i } \psi _ { i } ( \mathbf { x } ) ,
$$

vhere $\psi _ { i } \in \mathcal V$ denotes the basis functions and $f _ { i }$ the associated finite element coefficient.

The Vlasov equation is a linear hyperbolic equation which may be written in a simpler form,

$$
\frac { \partial f } { \partial t } + \mathbf { z } \cdot \nabla _ { \mathbf { q } } f = 0 ,
$$

where $\mathbf { q } = ( \mathbf { x } , \mathbf { v } )$ is the phase space variable and $\mathbf { z } = ( \mathbf { v } , - q _ { e } \mathbf { E } / m )$ is the combined force. The force term $- q _ { e } \mathbf { E } / m$ is independent of velocity, and therefore (20) may be written in the conservative form,

$$
\frac { \partial f } { \partial t } + \nabla _ { \mathbf { q } } \cdot ( \mathbf { z } f ) = 0 .
$$

Given this new advective form of the Vlasov equation, we can rewrite the equation for the characteristics $\mathbf { Q } = \left( \mathbf { X } , \mathbf { V } \right)$ ,

$$
\frac { d \mathbf { Q } } { d t } = \mathbf { z } ,
$$

which reexpressed with the original phase-space variables gives,

$$
\begin{array} { l } { \displaystyle \frac { d { \bf X } } { d t } = { \bf V } , } \\ { \displaystyle \frac { d { \bf V } } { d t } = - \frac { q _ { e } } { m } { \bf E } . } \end{array}
$$

Since particles follow characteristics, the Vlasov equation in the particle basis becomes

$$
\begin{array} { l } { \displaystyle \frac { d \mathbf { x } _ { p } } { d t } = \mathbf { v } _ { p } , } \\ { \displaystyle \frac { d \mathbf { v } _ { p } } { d t } = - \frac { q _ { e } } { m } \mathbf { E } . } \end{array}
$$

The equations of motion are stepped forward in time using structure-preserving symplectic integrators which have been we studied [15] Theelectric eld is solve cncurently at each step using fnitemt solver, discussed in the next section.

# 3.1 PETSc-FEM

At each step in the simulation, the Poisson equation is solved using the finite element method. The rae  the potentil theelecel  thenterpoatco eacceat the partiec. The interpolated electric field is then applied to the particles in the form of the Coulomb force.

The PETSc-FEM method is abstractly formalized by the Ciarlet triple [11, 20], such that a finite elemen is a triple $( \mathcal { T } , \mathcal { V } , \mathcal { V } ^ { \prime } )$ , where,

the domain $\tau$ is a bounded, closed subset of $\mathbb { R } ^ { d }$ (for $d = 1 , 2 , 3 , \ldots$ ) with nonempty interior and piecewise smooth boundary;   
the space $\nu = \nu ( \Omega )$ is a finite-dimensional function space on $\tau$ of dimension $n$ ;   
•the set of degrees of freedom (nodes) $\mathcal { V } ^ { \prime } = \{ l _ { 1 } , l _ { 2 } , . . . , l _ { n } \}$ is a basis for the dual space, that is, the space of bounded linear functionals on $\nu$ .

The cell $\tau$ together with the local function space $\nu$ and the set of rules for describing the functions in $\nu$ is the finiteelement. The discretization in PETSc is handled by the PETScFE object, which contains a PetscSpace $\nu$ ), PetscDualSpace ( $\mathcal { V } ^ { \prime }$ ), and DMPlex ( $\tau$ ). PETScFE supports simplicial elements, tensor cells, and some special cells such as pyramids.

In general, the finite element solve for the Poisson equation can be accomplished using the standard $H ^ { 1 }$ function space. In the $H ^ { 1 }$ space, the weak form of the Poisson equation is,

$$
\int _ { \Omega } \nabla \psi _ { i } \cdot \nabla \phi = \int _ { \Omega } \psi _ { i } ,
$$

where $\psi \in V$ and $V$ is the set of basis functions on the cell. The elements are then constructed such that the basis functions are continuous across the cell boundaries.

# 3.2 Conservative Projections

To preserve the conservation laws in a PIC simulation, a method must be constructed to conservatively project between the particle and grid representations. Weak equality of the representations,

$$
\int _ { \Omega } \psi _ { i } f _ { F E } = \int _ { \Omega } \psi _ { i } f _ { P }
$$

is enforced on the representations to achieve this [27, 32]. Restricting this equivalence to the finitedimensional analogues gives the matrix-vector form,

$$
M f _ { F E } = M _ { p } f _ { p } ,
$$

where M is the finite element mass matrix,

$$
M = \int _ { \Omega } \psi _ { i } \psi _ { j } ,
$$

$M _ { p }$ is the particle mass matrix,

$$
M _ { p } = \int _ { \omega } \psi _ { i } \delta ( { \bf x } - { \bf x _ { p } } ) ,
$$

$f _ { F E }$ is a vector containing the finite element coefficients and $f _ { p }$ is the vector of particle weights. The entries of $M _ { p }$ contain evaluations of the finite element basis functions at particle locations with rows being determined by the basis function index, and columns being determined by the particle indices. Moving from the particle basis to the mesh, we must invert the finite element mass matrix, which is easily accomplished with CG/Jacobi [37]. In the other direction, we must invert a rectangular particle mass matrix, usually with LSQR [32].

# 4 Numerical Results

In this section, the results of this numerical study are presented. We consider the one-dimensional (1X1V) case of the Vlasov-Poisson system. According to (15), derived in Section 2.1, and the zero finding data from Canosa [8], the damping rate should be $\gamma = - 0 . 1 5 3$ and the frequency of oscillations should be $\omega _ { r } = 1 . 4 1 6$ . All runs were conducted on a single 2.4 $[ G H z ]$ 8-Core Intel Core i9 processor with 64 $\left[ G B \right]$ of memory. The example code and packages/options required to run it are provided in Appendix A.

To begin, we show results from the densest run of the PETSc-PIC simulation with 160 spatial cells and $8 , 0 0 0$ particles per cell and a PIC timestep of $d t = 0 . 3$ . Figure 2 shows the maximum value of the electric field, $E _ { m a x } = \operatorname* { m a x } _ { \Omega } | E |$ , over time. The values for $\gamma$ and $\omega _ { r }$ were measured by fitting the peaks of the given dataTheequencyatinsdbesheequenc thelecrelcletinnuci. Since each oscillation of $E _ { m a x }$ $E _ { m a x }$ oscillations for each plasma oscillation. Values achieved by the PETSc-PIC algorithm, $\gamma = - 0 . 1 5 3 1$ and $\omega _ { r } = 1 . 4 1 2 4$ agree within $1 \%$ of the analytic values from Canosa and Shalaby et. al., which are assumed to be the most accurate for the case $k = 0 . 5$ .

![](images/2d5f6a6c2715d85233e26716f04b3bf23103c101be188e526dabf461b69768aa.jpg)  
Figure 2: The maximum value of the electric field as a function of time for the one-dimensional linea Landau damping problem.

The total error in the moments, shown in Figure 3, was shown to be stable over the entire runtime. At early times, the error in momentum and energy fluctuate but each converges by $t = 2 0$ . This convergence comes from the use of basic symplectic integrator in PETSc-PIC which guarantees the error does not grow over time. We also note that the error in the particle solve and the finite element solve is exactly equal, apartfromancease levl  noisthe partice olveThisconirms the efectivenes thecnseaive projector used in PETSc-PIC. More detailed tests of the conservative projector can be found in [32].

Convergence studies were conducted in which either the mesh or the number of particles per cell were increased while the other was held constant. In the case of mesh convergence, the number of particles per cell was held at $N _ { v } = 8 , 0 0 0$ while in the particle number convergence tests, the number of mesh cells was $N _ { x } = 1 0 0$ . We naively expect Monte Carlo convergence in particle number, $\mathcal { O } \left( 1 / \sqrt { N } \right)$ , and we indeed achieve this for $\omega _ { r }$ in the upper left of Figure 4. Since we have such a regular initial particle distribution, we might hope to see Quasi-Monte Carlo convergence, $\mathcal { O } \left( 1 / N \right)$ , and we do see this superconvergence in $\gamma$ in the upper right of Figure 4. We expect $\mathcal { O } \left( h ^ { 2 } \right)$ convergence in the mesh resolution $h$ since this controls the erro in theelectric fel, and we se this in gama in the lower right  Figure However, this conveence should quickly saturate as particle error begins to dominate, which we see in $\omega _ { r }$ in the lower left of Figure 4.

![](images/00309a25761c5fb3fc5adae9cb8350863bab40fbb6b4158e56d1cd1760ce55e2.jpg)  
Figure 3: The total mass, momentum and energy error for the particle and finite element solve. The moment errors all converge to zero given a long enough time.

# 4.1 Variations in Wavenumber and Charge Density

We have thus far shown that the PETSc-PIC algorithm is an accurate and structure-preserving method for modeling plasma systems. We next present results from tests in which the wavenumber $k$ , and consequently the domain size, and the charge density were varied.Varying either of these values impacts the value of the non-dimensional wavenumber $\hat { k }$ . In the case of the wavenumber $k$ , the calculated values for $\omega _ { r }$ and $\gamma$ were compared to the values obtained with the approximation equations (13) and (15), the numerical fit (17) and the zero finding results from Table 1. The results from PETSc-PIC, shown in Figure 5, clearly show a strong deviation of the approximation equations and the numerical fit for $\hat { k } > 0 . 5$ while closely matching the zero finding data. This demonstrates that these approximations quickly break down outside of the small parameter range typically chosen in numerical studies of Landau damping. When considering real plasma systems in which values for $k$ , $\omega _ { p }$ , etc. are more dynamic, it is far more effective to use zero finding methods to calculate expected values for $\omega _ { r }$ and $\gamma$ .

It may be naively assumed that data from numerical tests with varying charge densities will match the approximation equations (13), (15) and (17) or even the zero finding data from Canosa, however, these analytic results are based on an assumption of unchanging charge density. More specifically, these results are based on charge densities such that the plasma frequency, $\omega _ { p }$ , is always unity. Therefore, to accurately compare analytic results to our data we must resolve the dispersion relation for varying charge densities. A zero finding algorithm, using Newton's method [34], was employed to calculate new analytic values for $\omega _ { r }$ and $\gamma$ with charge densities ranging from 0.1 to 2.0. The zero finding algorithm calculates multiple values for $\omega _ { r }$ and $\gamma$ however we select the solution containing the largest $\gamma$ which corresponds to the smallest $\omega _ { r }$ . Other solutions found by the algorithm represent less dominant modes which can be ignored for the purposes of this study. Figure 6 contains the results from the new zero finding algorithm along with data from numerical tests which agree perfectly. We observe that when the charge density is increased, the frequency of oscillations also increases. This matches the expected physical behavior of an electrically charged plasma.

![](images/68f12460405cfd464770aaef9d934c812409a08673902bbbbcf20549f86454d1.jpg)  
Figure 4: (Top) Particle per cell number convergence plots for $N _ { x } = 1 0 0$ and (bottom) mesh convergence plots for $N _ { v } = 8 , 0 0 0$ .

We have extended our zero finding algorithm to the case where the charge density approaches zero $\hat { k } \to \infty$ ) to make note of an interesting phenomenon. At a charge density of zero, the dispersion relation has no solution. We capture this in Figure 6, where we observe that both $\omega _ { r }$ and $\gamma$ are asymptotic at $\int f _ { 0 } / V = 0$ . This can similarly be observed in the numerical results from our PETSc-PIC algorithm. As the charge density is decreased, the rate at which the electric feld oscillations becomes too large to resolve numerically. In the case of $\int f _ { 0 } / V = 0 . 2 5$ , shown in Figure 7, we can only observe two full oscillations of the elecric feld before the simulation becomes too noisy. Theoretically, the charge density could be decreased asymptotically in our simulations to observe the damping rate and frequency trends but in practice there is too much noise to resolve any real processes in the plasma.

![](images/3e64b34c29c4aa11fcfec862f2030cc29716f0f7e84b01a1fea87d5f00800ba9.jpg)  
Figure 5: A comparison of various approximations for $\omega _ { r }$ and $\gamma$ to root finding results and numerical results from PETSc-PIC. Plots on the right are zoomed in on the region $0 . 0 \leq \hat { k } \leq 0 . 7 5$ to show the accuracy of each approximation before they diverge from the data.

# 5 Conclusions and Future Work

We have presented PETSc-PIC, a structure-preserving Particle-in-Cell algorithm for solving the electrostatic Vlasov-Poisson systems. The accuracy of our algorithm has been demonstrated by comparing the frequency of electric field oscillations and the damping rate of the oscillations to analytic values. We have also shown that the approximations or the frequency and damping rate break down outside of narrow ranges for the wavenumber and charge density. These approximations are often cited in numerical Landau damping studies without further context or reference to the equations used to compute the parameters, which can lead to complications in reproducing results. We have sought to provide a complete picture of Landau damping and the numerical methods we have used to simulate this phenomenon.

Future work with the PETSc-PIC algorithm will fall into two primary categories: improvements to the algorithm and the extension of the Landau damping test to multi-dimensional and nonlinear cases. Improvements to the algorithm willfocus on reformulation using a mixed form of the Poisson equation and $H ( d i v )$ finite elements. We expect that $C ^ { 0 }$ electric fields will decrease the noise in our particle representation over time. While we have not observed any major negative impacts from using $H ^ { 1 }$ finite elements in the test problem chosen for this paper, Landau damping, using a mixed form makes a notable difference in the case of Two-Stream Instability. PETSc currently includes support for the $H ( d i v )$ conforming finite elements BrezziDouglas-Marini (BDM) and Raviart-Thomas (RT) on simplicial grids, however RT elements are currently the only element type supported on tensor cells. We will also replicate our tests in parall, allowing us to increase the number of particles per cll by several orders of magnitude, reducing thelargest source o error in the code.

![](images/0012327789b3b7dc6d4714b8f433623b31b0c7ac30cff02f8b89bd3d4e9f435d.jpg)  
Figure 6: Numerical results for varying charge densities compared to zero finding data. The charge density is represented on the x-axis as the integral of the initial distribution over the domain volume.

![](images/fcb9f10c2742574d5be4e12c679fd3f44edf0b05922ed86605b3ec345ef06d2c.jpg)  
Figure 7: A comparison of electric field oscillations given different charge densities.

Nonlinear Landau damping is a more complex in that non-damping phenomenon, such as plasma echo, are present. Vitally, the linearization of the Vlasov equation, used as the fundamental approximation in the study of linear Landau damping, does not guarantee that the asymptotic behavior of the linear Vlasov equation is an approximation of the asymptotic behavior of the nonlinear Vlasov equation [30]. There are reasons to doubt that the study of the linearized equations gives any hint on the long-time behavior of the nonlinear equations.Therefore fan algorithm is desired that can accurately capture the long-time behavior of a plasma, the nonlinear case of Landau damping must also be considered.

# A Appendix A

The data presented in this paper can be recreated with PETSc using the DMSwarm example ex9 (\$PETSC DIR/src/dm/impls/swarm/tests/ex9.c). Exact runtimes may vary depending on the architecture and compiler. The DMSwarm example can be run using the following options:

./ex9 -dm_plex_dim 2 -dm_plex_simplex 0 -dm_plex_box_bd periodic,none -dm_plex_box_faces 10,1 -dm_view -dm_plex_box_lower 0,-0.5 -dm_plex_box_upper 12.5664,0.5   
-dm_swarm_num_species 1 -dm_swarm_num_particles 50   
-vdm_plex_dim 1 -vdm_plex_simplex 0 -vdm_plex_box_faces 7500 -vdm_plex_box_lower -10 -vdm_plex_box_upper 10   
-petscspace_degree 1 -em_type primal -em_pc_type svd -em_snes_atol 1.e-12   
-ts_type basicsymplectic -ts_basicsymplectic_type 1 -ts_max_time 500 -ts_max_steps 1000 -ts_dt 0.03   
-fake_1D -cosine_coefficients 0.01,0.5 -charges -1.0,1.0 -perturbed_weights -periodic

This example uses a 100 square cell mesh on the domain $( x , v ) \in [ 0 , 4 \pi ] \times [ - 1 0 , 1 0 ]$ , with 8000 particles per cell. A first-order basic symplectic integrator is chosen as the time integration method and $H ^ { 1 }$ finite elements are chosen for the field solves.

# References

[1] Shrirang Abhyankar et al. PETSc/TS: A Modern Scalable DAE/ODE Solver Library. Preprint ANL/MCS-P5061-0114. ANL, Jan. 2014.   
[2] Satish Balay et. al. PETSc/TAO Users Manual. English. Version Revision 3.18. Argonne National Laboratory. 2022. 310 pp.   
[3] Satish Balay et al. PETSc Web page. https://petsc.org/. 2022. URL: https://petsc.org/.   
[4] Jeffrey W. Banks et al. "High-Order Accurate Conservative Finite Difference Methods for Vlasov Equations in 2D+2V". In: SIAM Journal on Scientific Computing 41.5 (2019), B953-B982. DOI: 10.1137/19M1238551. eprint: https://doi.org/10.1137/19M1238551. URL: https://doi.org/10. 1137/19M1238551.   
[5] D. Bohm and E. P. Gross. "Theory of Plasma Oscillations. A. Origin of Medium-Like Behavior". In: Phys. Rev. 75 (12 1949), pp. 18511864. DOI: 10.1103/PhysRev.75.1851. URL: https://1ink.aps. org/doi/10.1103/PhysRev.75.1851.   
[6] Jed Brown, Matthew G. Knepley, and Barry Smith. "Run-time extensibility and librarization of simulation software". In: IEEE Computing in Science and Engineering 17.1 (Jan. 2015), pp. 3845. DOI: 10.1109/MCSE.2014.95.   
[7] J.A. Byers. "Noise Suppression Techniques in Macroparticle Models of Collisionless Plasmas". In: Proceedings of the Fourth Conference of Numerical Simulation of Plasmas. Ed. by NTIS. NRL, Washington, D.C., 1970.   
[8] José Canosa. "Numerical solution of Landau's dispersion equation". In: Journal of Computational Physics 13.1 (1973), pp. 158160. IsSN: 0021-9991. D01I: https://doi.org/10.1016/0021-9991(73) 90131-9. URL: https://www.sciencedirect.com/science/article/pii/0021999173901319. [9] F.F. Chen. Introduction to Plasma Physics and Controlled Fusion. 2nd ed. Vol. 1. Plenum Press, New York, 1984.   
[10] C.Z Cheng and Georg Knor. "The Integration of the Vlasov Equation in Configuration Space". In: Journal of Computational Physics 22.3 (1976), pp. 330351. IssN: 0021-9991. DO1: https://doi.org/ 10.1016/0021-9991(76)90053-X. URL: https://www.sciencedirect.com/science/article/pii/ 002199917690053X.   
[11] Philippe G. Ciarlet. Numerical Analysis of the Finite Element jethod. Les Presses de L'Université de Montréal, 1976.   
[12] John Dawson. "On Landau Damping". In: The Physics of Fluids 4.7 (1961), pp. 869874. DO1: 10. 1063/1.1706419. eprint: https://aip.scitation.org/doi/pdf/10 .1063/1.1706419. URL: https://aip.scitation.org/doi/abs/10.1063/1.1706419.   
[13] J. Denavit and J.M. Walsh. Comments on Plasma Physics and Controlled Fusion. eng. New York, 1988.   
[14] J.P. Friedberg, R.L. Morse, and C.W. Nielson. In: Proceedings of the Third Conference on Numerical Simulation of Plasmas. Standford University, CA, 1969.   
[15] Ernst Hairer, Gerhard Wanner, and Christian Lubich. "Symplectic Integration of Hamiltonian Systems". In: Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations. Berlin, Heidelberg: Springer Berlin Heidelberg, 2006, pp. 179236. ISBN: 978-3-540-30666-5. DOI: 10.1007/3-540-30666-8_6. URL: https://doi.org/10.1007/3-540-30666-8_6.   
[16] Vaclav Hapla et al. "Fully Parallel Mesh I/O Using PETSc DMPlex with an Application to Waveform Modeling". In: SIAM Journal on Scientific Computing 43.2 (2021), pp. C127C153. ISSN: 1095-7197. DOI: 10.1137/20m1332748. URL: http://dx.doi.org/10.1137/20M1332748.   
[17] F.H. Harlow, M. Evans, and R.D. Richtmyer. A Machine Calculation Method for Hydrodynamic Problems. LAMS (Los Alamos Scientific Laboratory). Los Alamos Scientific Laboratory of the University of California, 1955. URL: https://books.google.com/books?id=rvM5zQEACAAJ.   
[18] Francis H Harlow. The Particle-In-Cell Method for Numerical Solution of Problems in Fluid Dynamics. LAMS (Los Alamos Scientific Laboratory). Los Alamos Scientific Laboratory of the University of California, 1962. DOI: 10.2172/4769185. URL: https://www.osti.gov/biblio/4769185.   
[19] J D Jackson. "Longitudinal plasma oscillations". In: Journal of Nuclear Energy. Part C, Plasma Physics, Accelerators, Thermonuclear Research 1.4 (1960), p. 171. D01: 10.1088/0368-3281/1/4/301. URL: https://dx.doi.org/10.1088/0368-3281/1/4/301.   
[20] Robert C. Kirby. "Algorithm 839: FIAT, a new paradigm for computing finite element basis functions". In: ACM Transactions on Mathematical Software 30.4 (2004), pp. 502516. ISSN: 0098-3500. DOI: 10.1145/1039813.1039820.   
[21] M. G. Knepley et al. "Achieving High Performance with Unified Residual Evaluation". In: ArXiv e-prints (Sept. 2013). arXiv: 1309.1204 [cs.MS].   
[22] Matthew G. Knepley and Dmitry A. Karpeev. "Mesh Algorithms for PDE with Sieve I: Mesh Distribution". In: Scientific Programming 17.3 (2009). http://arxiv.org/abs/0908.4427, pp. 215 230. DOI: 10.3233/SPR-2009-0249. URL: http://arxiv.org/abs/0908.4427.   
[23] Matthew G. Knepley, Michael Lange, and Gerard J. Gorman. Unstructured Overlapping Mesh Distribution in Parallel. 2017. eprint: 1506.06194.   
[24] L.D. Landau. "Die kinetische Gleichung für den Fall Coulombscher Wechselwirkung". In: Phys. Zs. Sowjetunion 10 (1936), pp. 154164.   
[25] L.D. Landau. "On the Vibrations of the Electronic Plasma". In: Journal of Physics USSR 10 (1946), pp. 2534.   
[26] Michael Lange et al. "Efficient mesh management in Firedrake using PETSc-DMPlex". In: SIAM Journal on Scientific Computing 38.5 (2016), S143S155. DO1: 10.1137/15M1026092. eprint: http: //arxiv.org/abs/1506.07749.   
[27] James R Maddison and Patrick E Farrell. "Directional integration on unstructured meshes via supermesh construction". In: Journal of Computational Physics 231.12 (2012), pp. 44224432.   
[28] Dave A May and Matthew G Knepley. "DMSwarm: Particles in PETSc". In: EGU General Assembly Conference Abstracts. Vol. 19. 2017, p. 10133.   
[29] C. J. McKinstrie, R. E. Giacone, and E. A. Startsev. "Accurate formulas for the Landau damping rates of electrostatic waves". In: Physics of Plasmas 6.2 (199), pp. 463466. D01: 10.1063/1.873212. eprint: https://doi.org/10.1063/1.873212. URL: https://doi.org/10.1063/1.873212.   
[30] Clement Mouhot and Cedric Villani. "On Landau Damping". In: Acta Mathematica 207.1 (2011), pp. 29-201. DOI: 10.1007/s11511-011-0068-9. URL: https://doi.org/10.1007/s11511-011- 0068-9.   
[31] Andrew Myers, Phillip Colella, and Brian Van Straalen. "A 4th-Order Particle-in-Cell Method with Phase-Space Remapping for the Vlasov-Poisson Equation". In: (2016). DOI: 10 .48550/ARXIV .1602. 00747. URL: https://arxiv.org/abs/1602.00747.   
[32] Joseph V. Pusztay, Matthew G. Knepley, and Mark F. Adams. "Conservative Projection Between Finite Element and Particle Bases". In: SIAM Journal on Scientific Computing 44.4 (2022), pp. C310 C319. DOI: 10.1137/21M1454079. eprint: https://doi .org/10.1137/21M1454079. URL: https: //doi.org/10.1137/21M1454079.   
[33] Mohamad Shalaby et al. "SHARP: A Spatially Higher-order, Relativistic Particle-in-cell Code". In: The Astrophysical Journal 841.1 (2017), p. 52. DO1: 10 .3847/1538-4357/aa6d13. URL: https : //dx.doi.org/10.3847/1538-4357/aa6d13.   
[34] E. Sonnendrücker. Numerical Methods for the VlasovMaxwell equations (Book in preperation).   
[35] N.G. Van Kampen. "On the Theory of Stationary Waves in Plasmas". In: Physica 21.6 (1955), pp. 949 963. ISSN: 0031-8914. DOI: https : //doi . org/10 . 1016/S0031-8914(55) 93068-8. URL: https : //www.sciencedirect.com/science/article/pii/S0031891455930688.   
[36] A. Vlasov. "Über die Schwingungseigenschaften des Elektronengases". Russian. In: Zh. Éksper. Teor. Fiz. 8 (1938), pp. 291318. ISSN: 0044-4510.   
[37] Andrew J Wathen. "Realistic eigenvalue bounds for the Galerkin mass matrix". In: IMA Journal of Numerical Analysis 7.4 (1987), pp. 449457.   
[38] Tie Zhou, Yan Guo, and Chi Wang Shu. "Numerical study on Landau damping". In: Physica D: Nonlinear Phenomena 157.4 (2001), pp. 322333. ISSN: 01672789. DO1: 10 . 1016/S0167-2789(01) 00289-5.