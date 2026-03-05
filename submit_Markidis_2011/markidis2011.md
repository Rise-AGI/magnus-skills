# The Energy Conserving Particle-in-Cell Method

Stefano Markidis and Giovanni Lapenta Centre for Plasma Astrophysics, Katholieke Universiteit Leuven, Celestijnenlaan 200B, B-3001 Leuven, Belgium

# Abstract

A new Particle-in-Cell (PIC) method, that conserves energy exactly, is presented. The particle equations of motion and the Maxwell's equations are differenced implicitly in time by the midpoint rule and solved concurrently by a Jacobian-free Newton Krylov (JFNK) solver. Several tests show that the finite grid instability is eliminated in energy conserving PIC simulations, and the method correctly describes the two-stream and Weibel instabilities, conserving exactly the total energy. The computational time of the energy conserving PIC method increases linearly with the number of particles, and it is rather insensitive to the number of grid points and time step. The kinetic enslavement technique can be effectively used to reduce the problem matrix size and the number of JFNK solver iterations.

Keywords: Energy Conserving Particle-in-Cell, Kinetic Plasma Simulations

# 1. Introduction

The Particle-in-Cell (PIC) method is one of the most used numerical methods for the solution of the collision-less kinetic equation of plasmas. The majority of PIC schemes has the property of conserving exactly the system total momentum [1, 2], while it does not conserve the system total energy. In fact typically PIC methods, that use explicit differentiation in time (explicit PIC schemes), tend to increase the total energy of the system by numerical heating [1, 2], while PIC methods, that use implicit differentiation in time (implicit PIC schemes), tend to decrease the total energy of the system by numerical cooling [3]. In the study of plasma physics instabilities, where one kind of energy is converted to another one, it is important to ensure that energy is not created spuriously by the numerical scheme in use. In fact, the numerical heating introduces spurious energy that can feed erroneously plasma instabilities, leading to unphysical results.

A new class of PIC methods were developed at beginning of the Seventies to address the problem of the total energy conservation. An "energy conserving" PIC method was first proposed by Lewis in 1970 [4]. However, the Lewis' "energy conserving" scheme conserves energy only in the limit of zero time steps, while accuracy errors preclude the possibility of having exact energy conservation with finite time steps [1, 2]. The scheme does not lead to the exact energy conservation, but instead to an improved energy conservation at a given time step if compared to momentum conserving PIC methods. The "energy conserving" method was derived from the Lagrangian formalism, and proved to be more robust against aliasing instabilities, such as the finite grid instability. An in-depth analysis of the "energy conserving" PIC method is presented by Langdon in Ref. [5] and in textbooks [1, 2].

Differently from previous " energy conserving" PIC methods, the proposed PIC scheme has the property of conserving exactly the total energy not only in the limit of zero time step, but also with finite time step. The new PIC method conserves the system total energy with a precision that depends only on the error tolerance value of the iterative solver in use. The scheme is not derived from the Lagrangian formalism but instead from the Newton approach. Energy conservation is achieved by using the midpoint integration rule for particle and field equations, and a convenient discretization of the discrete spatial operators.

The energy conserving PIC method requires the concurrent solution of the equation of motion for each particle and of the field equations for the electric and magnetic fields at each grid point, posing two challenges. First, the direct solution of the implicit numerical equations of PIC method implies the computation of a non-linear system. Non-linearity arises from the coupling between particles and fields variables through the interpolation functions of the PIC method. In the proposed PIC scheme, the non-linear equations are solved by a Newton Krylov solver [6, 7]. Despite the belief that the iterative solution of such equations could hardly converge [8], it has been proven that such PIC methods, that are called fully implicit, are convergent [9, 10]. The second challenge is that the energy conserving PIC scheme requires the solution of a very large matrix whose rank is of the order of the number of particles (the number of particles is considerably higher than the number of grid points in typical PIC simulations). For this reason, implementations of fully implicit PIC methods are based on the matrix-free Jacobian-free solvers to avoid the storage of the matrix, and Jacobian coefficients [9]. Previous implementations of fully implicit PIC methods, as those presented in Refs. [9, 10], were limited to electrostatic simulation and more importantly their formulation does not imply the total energy conservation. The new energy conserving PIC method is still based on the solution of coupled nonlinear equations by a Jacobian-free Newton Krylov (JFNK) solver, but on the contrary it is formulated for the electromagnetic case and conserves exactly the total energy.

This paper presents the algorithm, the properties, the implementation, the simulation and performance results of the energy conserving PIC method. It is organized as follows. Section 2 introduces the governing equations, shows their discretization in time and space and explains in detail the numerical algorithm. Sections 3,4,5 analyze the properties of the proposed method: the energy and momentum conservation, and the numerical stability. The implementation of an energy conserving PIC code is discussed in Section 6. Section 7 presents first the results of a Maxwellian plasma simulation, that is robust against the finite grid instability, and then of the two-stream and Weibel instabilities. The computational performance of the proposed method, and a technique to reduce the problem matrix size are shown in Section 8. Finally, Section 9 concludes the paper summarizing the algorithm and its properties. In Appendices A and B, a skeleton version of the energy conserving PIC code in Matlab/Octave programming language is provided.

# 2. Algorithm

In the PIC method $N _ { s }$ computational particles of the different $n _ { s }$ species with label $s$ mimic the real behavior of electrons and ions [11]. Each computational particle is characterized by a position $\mathbf { x } _ { p }$ and a velocity $\mathbf { v } _ { p }$ , whose evolution is described by the equation of motion (here and thereafter in CGS units):

$$
\left\{ \begin{array} { l } { d \mathbf { x } _ { p } / d t = \mathbf { v } _ { p } } \\ { d \mathbf { v } _ { p } / d t = q _ { s } / m _ { s } \left( \mathbf { E } _ { p } + \mathbf { v } _ { p } / c \times \mathbf { B } _ { p } \right) , } \end{array} \right.
$$

where $q _ { s } / m _ { s }$ are the charge to mass ratio of the species $s$ . $\mathbf { E } _ { p }$ , and $\mathbf { B } _ { p }$ are the electric and magnetic fields acting on the particle $p$ and they are calculated by interpolation from $\mathbf { E } _ { g }$ and $\mathbf { B } _ { g }$ , the values of the electric and magnetic field on the $N _ { g }$ grid points, through the use of the interpolation function

$W ( \mathbf { x } _ { g } - \mathbf { x } _ { p } )$

$$
\mathbf { E } _ { p } = \sum _ { g } ^ { N _ { g } } \mathbf { E } _ { g } W ( \mathbf { x } _ { g } - \mathbf { x } _ { p } ) \qquad \mathbf { B } _ { p } = \sum _ { g } ^ { N _ { g } } \mathbf { B } _ { g } W ( \mathbf { x } _ { g } - \mathbf { x } _ { p } ) .
$$

The Cloud-in-Cell interpolation functions [1, 2] are used:

$$
W ( \mathbf { x } _ { g } - \mathbf { x } _ { p } ) = \left\{ \begin{array} { l l } { 1 - | \mathbf { x } _ { g } - \mathbf { x } _ { p } | / \Delta x \quad \mathrm { i f } \quad | \mathbf { x } _ { g } - \mathbf { x } _ { p } | < \Delta x } \\ { 0 \quad \mathrm { o t h e r w i s e } . } \end{array} \right.
$$

Equations $1$ are differenced in time using the implicit midpoint integration rule [12]:

$$
\left\{ \begin{array} { l l } { \mathbf { v } _ { p } ^ { n + 1 } = \mathbf { v } _ { p } ^ { n } + q _ { s } / m _ { s } \Delta t ( \bar { \mathbf { E } } _ { p } + \bar { \mathbf { v } } _ { p } / c \times \bar { \mathbf { B } } _ { p } ) } \\ { \mathbf { x } _ { p } ^ { n + 1 } = \mathbf { x } _ { p } ^ { n } + \bar { \mathbf { v } } _ { p } \Delta t } \end{array} \right. ,
$$

where $n$ is the time level and the bar variables are the average in time of the quantities, and they are defined as:

$$
\begin{array} { c c } { { \displaystyle \bar { \mathbf { x } } _ { p } = ( \mathbf { x } _ { p } ^ { n + 1 } + \mathbf { x } _ { p } ^ { n } ) / 2 } } & { { \displaystyle \bar { \mathbf { v } } _ { p } = ( \mathbf { v } _ { p } ^ { n + 1 } + \mathbf { v } _ { p } ^ { n } ) / 2 } } \\ { { \displaystyle \bar { \mathbf { E } } _ { p } = \sum _ { g } ^ { N _ { g } } \bar { \mathbf { E } } _ { g } W ( \mathbf { x } _ { g } - \bar { \mathbf { x } } _ { p } ) } } & { { \displaystyle \bar { \mathbf { B } } _ { p } = \sum _ { g } ^ { N _ { g } } \bar { \mathbf { B } } _ { g } W ( \mathbf { x } _ { g } - \bar { \mathbf { x } } _ { p } ) } } \\ { { \displaystyle \bar { \mathbf { E } } _ { g } = ( \mathbf { E } _ { g } ^ { n + 1 } + \mathbf { E } _ { g } ^ { n } ) / 2 } } & { { \displaystyle \bar { \mathbf { B } } _ { g } = ( \mathbf { B } _ { g } ^ { n + 1 } + \mathbf { B } _ { g } ^ { n } ) / 2 . } } \end{array}
$$

It is possible to rewrite Equations 4 in terms of $\overline { { \mathbf { v } } } _ { p }$ after a series of algebraic manipulations [3, 9]:

$$
\begin{array} { r c l } { { \displaystyle \tilde { \mathbf { v } } _ { p } } } & { { = } } & { { \mathbf { v } _ { p } ^ { n } + \frac { q _ { s } \Delta t } { 2 m _ { s } } \bar { \mathbf { E } } _ { p } } } \\ { { } } & { { } } & { { } } \\ { { \displaystyle \bar { \mathbf { v } } _ { p } } } & { { = } } & { { \frac { \tilde { \mathbf { v } } _ { p } + \frac { q _ { s } \Delta t } { 2 m _ { s } c } ( \tilde { \mathbf { v } } _ { p } \times \bar { \mathbf { B } } _ { p } + \frac { q _ { s } \Delta t } { 2 m _ { s } c } ( \tilde { \mathbf { v } } _ { p } \cdot \bar { \mathbf { B } } _ { p } ) \bar { \mathbf { B } } _ { p } ) } { ( 1 + \frac { q _ { s } ^ { 2 } \Delta t ^ { 2 } } { 4 m _ { s } ^ { 2 } c ^ { 2 } } \bar { B } _ { p } ^ { 2 } ) } , } } \end{array}
$$

and the equation of motion becomes:

$$
\left\{ \begin{array} { l l } { \mathbf { v } _ { p } ^ { n + 1 } = 2 \bar { \mathbf { v } } _ { p } - \mathbf { v } _ { p } ^ { n } } \\ { \mathbf { x } _ { p } ^ { n + 1 } = \mathbf { x } _ { p } ^ { n } + \bar { \mathbf { v } } _ { p } \Delta t . } \end{array} \right.
$$

The evolution of the electric and magnetic fields is determined by solving the Maxwell's equations:

$$
\left\{ \begin{array} { l l } { \nabla \cdot \mathbf { E } = 4 \pi \rho } \\ { \nabla \cdot \mathbf { B } = 0 } \\ { 1 / c \partial \mathbf { E } / \partial t = \nabla \times \mathbf { B } - 4 \pi / c \mathbf { J } } \\ { 1 / c \partial \mathbf { B } / \partial t = - \nabla \times \mathbf { E } , } \end{array} \right.
$$

The implicit midpoint scheme is used to discretize the Maxwell's equations. The Faraday's and Ampere's laws are differenced in time as follows:

$$
\left\{ \begin{array} { l l } { \mathbf { E } _ { g } ^ { n + 1 } - \mathbf { E } _ { g } ^ { n } = c \nabla \times \bar { \mathbf { B } } _ { g } \Delta t - 4 \pi \bar { \mathbf { J } } _ { \mathbf { g } } \Delta t = c / 2 \nabla \times ( \mathbf { B } _ { g } ^ { n + 1 } + \mathbf { B } _ { g } ^ { n } ) \Delta t - 4 \pi \bar { \mathbf { J } } _ { \mathbf { g } } } \\ { \mathbf { B } _ { g } ^ { n + 1 } - \mathbf { B } _ { g } ^ { n } = - c \nabla \times \bar { \mathbf { E } } _ { g } \Delta t = - c / 2 \nabla \times ( \mathbf { E } _ { g } ^ { n + 1 } + \mathbf { E } _ { g } ^ { n } ) \Delta t . } \end{array} \right.
$$

The average current density $\mathbf { \bar { J } _ { g } }$ is calculated from the particle average positions and velocities by interpolation:

$$
\bar { \mathbf { J } } _ { g } = \sum _ { s } ^ { n _ { s } } \sum _ { p } ^ { N _ { s } } q _ { s } \bar { \mathbf { v } } _ { p } W ( \mathbf { x } _ { g } - \bar { \mathbf { x } } _ { p } ) / V _ { g } ,
$$

where $V _ { g }$ is the volume of cell $g$ .

The discretization of spatial operators must be chosen carefully to ensure that the vector identities, that are valid in the continuous space, hold on the discrete grid also. To achieve this, the Yee's lattice [13, 14] discretization of the spatial operators in Equations 12 is used. Taking a uniform rectangular grid for simplicity, the different components of the electromagnetic field and of the current densities are calculated on the cell center (half integer index) and on the cell vertices (integer index) according to the Yee's lattice configuration:

$$
\begin{array} { r l } & { \mathbf { E } _ { g } = ( E _ { i , j + 1 / 2 , k + 1 / 2 } ^ { x } , E _ { i + 1 / 2 , j , k } ^ { y } , E _ { i + 1 / 2 , j , k } ^ { z } ) } \\ & { \mathbf { B } _ { g } = ( B _ { i + 1 / 2 , j , k } ^ { x } , B _ { i , j + 1 / 2 , k + 1 / 2 } ^ { y } , B _ { i , j + 1 / 2 , k + 1 / 2 } ^ { z } ) } \\ & { \bar { \mathbf { J } } _ { g } = ( \bar { J } _ { i , j + 1 / 2 , k + 1 / 2 } ^ { x } , \bar { J } _ { i + 1 / 2 , j , k + 1 / 2 } ^ { y } , \bar { J } _ { i + 1 / 2 , j + 1 / 2 , k } ^ { z } ) . } \end{array}
$$

The discrete operator $\nabla$ is a centered difference in space (second order accurate):

$$
\begin{array} { r l } & { \nabla \mathbf { f } _ { i , j , k } = ( ( f _ { i + 1 / 2 , j , k } - f _ { i - 1 / 2 , j , k } ) / \Delta x , ( f _ { i , j + 1 / 2 , k } - f _ { i , j - 1 / 2 , k } ) / \Delta y , } \\ & { ( f _ { i , j , k + 1 / 2 } - f _ { i , j , k - 1 / 2 } ) / \Delta z ) . } \end{array}
$$

The properties hold for the discrete operator $\nabla$ in the chosen spatial discretization:

$$
\nabla \cdot \nabla \times \mathbf { f } = 0 \qquad \nabla \cdot ( \mathbf { f } \times \mathbf { h } ) = \mathbf { h } \cdot ( \nabla \times \mathbf { f } ) - \mathbf { f } \cdot ( \nabla \times \mathbf { h } )
$$

The energy conserving PIC method is based on the concurrent solution of the coupled Equations 9 and 12 by a non-linear solver. The algorithm is summarized in Figure 1. The simulation is initialized first, setting the particles positions and electromagnetic field self consistently. At each PIC computational cycle, Equations 9 and 12 (non-linearly coupled by Equation 13) are solved by a JFNK solver. The governing equations have been time differenced with the implicit midpoint method, and solved in terms of $\overline { { \mathbf { v } } } _ { p }$ , $\mathbf { E } ^ { n + 1 }$ , and $\mathbf { B } ^ { n + 1 }$ . The spatial operators of Equations 12 are differenced in space on the Yee's lattice. Once $\overline { { \mathbf { v } } } _ { p }$ has been calculated with the solver, the new particle positions and velocities are simply updated with Equation 10.

![](images/eeedfb890191347844281622005ff28bc6d81ff2701fdde0607ce5706afd9b4c.jpg)  
Figure 1: The energy conserving PIC algorithm. After the simulation has been initialized, the computational cycle (non-linear solver stage, and particle update) is repeated.

# 2.1. The Electrostatic Limit

The electrostatic formulation of the energy conserving PIC method can be derived by considering an unmagnetized plasma. In this case, the particle average velocity Equation 9 simply reduces to:

$$
\bar { \mathbf { v } } _ { p } = \mathbf { v } _ { p } ^ { n } + \frac { q _ { s } } { 2 m _ { s } } \bar { \mathbf { E } } _ { p } \Delta t = \mathbf { v } _ { p } ^ { n } + \frac { q _ { s } } { 4 m _ { s } } ( \mathbf { E } _ { p } ^ { n + 1 } ( \bar { \mathbf { x } } _ { p } ) + \mathbf { E } _ { p } ^ { n } ( \bar { \mathbf { x } } _ { p } ) ) \Delta t .
$$

The evolution of the electric field is determined by the Ampere's law:

$$
\mathbf { E } _ { g } ^ { n + 1 } - \mathbf { E } _ { g } ^ { n } = - 4 \pi \bar { \mathbf { J } } _ { \mathbf { g } } \Delta t .
$$

# 2.2. Divergence Equations

Only the Faraday's and Ampere's equations were considered so far, while the two divergence equations of the system 11 were not taken in account. It easy to show that the equation $\nabla \cdot \mathbf { B } = 0$ is always satisfied if it is initially [1]. Moreover, the Gauss' law $\nabla \cdot \mathbf { E } = 4 \pi \rho$ equation is automatically satisfied by Equation 12, if the charge continuity equation $\partial \rho / \partial t + \nabla \cdot \mathbf { J } = 0$ holds true. In Particle-in-Cell methods, the charge continuity equation is not always satisfied because discrepancies between the interpolation of charge and current densities into the grid [1]. In fact, the definition of the current density in the energy conserving PIC method (Equation 13) does not satisfy the charge density continuity equation. In this case, the method is said not to conserve the charge. As pointed out in Ref.[15], the Gauss' law can be regarded as a conservation principle: it is not a strictly necessary equation for describing the evolution of electromagnetic fields, but its violation introduces numerical errors in the simulation and might lead to unphysical behavior of the simulated plasma. All the energy conserving PIC simulations, based on the non conservative current density definition, have been first initialized solving the Gauss' law, ensuring there is no error due to the violation of the charge continuity equation initially. Then, the charge conservation has been constantly checked. The error did not grow considerably or lead to unphysical behavior of the simulated plasma. The numerical error can be reduced using the pseudo current method as in Refs.[15, 16]. In this approach, $F _ { g } ^ { n }$ is defined as the violation of the Gauss' law in the cell $g$ at the time level $n$ , $F _ { g } ^ { n } = \nabla \cdot \mathbf { E } _ { g } ^ { n } - 4 \pi \rho _ { g } ^ { n }$ , and its gradient added to the Ampere's equation:

$$
\frac { \mathbf { E } _ { g } ^ { n + 1 } - \mathbf { E } _ { g } ^ { n } } { \Delta t } = c \nabla \times \bar { \mathbf { B } } _ { g } - 4 \pi \bar { \mathbf { J } } _ { \mathbf { g } } + 4 \pi d \nabla \bar { F } _ { g } ,
$$

where $d$ is a parameter that regulates the charge conservation, and $\bar { F } _ { g }$ is $1 / 2 ( F _ { g } ^ { n + 1 } + F _ { g } ^ { n } )$ (differently from Ref.[15], F is calculated at $n + 1 / 2$ and solved implicitly). Taking the divergence of Equation 19:

$$
\frac { F _ { g } ^ { n + 1 } - F _ { g } ^ { n } } { d \Delta t } - \nabla ^ { 2 } \bar { F } _ { g } = - ( \frac { \rho _ { g } ^ { n + 1 } + \rho _ { g } ^ { n } } { \Delta t } + \nabla \cdot \bar { \mathbf { J } } _ { \mathbf { g } } )
$$

The left side of the equation above is the heat equation. If $F _ { g } ^ { \prime }$ is fixed to zero at the boundaries initially and at each simulation time step, the error due to non conservation of charge diffuses away with a rate determined by the parameter $d$ . Figure 2 shows the error due to the violation of Gauss' law using $d$ equal to 0 (no pseudo current correction), and to $0 . 1 c ^ { 2 } / \omega _ { p e }$ , $0 . 5 c ^ { 2 } / \omega _ { p e }$ in a simulation of the two-stream instability. It has been found that the pseudo current correction decreases the error related to the non conservation of charge. No major differences appeared in runs with and without the pseudo current correction: the plasma instabilities under study started at the same time and presented the same growth rate. It is clear from Figure 2 that the numerical error builds up slowly in the case of no pseudo current correction (blue line).

![](images/7ae5b3a203e2f4d1a3aef5b29c71ef24142f3929dc71a8af955c91a2454a129b.jpg)  
Figure 2: Maximum violation of Gauss' law divided by the maximum charge density in a two-stream instability simulation. The blue line represents the numerical error without the pseudo current correction, while the red and green lines show its evolution with a pseudo current correction with $d = 0 . 1 c ^ { 2 } / \omega _ { p e }$ and $0 . 5 c ^ { 2 } / \omega _ { p e }$ respectively.

# 3. Energy Conservation

The discretized equations of the proposed PIC method have the property of conserving exactly the total energy. In fact, the variation of the magnetic

and electric fields energies ( $W _ { E }$ , $W _ { B }$ ) during the time step $\Delta t$ is:

$$
\Delta ( W _ { E } + W _ { B } ) = \frac { 1 } { 8 \pi } \sum _ { g } ^ { N _ { g } } \left( ( \mathbf { E } _ { g } ^ { n + 1 } ) ^ { 2 } - ( \mathbf { E } _ { g } ^ { n } ) ^ { 2 } + ( \mathbf { B } _ { g } ^ { n + 1 } ) ^ { 2 } - ( \mathbf { B } _ { g } ^ { n } ) ^ { 2 } \right) V _ { g }
$$

Substituting Equations 12 in the formula above:

$$
\begin{array} { r l } & { \Delta ( W _ { E } + W _ { B } ) = \frac { 1 } { 4 \pi } \sum _ { g } ^ { N _ { g } } \Big ( \bar { \mathbf { E } } _ { g } \cdot ( c \nabla \times \bar { \mathbf { B } } _ { g } - 4 \pi \bar { \mathbf { J } } _ { g } ) - \bar { \mathbf { B } } _ { g } \cdot ( c \nabla \times \bar { \mathbf { E } } _ { g } \Delta t ) \Big ) \Delta t V _ { g } } \\ & { = \sum _ { g } ^ { N _ { g } } \Big ( - \bar { \mathbf { E } } _ { g } \cdot \bar { \mathbf { J } } _ { g } - \nabla \cdot \bar { \mathbf { S } } _ { g } \Big ) \Delta t V _ { g } . } \end{array}
$$

The vector identities 16, that hold in the chosen spatial discretization, have been used. The first term of the equation above represents the work of the field on the particles, while the second term $\bar { \mathbf { S } } _ { g } = c / 4 \pi ( \bar { \mathbf { E } } _ { g } { \times } \bar { \mathbf { B } } _ { g } )$ is the average Poynting flux. Its contribution to the energy variation is zero in an isolated system. If the expression for the average current density ${ \bar { \mathbf { J } } } _ { g }$ (Equation 13) is substituted in the formula above, then:

$$
\begin{array} { r l } & { \boldsymbol { \Delta } ( W _ { E } + W _ { B } ) = \sum _ { g } ^ { N _ { g } } \Big ( - \bar { \mathbf { E } } _ { g } \cdot ( \sum _ { s } ^ { n _ { s } } \sum _ { p } ^ { N _ { s } } q _ { s } \bar { \mathbf { v } } _ { p } W ( \mathbf { x } _ { g } - \bar { \mathbf { x } } _ { p } ) / V _ { g } ) \Big ) \boldsymbol { \Delta } t V _ { g } } \\ & { = - \sum _ { s } ^ { n _ { s } } \sum _ { p } ^ { N _ { s } } q _ { s } \boldsymbol { \Delta } t \bar { \mathbf { v } } _ { p } \cdot \sum _ { g } ^ { N _ { g } } \bar { \mathbf { E } } _ { g } W ( \mathbf { x } _ { g } - \bar { \mathbf { x } } _ { p } ) } \\ & { = - \sum _ { s } ^ { n _ { s } } \sum _ { p } ^ { N _ { s } } m _ { s } \bar { \mathbf { v } } _ { p } \cdot ( \mathbf { v } _ { p } ^ { n + 1 } - \mathbf { v } _ { p } ^ { n + 1 } - \sum _ { g } ^ { N _ { g } } \bar { \mathbf { v } } _ { p } / c \times \bar { \mathbf { B } } _ { g } ) } \\ & { = - \sum _ { s } ^ { n _ { s } } \sum _ { p } ^ { N _ { s } } 1 / 2 m _ { s } ( ( \mathbf { v } _ { p } ^ { n + 1 } ) ^ { 2 } - ( \mathbf { v } _ { p } ^ { n } ) ^ { 2 } ) . } \end{array}
$$

Thus the numerical scheme presented in Section 2 conserves the total energy. It must be pointed out that a different definition of $\overline { { \mathbf { v } } } _ { p }$ in Eq. 9, using $\mathbf { B } _ { p } ^ { n }$ instead of $\bar { \mathbf { B } } _ { p }$ , still leads to the conservation energy because the magnetic field does no work on a charged particle. In addition, the conservation of energy is a consequence of the definition of the current density, Equation 13, and PIC schemes with different techniques for the computation of current density to ensure charge conservation do not conserve the total energy.

# 4. Momentum Conservation

The energy conserving PIC method does not conserve momentum, as reported by Langdon in Ref.[5] for other "energy conserving" PIC schemes.The non conservation of momentum in energy conserving PIC schemes is due to spurious particle self-forces arising from the non smoothness of the current density deposition to the grid. The relevance of the particle self-forces depends largely on the interpolation functions in use, on the number of particles per cell, and on the grid spacing [5]. Self-forces might trigger a macroscopic instability, as shown in Ref.[5]. For instance in Figure 3, two cold electron beams, composed by 10000 particles (in blue in the Figure 3) are moving initially at opposite velocities $\pm 0 . 2 c$ in a one dimensional 2.053 $c / \omega _ { p e }$ long periodic domain with simulation time step $\Delta t = 0 . 5$ and 64 grid cells. The simulation is completed with the energy conserving PIC method. An aliasing instability develops phase space vortices (red dots in Figure 3) later in time. A parametric study has been completed varying the grid spacing, the time step and $v _ { c }$ , the characteristic velocity of the plasma (e.g. the thermal velocity in a Maxwellian plasma or the drift velocity in beams). It has been found empirically that this aliasing instability disappears if the condition

![](images/a84cf148e250b355d18b978d2191138185cf40080fa056521208450cd806399d.jpg)  
Figure 3: Phase space of simulation of two cold counter-streaming electron beams with the energy conserving PIC method. The blue dots represent particles in the initial config$t = 1 3 \omega _ { p e } ^ { - 1 }$ . An instability due to the non conservation of momentum is visible as red vortices in the phase space.

$$
v _ { c } \Delta t / \Delta x < 1 . 5
$$

is satisfied. Thus, a modified Courant-Friederics condition, that restricts the maximum particle motion to one and half cell per time step, must be satisfied. The velocity $v _ { c }$ can be connected to the fast electron scales, but no significant limitation arises when one desires the simulation of ion dynamics. It is important to note that, given a $v _ { c }$ that is related to the electron scales (e.g., electron thermal velocity), the ratio $\Delta t / \Delta x$ can be still adjusted to satisfy Eq. 24, to follow the slower ion scales. In fact, such numerical condition does not preclude the possibility of simulations with large $\Delta t$ , but it only requires that the grid spacing is chosen accordingly to satisfy the condition 24. It must be pointed out that this constraint, that arises from the non conservation of momentum, is very similar to the one of the implicit moment PIC method [3] and for this reason the energy conserving PIC method allows simulation time steps that are as large as the one allowed by the implicit moment PIC method.

# 5. Numerical Stability

The numerical stability of PIC methods can be determined by studying the plasma numerical dispersion relation [17], following the examples of Langdon in Ref.[8], and of Brackbill and Forslund in Ref.[3] for the electrostatic limit. In this approach, the particle equation of motion is linearized, and the electric field is assumed to have an $\exp ( i \omega t )$ dependence. The linearized equation of motion is Fourier transformed in $x$ , and then the perturbed charge density and the plasma susceptibility are calculated. Following this approach, the numerical dispersion of the energy conserving PIC method results:

$$
1 - ( \frac { \omega _ { p e } \Delta t } { 2 } ) ^ { 2 } \int _ { - \infty } ^ { + \infty } f _ { 0 } ( { \bf { v } } ) \frac { \cos ( ( \omega - { \bf { k } } \cdot { \bf { v } } ) \frac { \Delta t } { 2 } ) } { \sin ^ { 2 } ( ( \omega - { \bf { k } } \cdot { \bf { v } } ) \frac { \Delta t } { 2 } ) } d { \bf { v } } = 0 ,
$$

where $f _ { 0 } ( \mathbf { v } )$ is the equilibrium distribution function, $\omega _ { p e } = \sqrt { 4 \pi n _ { e } q _ { e } ^ { 2 } / m _ { e } }$ is the plasma frequency, $n _ { e }$ is the plasma density, and $\mathbf { k }$ is the wave vector. In the case of cold plasma with $f _ { 0 } ( \mathbf { v } ) = \delta ( \mathbf { v } )$ , the numerical dispersion relation reduces to:

$$
\tan ( \frac { \omega \Delta t } { 2 } ) \sin ( \frac { \omega \Delta t } { 2 } ) = ( \frac { \omega _ { p e } \Delta t } { 2 } ) ^ { 2 } .
$$

The roots of the dispersion relation are always real and therefore neither exponential growth nor damping is present at any choice of $\Delta t$ . For this reason, the numerical scheme is linearly unconditional stable. For comparison, the numerical dispersion relation of explicit PIC method for a cold plasma leads to exponential growth for the well known condition $\omega _ { p e } \Delta t > 2$ [1]. The dispersion analysis that includes grid effects is not carried out in the present paper, and it will part of a future work.

In the case of implicit moment PIC methods, a $\theta$ parameter is introduced in the numerical scheme to decenter in time the discretization of the field equations [3]. The quantities $q$ at time level $\theta$ are defined as $( 1 - \theta ) q ^ { n } + \theta q ^ { n + 1 }$ . The numerical dispersion relation of the implicit moment PIC method has growing solutions for $\theta < 1 / 2$ , damped solutions for $\theta > 1 / 2$ , and neither damping nor growth for $\theta = 1 / 2$ [3]. For $\theta > 0 . 5$ , the implicit moment PIC method damps high frequency waves, that are not resolved by the time step, as shown in Ref.[3]. On the contrary, in energy conserving PIC simulations unresolved waves are not artificially damped and hold over the simulation. To compare the behavior of other fully implicit PIC schemes with the energy conserving PIC method, a $\theta$ parameter has been introduced in the numerical scheme as follows:

$$
\begin{array} { r l } & { { \tilde { \bf v } } _ { p } = { \bf v } _ { p } ^ { n } + \frac { q _ { s } \Delta t } { 2 m _ { s } } { \bf E } _ { p } ^ { \theta } } \\ & { { \bf v } _ { p } ^ { \theta } = { \tilde { \bf v } } _ { p } + \frac { q _ { s } \Delta t } { 2 m _ { s } c } ( { \tilde { \bf v } } _ { p } \times { \bf B } _ { p } ^ { \theta } + \frac { q _ { s } \Delta t } { 2 m _ { s } c } ( { \tilde { \bf v } } _ { p } \cdot { \bf B } _ { p } ^ { \theta } ) { \bf B } _ { p } ^ { \theta } ) / ( 1 + \frac { q _ { s } ^ { 2 } \Delta t ^ { 2 } } { 4 m _ { s } ^ { 2 } c ^ { 2 } } B _ { p } ^ { \theta 2 } ) } \\ & { { \bf E } _ { g } ^ { n + 1 } - { \bf E } _ { g } ^ { n } = c \nabla \times { \bf B } _ { g } ^ { \theta } \Delta t - 4 \pi { \bf J } _ { g } ^ { \theta } \Delta t } \\ & { { \bf B } _ { g } ^ { n + 1 } - { \bf B } _ { g } ^ { n } = - c \nabla \times { \bf E } _ { g } ^ { \theta } \Delta t . } \end{array}
$$

The equations above are solved concurrently by a JFNK solver. For $\theta = 0 . 5$ , the method reduces to energy conserving PIC scheme, while for $\theta > 0 . 5$ does not conserve the total energy. In the latter case, the method artificially cools the plasma damping unresolved waves. This is clear from Figure 4, where the numerical dispersion relation of a plasma undergoing Weibel instability [18] is shown. The numerical dispersion relation has been calculated by applying the fast Fourier transform in space and in time to the $z$ component of the magnetic field, as shown in Ref. [19]. The Weibel instability can be seen as a vertical red line in both panels of Figure 4. In addition, thermal noise is visible as a red line, that is diagonal for small $k$ and becomes horizontal for high $k$ , only in the left panel of Figure 4 for the energy conserving PIC code with $\theta = 0 . 5$ . The radiation field is due to aliasing errors, a feature of the quadratically conserving schemes [20]. On the contrary, the radiation field noise is damped in the fully implicit PIC simulation with $\theta = 0 . 6$ and not visible in the right panel.

![](images/c2331a23683f827b25ab0b82f0322b7eaa0804085324fefdc9a676471db52011.jpg)  
Figure 4: Dispersion relation in the Weibel instability simulation for $\theta = 0 . 5$ (energy conserving PIC code) and for $\theta = 0 . 6$ (fully implicit PIC code with numerical damping). The Weibel instability is visible as vertical red lines in both simulations. The noise in the radiation field (a red line that is diagonal for small $k$ and becomes horizontal for high $k$ ) is present in the energy conserving simulation (left panel, $\theta = 0 . 5$ ). Instead it is damped, and therefore not present, in the $\theta = 0 . 6$ PIC simulation (right panel).

# 6. Implementation

An energy conserving PIC code has been developed in the Matlab/Octave programming language. For implementation simplicity, the code is 1D3V [1] with uniform grid: the space is one dimensional while the particle velocities, the electric and magnetic fields have three components. The extension of the code to the three dimensional case is straight-forward.

After the simulation has been initialized setting self-consistently particle quantities and electromagnetic fields, and ensuring that the Gauss' law is satisfied initially, two steps are completed at each computational cycle, as shown in Figure 1:

1. The values of the dependent variables h he N $\mathbf { E } _ { g } ^ { n + 1 }$ b vn. $\mathbf { v } _ { p } ^ { \pi }$ , ${ \bf E } _ { a } ^ { n }$ $\mathbf { B } _ { g } ^ { n + 1 }$ and and $\mathbf { B } _ { g } ^ { n }$ from the pre- $\overline { { \mathbf { v } } } _ { p }$ are detervious computational cycle, $\overline { { \mathbf { v } } } _ { p }$ , $\mathbf { E } _ { g } ^ { n + 1 }$ and $\mathbf { \Delta B } _ { g } ^ { n + 1 }$ are calculated solving concurrently Equations 9 and 12 by a JFNK solver [21].

Th $\mathbf { X } _ { p } ^ { n + 1 }$ , $\mathbf { v } _ { p } ^ { n + 1 }$ are updated with Equations 10.

The JFNK method solves the non-linear system $\mathbf { G } ( \mathbf { x } ) = \mathbf { 0 }$ iteratively by computing successive linear systems:

$$
\frac { \partial \mathbf { G } ( \mathbf { x } ) } { \partial \mathbf { x } } \bigg \rvert _ { i } \delta \mathbf { x } _ { i } = - \mathbf { G } ( \mathbf { x } _ { i } ) .
$$

the solution guess $\mathbf { x _ { i } }$ at the iteration $i$ of the initial non-linear system ${ \bf G } ( { \bf x } ) =$ $\mathbf { 0 }$ is calculated as:

$$
\mathbf { x } _ { i + 1 } = \mathbf { x } _ { i } + \delta \mathbf { x } _ { i } .
$$

The solution of the linear system 28 and the solution update (Eq. 29) compose the Newton iteration, that stops when:

$$
\begin{array} { r } { \| \textbf { G } ( \mathbf { x } _ { i } ) \ \| < \epsilon _ { a } + \epsilon _ { r } \ \| \textbf { G } ( \mathbf { x } _ { 0 } ) \ \| , } \end{array}
$$

where $\| \cdot \|$ is the Euclidean norm, and $\epsilon _ { a }$ and $\epsilon _ { r }$ are the absolute and relative error tolerance values. The successive linear systems 28 are solved iteratively by a Krylov method, the Generalized Minimal Residual (GMRes) solver [21] in the present study. The Krylov method convergence is adjusted at each Newton iteration as follow:

$$
\parallel \left. \frac { \partial \mathbf { G } ( \mathbf { x } ) } { \partial \mathbf { x } } \right| _ { i } \delta \mathbf { x } _ { i } + \mathbf { G } ( \mathbf { x } _ { i } ) \parallel < \zeta _ { i } \parallel \mathbf { G } ( \mathbf { x } _ { i } ) \parallel ,
$$

where $\zeta _ { i }$ is the inexact Newton parameter [21]. The number of iterations of the GMRes solver are called Krylov iterations. The Jacobian $\frac { \partial \mathbf { G } ( \mathbf { x } ) } { \partial \mathbf { x } }$ is not calculated directly, but instead the Gateaux derivative is used to compute:

$$
\left. \frac { \partial \mathbf { G } ( \mathbf { x } ) } { \partial \mathbf { x } } \right| _ { i } \delta \mathbf { x } _ { i } = \operatorname* { l i m } _ { \epsilon \to 0 } \frac { \mathbf { G } ( \mathbf { x } _ { i } + \epsilon \delta \mathbf { x } _ { i } ) - \mathbf { G } ( \mathbf { x } _ { i } ) } { \epsilon } ,
$$

where $\epsilon$ is a small but finite number. Because the the Jacobian does not need to be formed and calculated explicitly, the method is said " Jacobian-free".

A guess of the particle average velocity, of the new electric and magnetic fields (the dependent variables) is given initially as vector $\mathbf { x } _ { K R }$ to the GMRes solver. A better estimate of these values is calculated by minimizing through successive solver iterations the residual $\mathbf { r }$ (the difference between the known term $\mathbf { b }$ and $A ( \mathbf { x } _ { K R } )$ , the non-linear system $A$ applied to the solution guess $\mathbf { x } _ { K R }$ )

$$
\mathbf { r } = \mathbf { b } - A ( \mathbf { x } _ { K R } ) .
$$

A function where the residual $\mathbf { r }$ is calculated, must be provided to the JFNK solver. The residual is computed by solving the Equations 9 and 12 for the problem unknowns, $\bar { \mathbf { v } } _ { p } ~ \mathbf { E } _ { g } ^ { n + 1 }$ n Bn+1,  v:

1. Given a $\overline { { \mathbf { v } } } _ { p }$ estimate in $\mathbf { x } _ { K R }$ , $\bar { \bf x }$ is calculated as $\overline { { \mathbf { v } } } _ { p } \Delta t / 2$ , and ${ \bar { \mathbf { J } } } _ { g }$ is computed with Equation 13.   
2. Given $\mathbf { E } _ { g } ^ { n + 1 }$ and $\mathbf { B } _ { g } ^ { n + 1 }$ estimates in $\mathbf { x } _ { K R }$ , $\bar { \mathbf { E } } _ { p }$ and $\bar { \mathbf { B } } _ { p }$ are calculated with Equation 6.   
3. The residual $\mathbf { r }$ is computed with Equations 9 and 12.

The solution of the particle equations of motion and field equations, and the current deposition are completed at each Krylov iteration.

A skeleton version of the Matlab/Octave code for the electrostatic limit with electrons and motionless ions is presented in Appendices A and B to show the simplicity of the proposed PIC method. The software implementation of the Newton Krylov GMRes solver is from the Kelley's textbook [21] and available at the website [22]. In all the simulations, the solver maximum number of Newton and Krylov iterations is set to 40, and the inexact Newton parameter $\zeta _ { i }$ is determined by the Eisenstat-Walker formula [21]. The Eisenstat-Walker parameter is chosen as 0.9. In the energy conserving PIC method, smaller error tolerance values lead to simulations with increased energy conservation. This is clearly visible in Figure 5, where the energy history of the same simulation of the two-stream instability with different absolute and relative solver tolerance values is plotted. On the contrary, it has been found that decreasing the error tolerances does not have any effect in the conservation of the momentum.

# 7. Simulation Results

The energy conserving PIC codes have been tested throughly. The goal of these tests is first to verify the new PIC method through comparison of the simulation results with analytical theory, and second to show the exact energy conservation. In this paper, the energy conserving PIC code is first run in the electrostatic formulation for the problems of the finite grid and two-stream instabilities, and then in the fully electromagnetic case for the Weibel instability test.

# ABS $=$ 1 RE $=$ 1 $=$ $=$ $=$ $=$

![](images/b00d76d0c61aa88716fd835affa5bd10679ee7a8ed078e8268bfca7d0c411373.jpg)  
Figure 5: Comparison of energy histories in a two-stream instability simulation with different absolute and relative error tolerance values. Smaller tolerance errors lead to an increased energy conservation.

# 7.1. Finite Grid Instability

An aliasing instability, called finite grid instability, arises in explicit momentum conserving PIC methods when the simulation grid spacing is approximately two times larger than Debye length $\Lambda _ { D }$ [1, 19]. The finite-grid instability heats non physically the plasma, until the Debye length reaches a value comparable to half the grid spacing. Because this instability introduces numerical heat in the system, it appears as a macroscopic violation of the energy conservation.

To test the energy conserving electrostatic PIC method against the finite grid instability, a Maxwellian plasma is initialized with thermal velocity $v _ { t h e } = 0 . 2 c$ in a simulation box $5 0 \pi c / \omega _ { p e }$ long with 64 grid cells and 50,000 particles. The Debye length $\Lambda _ { D } = v _ { t h e } / \omega _ { p e } = 0 . 2 c / \omega _ { p e }$ results approximately ten times smaller than the grid spacing $\Delta x = 2 . 4 5 c / \omega _ { p e }$ . This geometrical set-up leads to the finite grid instability if an explicit momentum conserving PIC code is used. The simulation evolves over 200 computational cycles with time step equal to $0 . 5 \omega _ { p e } ^ { - 1 }$ . Therefore, the numerical constraint $v _ { t h e } \Delta t / \Delta x = 0 . 2 4 5 < 1 . 5$ is satisfied and numerical instability does not arise because of the non conservation of momentum. The absolute and relative solver error tolerance values are both set to $1 0 ^ { - 8 }$ . In addition, a simulation with an explicit momentum conserving PIC code, starting from the same initial configuration, has been run to compare the results. Figure 6 represents the phase space (each dot represents a particle in the position-velocity space) of the system at $t = 1 0 0 \omega _ { p e } ^ { - 1 }$ forthe eney cnering (red ot) and explt momentum conserving (blue dots) methods. The finite grid instability in the simulation with the explicit PIC code is visible from the blue peaks in the phase space. Instead the plasma retains the initial Maxwellian distribution in the energy conserving PIC simulation. Figure 7 shows the total energy history for the two simulations. The finite grid instability produces a 2% energy increase in the explicit PIC simulation, while the variation with the energy conserving PIC method is very low, $1 0 ^ { - 7 } \%$ . Because the energy conserving PIC does not undergo finite grid instability, and does not require to resolve the Debye length, the proposed PIC method is well suited for simulations with large domains and/or few grid points. A complete linear and non-linear analysis of the finite grid instability in the energy conserving PIC has not been carried out and it will be a topic of a future work.

![](images/8876fdead317d6372065fdce4ea4b55b36d3dbe9a3fb3ed8aee1004a7791d88b.jpg)  
Figure 6: Phase space plot of a Maxwellian plasma simulation at $t = 1 0 0 \omega _ { p e } ^ { - 1 }$ with the energy conserving (red dots) and explicit (blue dots) PIC codes. The finite grid instability appears as peaks of the electron distribution in the phase space in the explicit PIC simulation, while is not present in the energy conserving PIC simulation, where the plasma retains the initial Maxwellian distribution.

![](images/27f95abf5aed9e90144ed2329695f5d5d37583661793a8ca0f0f75007ebfc2cc.jpg)  
Figure 7: Total energy history of a Maxwellian plasma simulation with the energy conserving PIC code (red line and left $y$ axis) and with the explicit momentum conserving PIC code (blue line and right $y$ axis). Energy is in $n _ { e } m _ { e } c ^ { 2 } / 2$ units. Because of the finite grid instability, the total energy in the explicit PIC simulation increases $2 \%$ . On the contrary, the energy variation is limited to $1 0 ^ { - 7 } \%$ in the case of the energy conserving PIC simulation.

# 7.2. Two-stream Instability

The two-stream instability is an important phenomenon occurring in space physics, in the injection systems for nuclear fusion machine, and in particle accelerators [23]. In this problem, two electron beams move initially in opposite directions. The two beams extinguish as result of the beams instability. A simulation of the two-stream instability has been completed with the energy conserving PIC code. The drift velocity of the counter-streaming electron cold beams is $\pm 0 . 2 c$ ; the simulation box is $2 . 0 5 3 ~ c / \omega _ { p e }$ long with 64 grid points and periodic boundaries. The number of electrons and ions is 200,000. The charge to mass ratio for electrons and ions is one and 1836 repectively. The smulation timetep $\Delta t$ is 0.1 $\omega _ { p e } ^ { - 1 }$ . This set-up leads to $v _ { c } \Delta t / \Delta x = 0 . 6 2 3 < 1 . 5$ , and thus the instability due to the non conservation of the total momentum does not occur. The absolute and relative tolerance are both set to $1 0 ^ { - 8 }$ .

The linear theory predicts a growth rate of instability for the spectral component $k = 1 \omega _ { p e } / c$ equal to 0.35355 $\omega _ { p e }$ [23] in this system configuration. Figure 8 shows an excellent agreement between the linear theory in red line and the simulation results in blue line in the linear stage of the instability. The energy variation for an explicit momentum conserving, and an energy conserving PIC code, simulating the two-stream instability are compared in Figure 9. The total energy in the explicit PIC code shows an approximately $5 \%$ variation, while the variation is limited to $1 0 ^ { - 4 } ~ \%$ in the energy conserving PIC code. An increased energy conservation can be achieved, decreasing the error tolerance values, as shown in Figure 5.

# 7.3. Weibel Instability

The Weibel instability is triggered by an anisotropic temperature in the plasma [18, 24]. The instability converts the plasma temperature anisotropy into magnetic field relaxing the initial particle distribution function to an isotropic Maxwellian. Because plasma temperature anisotropies are very common in laboratory, space and astrophysical plasmas, the Weibel instability has been thoroughly studied with PIC methods [24]. The energy conserving PIC code has been tested in a simulation box $2 \pi c / \omega _ { p e }$ long, with 64 grid points and periodic boundaries. The time step $\Delta t$ is 0.25 $\omega _ { p e } ^ { - 1 }$ and simulation is run for 400 computational cycles. 100,000 electrons are initialized with uniform distribution in space and bi-Maxwellian distribution with thermal velocity $v _ { t h y } = ~ 0 . 4$ , and anisotropy $a = ( v _ { t h y } ^ { 2 } / v _ { t h x } ^ { 2 } - 1 ) = 1 5$ . Thus, the condition $v _ { t h x } \Delta t / \Delta x = 1 . 0 2$ to avoid the aliasing instability due to the non conservation of momentum, is satisfied. The mass ratio between ions and electron is 1836, and ions are initialized with same temperature of electrons. The electric and magnetic field are initially zero. The solver absolute and relative error tolerance values are both set to $1 0 ^ { - 5 }$ .

![](images/9a610d56a439d4f7160cee3860e5cebea612c94e624cf40cdfd18ed213db1482.jpg)  
Figure 8: Comparison between linear theory and energy conserving PIC simulation of the two-stream instability. The $k = 1 \omega _ { p e } / c$ spectral component of the electric field (in blue color) grows as predicted by the linear theory (red line). The electric field is in $\sqrt { 4 \pi n _ { e } m _ { e } c ^ { 2 } }$ units.

![](images/9d2d128681eec691eaaba973375e1009a09a6e2fe5c538db8a0ebe328e243f72.jpg)  
Figure 9: Comparison of the total energy history of the explicit momentum conserving (blue line, right $y$ axis) and energy conserving PIC (red line, left $y$ axis) codes for the two-stream instability. Energy is in $n _ { e } m _ { e } c ^ { 2 } / 2$ units. The plot shows $5 \%$ and $1 0 ^ { - 4 } \ \%$ variations for the explicit and energy conserving PIC codes respectively.

The linear theory predicts a growth rate of the $B _ { z }$ component for the spectral component $k = 1 \omega _ { p e } / c$ equal to $0 . 2 2 \omega _ { p e }$ [23]. Figure 10 compares the Weibel instability simulation results with the analytical calculations. Linear theory and simulation results are in good agreement in the linear regime of the instability. However, the growth of the magnetic field is not purely exponential as predicted [23], but presents an oscillation in time. This oscillation is caused by the radiation field noise [24, 25]. As shown in Section 5 and in the left panel of Figure 4, the numerical dispersion relation of the Weibel instability with the energy conserving PIC code shows the presence of waves due to the thermal noise [25]. Figure 11 shows the total energy and momentum history in a Weibel instability simulation with the energy conserving PIC method. The total energy variation is within $1 0 ^ { - 3 } ~ \%$ . Momentum is not conserved and oscillates between $- 0 . 0 2 m _ { e } c$ and $0 . 0 2 m _ { e } c$ .

# 8. Performance Results

A study of the computational performance of the proposed PIC scheme has been completed. A Maxwellian plasma, composed of electron and a background ions, is simulated with an electrostatic energy conserving PIC in a $6 . 4 c / \omega _ { p e }$ long box over a period of 1000 cycles. The other simulation settings vary to perform a parametric study. Table 1 shows the average of Newton and Krylov iterations, and the execution time for different number of grid points, time step, number of particles, absolute and relative error tolerance values. The tests have been completed on a 2.4 GHz Intel Core Duo, 2 GB RAM memory, Mac OS X 10.6.6 using the Matlab 7.5 and the code in Appendices A and B. The computational performance of the energy conserving weakly depends on the number of cells: increasing the number of grid point from 64 to 512 leads to a 20% computational time increase. The time step has a similar effect on performance also. The simulation with $d t =$ $0 . 8 \omega _ { p e } ^ { - 1 }$ takes only an additional 33% computational time of the simulation with $d t = 0 . 1 \omega _ { p e } ^ { - 1 }$ . Instead, the computational performance strongly depends on the number of particles. In fact, doubling the number of particles leads to doubling the computational time. In addition, increasing the number of particles reduces the statistical noise, and consequently the convergence in the Newton step. The decrease JFNK error tolerance increases the degree of conservation energy as shown in Figure 5 at the cost of an increased number of Newton iterations and computational time.

![](images/e1fcb74333408b252dc052f4fcde4b6fb599ab27d48319d9c5a9975c816dbbd1.jpg)  
Figure 10: Comparison between the spectral component $k = 1 \omega _ { p e } / c$ of $B _ { z }$ and linear theory. Simulation and analytical results are in good agreement in the linear regime of the Weibel instability. The magnetic field is in $\sqrt { 4 \pi n _ { e } m _ { e } c ^ { 2 } }$ units.

![](images/c1951c80cfb9ee99119ceb165a85e02da26201b0a9c20fd59ba8e1252cf35266.jpg)  
Figure 11: Total energy and momentum histories in the Weibel instability simulation. Energy and momentum are in $n _ { e } m _ { e } c ^ { 2 } / 2$ and $m _ { e } c$ units. The left red $y$ and the right blue $y$ axes correspond to the total energy and momentum respectively. The total energy is conserved within $1 0 ^ { - 3 } ~ \%$ variation, while the momentum is not conserved and oscillates between $- 0 . 0 2 m _ { e } c$ and $0 . 0 2 m _ { e } c$ .

Table 1: Average number of Newton and Krylov (GMRes) iterations and execution time for different number grid points, time step, number of particles, and solver error tolerances $( \epsilon _ { a } , \epsilon _ { r } )$ for a simulation of Maxwellian plasma with an electrostatic energy conserving PIC code.   

<table><tr><td>Ng</td><td>dt (ωpe 1)</td><td>Ns</td><td>ar</td><td>Newton</td><td>Krylov</td><td>Time (s)</td></tr><tr><td>64</td><td>0.1</td><td>100000</td><td>10-7, 10-7</td><td>3.62</td><td>2.32</td><td>341.86</td></tr><tr><td>128</td><td>0.1</td><td>100000</td><td>10-7, 10-7</td><td>4.04</td><td>2.33</td><td>401.28</td></tr><tr><td>256</td><td>0.1</td><td>100000</td><td>10-7, 10-7</td><td>3.98</td><td>2.50</td><td>416.66</td></tr><tr><td>512</td><td>0.1</td><td>100000</td><td>10-7, 10-7</td><td>4.02</td><td>2.49</td><td>413.18</td></tr><tr><td>64</td><td>0.2</td><td>100000</td><td>10-7, 10-7</td><td>4.15</td><td>2.45</td><td>408.38</td></tr><tr><td>64</td><td>0.4</td><td>100000</td><td>10-7, 10-7</td><td>4.03</td><td>2.58</td><td>424.98</td></tr><tr><td>64</td><td>0.8</td><td>100000</td><td>10-7, 10-7</td><td>4.31</td><td>2.69</td><td>455.57</td></tr><tr><td>64</td><td>0.1</td><td>200000</td><td>10-7, 10-7</td><td>3.25</td><td>2.37</td><td>617.13</td></tr><tr><td>64</td><td>0.1</td><td>400000</td><td>10-7, 10-7</td><td>2.91</td><td>2.39</td><td>1,144.90</td></tr><tr><td>64</td><td>0.1</td><td>800000</td><td>10-7, 10-7</td><td>2.79</td><td>2.42</td><td>2,194.00</td></tr><tr><td>64</td><td>0.1</td><td>100000</td><td>10-8, 10-8</td><td>4.55</td><td>2.34</td><td>402.19</td></tr><tr><td>64</td><td>0.1</td><td>100000</td><td>10-9, 10-9</td><td>4.82</td><td>2.42</td><td>433.95</td></tr><tr><td>64</td><td>0.1</td><td>100000</td><td>10-10&#x27; 10-10 ,</td><td>4.92</td><td>2.43</td><td>447.01</td></tr></table>

# 8.1. Kinetic Enslavement

The main disadvantage of the energy conserving PIC method is that it requires the solution of a very large system whose size increases with the number of particles. The number of particles is easily more million in a typical PIC simulation, leading to matrix to be inverted whose rank size is of the order of million. To reduce the size of this matrix, it is possible to use a technique, called kinetic enslavement, following Ref.[26]. In this method, the unknowns of the problem are only the value of the electric magnetic field on the grid points at the new time level and the JFNK solver computes only the field Equations 12. The particle equations of motion are calculated by a Newton-Raphson method and embedded in the field solver as function evaluations [26]. Figure 12 shows the the computational cycle of the energy conserving PIC method with kinetic enslavement. A guess of the electromagnetic field ( $\tilde { \bf E }$ , $\tilde { \mathbf { B } }$ in Figure 12) is provided at each Newton iteration. Particle positions and velocities ( $\ddot { \bf x }$ , $\tilde { \bf v }$ in Figure 12) are computed consistently with the the electromagnetic field guess by the Newton-Raphson method.

An energy conserving PIC code with kinetic enslavement technique has been developed to test its effectiveness. Figures 13, 14 show a comparisons of the solver iterations for the energy conserving PIC code with and without the kinetic enslavement. The plots represent the number of Newton and average Krylov iterations per Newton step in the simulation of the two stream instability with an electrostatic energy conserving PIC code with and without kinetic enslavement. The two stream instability is simulated for 500 cycles, with $d t ~ = ~ 0 . 1$ , 64 grid points and 1000 particles. The Newton-Raphson convergence tolerance is set to $1 0 ^ { - 4 }$ , while the JFNK error tolerance values are set to $1 0 ^ { - 7 }$ . The size of the systems is reduced from $1 0 6 4 \times 1 0 6 4$ (energy conserving PIC) to $6 4 \times 6 4$ (energy conserving PIC with kinetic enslavement). In addition, it is clear from Figure 13 that the use of kinetic enslavement technique reduces the number of Newton iterations: the average number of Newton iterations is 1.91 and 4.33 in the simulation with and without kinetic enslavement method. The number of Krylov iterations remains almost unchanged in the two cases (Figure 14). It must be noted that an average of

![](images/f8bfac69c491491dccc65b59c80651b8ad0c95b398f4826549d3abe056fac733.jpg)  
Figure 12: Computation cycle for the energy conserving PIC with kinetic enslavement. Particle positions and velocities, $\tilde { \bf x }$ , $\tilde { \mathbf { v } }$ , are calculated consistently with the fields, $\tilde { \bf E }$ , $\tilde { \mathbf { B } }$ by Newton-Raphson method at each JFNK iteration.

3.95 Newton-Raphson iterations have been completed at each solver iteration.

# 9. Conclusions

The algorithm, the properties, the implementation, the simulation and performance results of the energy conserving Particle-in-Cell method have been presented. The proposed PIC method has been tested against the finite grid, two-stream and Weibel instabilities to prove the algorithm correctness and its property of conserving exactly the total energy. The method is based on a non conservative definition of the current density. The numerical error due to the violation of the Gauss' law built up slowly in all the completed simulations and did not affect the results. The new method does not conserve the momentum, and a condition on the maximum number of cells a particle can cross per time step, must be satisfied to avoid an aliasing instability. The new PIC scheme is based on the implicit discretization of the governing equation, and therefore results linearly unconditionally stable.

The energy conserving PIC method is a fully implicit PIC method[10, 9], where the particle average velocities and the electromagnetic field are calculated self-consistently in the JFNK solver to preserve the system total energy. The performance of the fully implicit PIC methods, and a comparison between fully implicit and implicit moment PIC methods in terms of efficiency and required computational resources, have been presented in Refs.[9, 10]. It has been shown in this paper that the energy conserving PIC performance critically depends on two factors: the number of computational particles and the solver error tolerance values. In fact, the computational time increases linearly with the number of particles, and it is rather insensitive to the number of grid points and the time step. In addition, a smaller error tolerance value leads to a larger number of iterations, and therefore to a larger computational time. The use of the kinetic enslavement technique [26] reduces the size of problem matrix and has the beneficial effect of decreasing the number of Newton iterations.

![](images/d49526d864ffc07b93795120c375a1e50854ee22a89cf780ea13a12e23dfb45f.jpg)  
Figure 13: Comparison of number of Newton iterations in the energy conserving PIC code with and without kinetic enslavement for the two stream instability test.

![](images/f320501b6faee95dfe99b5e08005d69b7d4d6440bdcd31fdc40e29a0172aa137.jpg)  
Figure 14: Comparison of the average number of Krylov iterations per Newton step in the energy conserving PIC code with and without kinetic enslavement for the two stream instability test.

# Acknowledgment

The authors are grateful to Gianni Coppa for his clever implementation of explicit electrostatic PIC in Matlab/Octave as well as for the fruitful discussions on the mathematical foundations of the PIC method. The authors are also grateful to Jerry Brackbill for the stimulating exchanges of ideas on the implicit PIC method. The present work is supported by the Onderzoekfonds KU Leuven (Research Fund KU Leuven) and by the European Commissions Seventh Framework Programme (FP7/2007-2013) under the grant agreement no. 218816 (soteria-space.eu).

# References

[1] C. K. Birdsall, A. B. Langdon, Plasma Physics via Computer Simulation, McGraw-Hill, New York, 1985.   
[2] R. Hockney, J. Eastwood, Computer Simulation Using Particles, Taylor & Francis, 1988.   
[3] J. Brackbill, D. Forslund, Simulation of low-frequency, electromagnetic phenomena in plasmas, in: J. U. Brackbill, B. Cohen (Eds.), Multiple Time Scales, Academic Press Orlando, 1985.   
[4] H. R. Lewis, Energy-conserving numerical approximations for Vlasov plasmas, Journal of Computational Physics 6 (1) (1970) 136 - 141.   
[5] A. B. Langdon, 'Energy-conserving' plasma simulation algorithms, Journal of Computational Physics 12 (2) (1973) 247  268.   
[6] C. Kelley, Solving Nonlinear Equations with Newton's Method, Fundamentals of Algorithms, SIAM, Philadelphia, 2003.   
[7] D. A. Knoll, D. E. Keyes, Jacobian-free Newton-Krylov methods: a survey of approaches and applications, Journal of Computational Physics 193 (2) (2004) 357397.   
[8] A. Langdon, Analysis of the time integration in plasma simulation, Journal of Computational Physics 30 (2) (1979) 202  221.   
[9] S. Markidis, Development of implicit kinetic simulation methods, and their application to ion beam propagation in current and future neutralized drift compression experiments., Ph.D. thesis, University of Illinois at Urbana-Champaign (2010).   
[10] H. Kim, L. Chacon, G. Lapenta, Fully implicit particle in cell algorithm, in: Bull. Am. Phys. Soc, no. 2913, Denver, 2005.   
[11] J. Dawson, Particle simulation of plasmas, Reviews of Modern Physics 55 (1983) 403447.   
[12] E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration, Springer-Verlag, 2002.   
[13] J. M. Hyman, M. Shashkov, Mimetic discretizations for maxwell's equations, Journal of Computational Physics 151 (2) (1999) 881  909.   
[14] K. Yee, Numerical solution of initial boundary value problems involving maxwell's equations in isotropic media, IEEE Transactions on Antennas and Propagation 14 (1966) 302307.   
[15] B. Marder, A method for incorporating Gauss' law into electromagnetic PIC codes, Journal of Computational Physics 68 (1) (1987) 48  55.   
[16] A. B. Langdon, On enforcing Gauss' law in electromagnetic particlein-cell codes, Computer Physics Communications 70 (3) (1992) 447 - 450.   
[17] E. Lindman, Dispersion relation for computer-simulated plasmas, Journal of Computational Physics 5 (1970) 1322.   
[18] E. S. Weibel, Spontaneously growing transverse waves in a plasma due to an anisotropic velocity distribution, Phys. Rev. Lett. 2 (3) (1959) 8384.   
[19] H. Matsumoto, Y. Omura, Particle simulation of electromagnetic waves and its application to space plasmas, 1985, pp. 43102.   
[20] A. Arakawa, Y. Hsu, Energy conserving and potential-enstrophy dissipating schemes for the shallow water equations, Monthly Weather Review 118 (1990) 19601969.   
[21] C. Kelley, Iterative Methods for Linear and Nonlinear Equations, Frontiers in Applied Mathematics, SIAM, Philadelphia, 1995.   
[22] C. Kelley, Matlab codes from iterative methods for linear and nonlinear equations, http://www.siam.org/books/kelley/fr16/matlabcode. php (August 2010).   
[23] N. A. Krall, A. W. Trivelpiece, Principles of Plasma Physics, McGrawHill New York, 1973.   
[24] R. L. Morse, C. W. Nielson, Numerical simulation of the Weibel instability in one and two dimensions, Physics of Fluids 14 (4) (1971) 830840.   
[25] A. B. Langdon, Some electromagnetic plasma simulation methods and their noise properties, Physics of Fluids 15 (6) (1972) 11491151.   
[26] W. T. Taitano, Development of a Jacobian-free-Newton-Krylov method with kinetic enslavement to implicitly solve Vlasov-Poisson system in plasma physics (2010).   
[27] S. Markidis, G. Lapenta, Matlab energy conserving PIC code, https: //perswww.kuleuven.be/\~u0070495/Site/EC.html (March 2011).

# Appendix A. ECpicES.m

A skeleton version of the energy conserving PIC in Matlab/Octave programming language is here presented for the electrostatic limit (Section 2.1) with a plasma of electrons and motionless background ions [27]. The energy conserving PIC code requires additional files (nsolgm.m, fdgmres.m, givapp.m, and dirder.m), that are available at the Kelley's textbook website [22].

After the initial parameters are set and the electric field is calculated solving the Gauss' law to ensure the charge conservation, the average particle velocities and the new electric field are calculated at line 52, and the new particle positions and velocities are updated at lines 55 and 56 at each computational cycle.

1 % Energy Conserving PIC code   
2 global L; global dx; global NG; global DT; global N; global WP;   
3 global QM; global Q; global rho_back;   
4 global x0; global v0; global E0;   
5   
6 % simulation parameters   
7 $\mathrm { L } { = } 2 { * } \mathrm { p } \mathrm { i }$ /3.0600; % Simulation box length   
8 $\mathrm { D T } { = } 0 . 1$ ; % time step   
9 $\mathrm { N T } { = } 8 0 0$ ; % number of computational cycles   
10 $\mathrm { N G } { = } 1 2 8$ ; % number of cells   
11 $\mathrm { N } { = } 5 0 0 0 0 0$ ; % number of particles   
12 $\mathrm { W P } { = } 1$ ; % plasma frequency   
13 $\mathrm { Q M = - 1 }$ ; % electron charge to mass ratio   
14 $\mathrm { V 0 = 0 . 2 }$ ; % beam velocity   
15 $\mathrm { V T } { = } 0 . 0$ ; % thermal velocity   
16 $\mathrm { t o l ~ } = \ [ 1 \mathrm { E } \mathrm { - } 7 , \ 1 \mathrm { E } \mathrm { - } 7 ]$ ; % absolute and relative error tolerance   
17 $\mathrm { { d x } \mathrm { { = } L / N G } }$ ; % grid spacing   
18 QWP^2/(QM\*N/L); % computational particle charge   
19 rho_back= Q\*N/L; % background ion charge density   
20 histEnergy $\qquad = \ [ \ ]$ ; %total energy history   
21 % two-stream instability   
22 % initial particle positions   
23 $\mathrm { x 0 { = } 1 }$ inspace $\left( 0 , \mathrm { L - L } / \mathrm { N } , \mathrm { N } \right)$ '; % uniform in space   
24 % initial particle velocities   
25 v0=VT\*randn(N,1); % two counterstreaming beams   
26 $\mathrm { { p m } = [ 1 { : N } ] }$ '; $\mathrm { p m } { = } 1 { - } 2 { * } \mathrm { m o d } \left( \mathrm { p m } , 2 \right)$ ; $\mathrm { \ v 0 = v 0 + p m . * V }$ 0;   
27   
28 % Perturbation   
29 $\mathrm { X P 1 = 1 }$ ; $\mathrm { V 1 { = } 0 . 0 }$ ; mode $= 1$ ;   
30 $\scriptstyle \mathrm { v 0 = v 0 + V 1 * \mathbf { s } i n }$ (2\*pi\*x0/L\*mode);   
31 x0=x0+XP1\*(L/N)\* sin(2\*pi\*x0/L\*mode ) ;   
32 $\mathrm { { o u t } = ( x 0 < 0 ) ; ~ x 0 \left( \ o u t \right) = x 0 \left( \ o u t \right) + L }$ ;   
33 $\mathrm { \ o u t { = } } ( \mathrm { x 0 { > } { = } L } ) ,$ ; $\mathrm { { x 0 } ( \mathrm { { o u t } ) { = x 0 } ( \mathrm { { o u t } ) \mathrm { { - L } } } } }$ ;   
34   
35 % calculate E0, satisfying the Gauss' Law   
36 solving the Poisson equation   
37 $\mathrm { p } { = } 1 { : } \mathrm { N }$ ; $\mathrm { { p = } [ p ^ { \Delta } \mathrm { { p } ] } }$ ; un=ones(NG-1,1);   
38 Poisson $=$ $\mathbf { \bar { s } } \mathbf { p } \mathbf { d i a g s } \left( \left[ \mathrm { u n } \ - 2 \mathrm { * u n } \ \mathrm { u n } \right] , \left[ - 1 \ 0 1 \right] , \mathrm { N G - 1 , N G - 1 } \right) ;$   
39 g1=floor $( \mathrm { x 0 } / \mathrm { d x } - . 5 ) + 1$ b $\mathbf { g } = [ \mathrm { g } 1 ; \mathrm { g } 1 + 1 ]$ ;   
40 $\mathrm { f r a z 1 { = } 1 { - } a b s \left( x 0 / d x { - } g 1 + . 5 \right) ; }$ fr a z =[ f r a z l;1 − f r a z 1 ];   
41 $\mathrm { { o u t } = ( g < 1 ) ; g ( o u t ) = g ( o u t ) + N G ; }$   
42 $\mathrm { o u t } = \left( \mathrm { g } > \mathrm { N G } \right) ; \mathrm { g } \left( \mathrm { o u t } \right) = \mathrm { g } \left( \mathrm { o u t } \right) - \mathrm { N G }$   
43 mat=sparse(p,g,fraz ,N,NG);   
44 rho=full $( \mathrm { ( Q / d x ) * s u m ( m a t ) } ) ^ { \prime } + \mathrm { r h o \_ b a c k } ;$   
45 $\mathrm { P h i { = } P o i s s o n \setminus ( - r h o \left( 1 : N G { - } 1 \right) * d x \setminus 2 \sqrt { \Omega } ) }$ ; $\mathrm { { P h i = } [ \mathrm { { P h i } } ; 0 ] }$ ;   
46 $\begin{array} { r l } { \mathrm { E 0 } } & { = \left( \left[ \mathrm { P h i } \left( \mathrm { N G } \right) ; \mathrm { \nabla { P h i } } \left( 1 : \mathrm { N G } - 1 \right) \right] - \left[ \mathrm { P h i } \left( 2 : \mathrm { N G } \right) ; \mathrm { P h i } \left( 1 \right) \right] \right) / \left( 2 * \mathrm { d } \mathbf { x } \right) ; } \end{array}$   
47   
48 for $\mathrm { i t } = 1 { : } \mathrm { N T }$   
49 % start computational cycle   
50 xkrylov $= \ [ \mathrm { v 0 } ; \ \mathrm { E 0 } ]$ ;   
51 % Newton Krylov GMRes solver   
52 [sol, it_hist , ierr] = nsolgm(xkrylov, 'residueEC',tol);   
53 v_average $=$ sol(1:N);   
54 $\%$ update particle positions and velocities   
55 $\mathrm { \Delta v 0 ~ = ~ 2 * v }$ average - v0;   
56 $\mathrm { x 0 ~ = ~ \mathrm { x 0 ~ + } }$ v_average $\ast \mathrm { D T }$ ;   
57 % check if particle are out of the periodic boundaries   
58 $\scriptstyle \mathrm { o u t } = \left( \mathrm { x 0 } < 0 \right)$ ; x0(out) ${ \boldsymbol { \mathbf { \mathit { \sigma } } } } = \mathbf { \boldsymbol { \mathbf { \mathit { \mathbf { x } } } } } \mathbf { \boldsymbol { 0 } }$ (out) $+ \mathrm { L }$ ;   
59 out= $( \mathrm { x 0 } \mathrm { > } \mathrm { = } \mathrm { L }$ ) ; x0(out) ${ \boldsymbol { \mathbf { \mathit { \sigma } } } } = \mathbf { \boldsymbol { \mathbf { \mathit { \mathbf { x } } } } } \mathbf { \boldsymbol { 0 } }$ (out)−L;   
60 % new electric field   
61 $\mathrm { E 0 ~ = }$ sol $( \mathrm { N } + 1 ) \colon ( \mathrm { N } + \mathrm { N G } ) )$ ;   
62 % calculate the total energy   
63 Etot $\mathbf { \Sigma } = \mathbf { \Sigma } 0 . 5 * \mathbf { a b s } \left( \mathbf { Q } \right) * \mathbf { s u m } ( \mathbf { v } 0 \cdot \mathbf { \hat { \mu } } 2 ) \mathbf { \Sigma } + \mathbf { \Sigma } 0 . 5 * \mathbf { s u m } ( \mathbf { E } 0 \cdot \mathbf { \hat { \mu } } 2 ) * \mathbf { d x } ;$   
64 % save the total energy   
65 histEnergy $=$ [histEnergy Etot];   
66 % end computational cycle   
67 end

# Appendix B. residueEC.m

The Newton Krylov GMRes solver (nsolgm.m) requires the definition of a residue function (residueEC.m), where the discretized equations of the energy conserving PIC method are formulated. The following Matlab/Octave code presents the residue function for the energy conserving PIC code. The particle average velocity equations 17 are defined at line 22, while the field equations 18 at line 30.

% residual calculation for the EC PIC   
2 function res $=$ residueEC(xkrylov)   
  
global L; global dx; global NG; global DT; global N; global WP;   
5 global QM; global Q; global rho_back;   
6 global $\mathrm { x 0 }$ ; global v0; global E0;   
8 % calculate the x at $n { + } \bar { \bot } \bar { \angle }$ time level   
9 x_average $\qquad = \ \mathrm { \bf { x } } 0$ $^ +$ xkrylov $\left( 1 : \mathrm { N } \right) * \mathrm { D T } / 2$ ;   
10 % check if particle are out $o f$ the periodic boundaries   
11 out=(x_average $< 0$ );x_average(out) $= \mathrm { x }$ -average(out) $+ \mathrm { L }$ ;   
12 $\scriptstyle \mathrm { o u t } = ( \mathrm { x \_ a v e r a g e } > = \mathrm { L } )$ );x_average(out) $=$ x_average(out)−L;   
13 % interpolation   
14 $\operatorname { p = 1 } { : \mathrm { N } } ; \operatorname { p = } [ \operatorname { p } \ \operatorname { p } ]$ ; $\mathrm { g 1 = }$ floor(x_average $^ { \prime } \mathrm { d } \mathrm { x } - . 5 ) + 1$ ; $\mathbf { g } = [ \mathrm { g } 1 ; \mathrm { g } 1 + 1 ]$ ;   
15 fraz $1 { = } 1 -$ abs(x_average $( 1 : \mathrm { N } ) / \mathrm { d } \mathrm { x } { - } \mathrm { g } 1 + . 5 )$ ; fra $z =$ [(fraz1);1 − fr az1 ];   
16 $\mathrm { { o u t } = ( g < 1 ) ; g ( \ o u t ) = g ( \ o u t ) \Sigma + N G }$ ;   
17 $\mathrm { o u t } = \left( \mathrm { g } \mathrm { > N G } \right) ; \mathrm { g } \left( \mathrm { o u t } \right) = \mathrm { g } \left( \mathrm { o u t } \right) - \mathrm { N C }$ ;   
18 mat=sparse(p,g, fraz ,N,NG);   
19   
20 res = zeros $( \mathrm { N } + \mathrm { N G } , 1 )$ ;   
21 % average velocity residual   
22 res $( 1 { : } \mathrm { N } , 1 ) =$ xkrylov (1:N)−v0 −0.25\*mat\*QM\*(E0+xkrylov ((N+1):(N+NG)))\*DT;   
23   
24 % calculate the average $\boldsymbol { \mathcal { T } }$   
25 fraz $=$ [( fr az 1 ) . $^ *$ xkrylov (1:N);(1 −fraz1 ). $^ *$ xkrylov (1:N)];   
26 mat=sparse(p,g, fraz ,N,NG);   
27 $\mathrm { ~ J ~ } =$ full $\mathrm { ' ( Q / d x ) * s u m ( m a t ) ) } ^ { \prime }$ ;   
28   
29 electric field residual   
30 res((N+1):(N+NG)) $=$ xkrylov ((N+1):(N+NG))−E0+J\*DT;