# Mie Scattering in the Macroscopic Response and the Photonic Bands of Metamaterials

Lucila Juárez\*1, Bernardo S. Mendoza', and W. Luis Mochán²2

1Centro de Investigaciones en ptica, Lomas del Bosque 115, Lomas del Campestre, 37150 Lón, Guanajuato, Méic. vcaiv Morelos, México

Key words: metamaterials; homogenization; Mie resonances; non-local optics \* Corresponding author: lucilajr@icf.unam.mx

We present a general approach for the numerical calculation of the effective dielectric tensor of metamaterials and show that our formalism can be used to study metamaterials beyond the long wavelength limit. We consider a system composed of high refractive index cylindric inclusions and show that our method reproduces the Mie resonant features and photonic band structure obtained from a multiple scattering approach, hence opening the possibility to study arbitrarily complex geometries for the design of resonance-based negative refractive metamaterials at optical wavelengths.

1 Introduction Metamaterials are usually composed of periodic arrangements of micro- or nano- structures form ing ordered patterns designed to control the propagation of light. They can display exotic optical properties with interesting aplications, such as a negative refractive index which can be achieved exploiting underlying resonant features of the microstructure [1,2,3,4,5]. Resonant behavior has been vastly studied in split-ring-resonator (SRR) structures, which are generally composed of metallic ringlike structures which give rise to a significant magnetic response when excited by an external inhomogeneous electric field. However, metamaterials based on SRRs are not well suited to visible frequencies, due to size limitations for its fabrication and high losses of the metallic components [6,7,8]. Recently, high refractive index inclusions have been investigated as potential components for the fabrication of low-loss metamaterials displaying strong magnetic properties at optical frequencies [6,9, 10, 11, 12, 13, 14, 15]. In dielectric materials, the displacement current increases with increasing permittivity. In high refractive index particles the displacement current can become large and give rise to the so-called Mie resonances at wavelength comparable with the size of the particles. Thus, small highindex inclusions which display strong electric and magnetic multipole Mie resonances in the optical region can be potentially used as resonators for the fabrication of low loss negative refraction metamaterials. Negative refraction has been reported for microstructure geometries as simple as cylinders [12, 16]. Indeed, the number and frequencies of the resonant modes depend strongly on the size and geometry of the particles. Therefore, the generation and interference of such modes can be controlled by manipulating the composition, size and shape of the particles, allowing a wide range of possibilities for the design of optical metamaterials. However, general methods to compute the effective response of metamaterials of arbitrary shape and composition are often limited to computationally expensive numerical approaches or approximations. If the wavelength of the incident light $\lambda$ is large compared to the microstructure of the metamaterial, its response can be described by an effective macroscopic dielectric function efficiently computed within the so-called long wavelength limit [17, 18]. Nevertheless, the description of Mie scattering lies by definition outside the validity range of the long wavelength limit generally used to study metamaterials. When the incident field varies in space on a length scale comparable with the microscale of of the metamaterial, the effects of retardation and non locality become important. It has been shown that even when retardation effects are important, the system may be characterized by a macroscopic dielectric response [19], which must be described by a non-local tensor $\epsilon ( r , r ^ { \prime } , t , t ^ { \prime } )$ [20]. This non-local retarded response results in a frequency $\omega$ and wavevector $\boldsymbol { k }$ dependence of the effective dielectric tensor in Fourier space which leads to the constitutive equation $D ( \omega , \pmb { k } ) =$ $\epsilon ( \omega , \boldsymbol { k } ) E ( \omega , \boldsymbol { k } )$ .

In this paper, we use a general efficient formalism for the numerical computation of the effective dielectric response of metamaterials at arbitrary wavelenghts. We consider the case of a metamaterial composed by high index dielectric cylinders and show that our approach can reproduce the analytical results obtained from Mie theory at waveleghts comparable with the size of the particle, thereby demonstrating that our method can be used to study metamaterials based on Mie resonances.

The structure of the paper is the following. In Section 2 we present our theory for the calculation of the electromagnetic response of the metamaterial. First in Subsection 2.1 we develop an efficient computational method based on the calculation of the macroscopic dielectric response through a recursive procedure. This method is applilcable to arbitrary materials and geometries. In order to interpret and test its results, in Subsection 2.2 we develop a multiple scattering approach applicable to an array of dielectric cylinders. Results on the dielectric response of high-index dielectric cylinders are presented in Section 3, compared to the results of the multiple scattering approach and interpreted in terms of coupled Mie resonances. Finally, our conclusions are presented in Section 4.

# 2 Theoretical methods

2.1 Macroscopic response Recall that any physical vector field has a longitudinal part which can be obtained by the projection operator $\bar { \mathcal { P } } _ { L } = \nabla \nabla ^ { - 2 } \nabla$ where $\nabla ^ { - 2 }$ is the inverse of the Laplacian operator. In reciprocal space this longitudinal projector can be written as $\bar { \mathcal { P } } _ { L } =$ $\bar { k } k / k ^ { 2 }$ , where $\boldsymbol { k }$ is the wavevector of magnitude $k$ . Furthermore, for a periodic system with Bravais lattice $\{ R \} =$ $\{ \sum _ { i } ^ { D } n _ { i } d _ { i } \}$ where $n _ { i }$ are inegers, $\mathbf { \mathbf { { \alpha } } } d _ { i }$ are primitive lattice vectors and is the number of dimensions, the electric field can be expressed using Bloch's theorem as

$$
E _ { \boldsymbol { k } } ( \boldsymbol { r } ) = \sum _ { \boldsymbol { G } } E _ { \boldsymbol { G } } \exp ( \boldsymbol { i } ( \boldsymbol { k } + \boldsymbol { G } ) \cdot \boldsymbol { r } ) ,
$$

where $E _ { G }$ is the amplitude of a plane wave with wavevector $\pm G$ , with $G$ a vector of the reciprocal lattice defined by $\exp ( i G \cdot R ) = 1$ and $\boldsymbol { k }$ a vector within the first Brillouin zone. Here, the long wavelength limit corresponds to $k \ \ll \ G$ [17]. Notice that the field components with wavevectors $\pm G$ fluctuate over distances of the order $d _ { i }$ , except for the term $G = 0$ which we identify as the average of the field. With this definition we can write an expression for the average operator as

$$
\hat { \mathcal P } _ { G G ^ { \prime } } ^ { a } = \delta _ { G 0 } \delta _ { G ^ { \prime } 0 } .
$$

Consider a binary metamaterial composed of a periodic lattice of microstructures of arbitrary shape embedded in an homogeneous medium. We define the structure function $\boldsymbol { B }$ which contains the information on the shape of the inclusions

$$
\begin{array} { r } { B ( r ) = \left\{ \begin{array} { l l } { 1 , r \in B , } \\ { 0 , r \notin B . } \end{array} \right. } \end{array}
$$

For a system of two components, say, a host $( A )$ and inclusions $( B )$ , having each well defined dielectric functions $\epsilon _ { A }$ and $\epsilon _ { B }$ , the microscopic dielectric function can be written as

$$
\epsilon ( r ) = \left\{ \epsilon _ { A } , r \in A , \right.
$$

which we abbreviate in terms of the structure function as $\begin{array} { r } { \epsilon ( r ) = \frac { \epsilon _ { A } } { u } \left( u - B ( r ) \right) } \end{array}$ ,where $\begin{array} { r } { u = 1 / ( 1 - \frac { \epsilon _ { A } } { \epsilon _ { B } } ) } \end{array}$ is the spectral variable.

We begin by writing the wave equation for the electric field at frequency $\omega$ using the free wavevector $q \equiv \omega / c$ as

$$
\hat { \mathcal { W } } { \cal E } = \left( \hat { \epsilon } + \frac { 1 } { q ^ { 2 } } \nabla ^ { 2 } \hat { \mathcal { P } } _ { T } \right) { \cal E } = \frac { 4 \pi } { i \omega } { \cal J } _ { \mathrm { e x t } } ,
$$

which we solve formally as

$$
{ \pmb E } = \frac { 4 \pi } { i \omega } \hat { \mathcal { W } } ^ { - 1 } J _ { \mathrm { e x t } } ,
$$

where we have introduced a wave operator $\hat { \mathcal W }$ written in terms of the transverse projector $\hat { \mathcal { P } } _ { T } = \mathbb { 1 } - \hat { \mathcal { P } } _ { L }$ .Here, we identify the average of the field $\pmb { { \cal E } } ^ { a } = \hat { \mathcal { P } } ^ { a } \pmb { { \cal E } }$ in Eq. (6) with the macroscopic field ${ \pmb E } ^ { M }$ .

Since $\mathbf { \delta } J _ { \mathrm { e x t } }$ is an external current, it does not have spa$J _ { \mathrm { e x t } } ^ { a } = J _ { \mathrm { e x t } }$ Thus, we average Eq. (6) to obtain

$$
{ \pmb E } _ { M } = \frac { 4 \pi } { i \omega } \ \hat { \mathcal { W } } _ { M } ^ { - 1 } { \pmb J } _ { \mathrm { e x t } } ,
$$

where the macroscopic wave operator $\hat { \mathcal { W } } _ { M }$ is given by $\hat { \mathcal { W } } _ { M } ^ { - 1 } = \hat { \mathcal { W } } _ { a a } ^ { - 1 } = \hat { \mathcal { P } } ^ { a } \hat { \mathcal { W } } ^ { - 1 } \hat { \mathcal { P } } ^ { a }$ i i n i theavge of the inverse of the microscopic wave operator.

Substitution of the permittivity tensor in terms of the structure function Eq. (3) and the spectral variable $u$ , leads to the wave-operator

$$
\hat { \mathcal { W } } = \frac { \epsilon _ { A } } { u } \left( u - \hat { \mathcal { B } } \right) + \frac { 1 } { q ^ { 2 } } \nabla ^ { 2 } \hat { \mathcal { P } } _ { T } ,
$$

which we rewrite as

$$
\hat { \mathcal { W } } = \frac { \epsilon _ { A } } { u } \left( u \hat { g } ^ { - 1 } - \hat { \mathcal { B } } \right)
$$

by introducing a metric operator

$$
\hat { g } = \left( \mathbf { 1 } + \hat { \mathcal { P } } _ { T } \frac { \nabla ^ { 2 } } { q ^ { 2 } \epsilon _ { A } } \right) ^ { - 1 } .
$$

Inverting the wave operator and taking the average we obtain

$$
\hat { \mathcal { W } } _ { M } ^ { - 1 } = \hat { \mathcal { W } } _ { a a } ^ { - 1 } = \frac { u } { \epsilon _ { A } } \hat { g } _ { a a } \left( u - \hat { \mathcal { B } } \hat { g } \right) _ { a a } ^ { - 1 } ,
$$

where we used the fact that the metric doesn't couple average to fluctuating fields.

Finally, we extract the macroscopic dielectric tensor from the corresponding wave operator

$$
\epsilon ^ { M } ( \omega , \boldsymbol { k } ) = \frac { 1 } { q ^ { 2 } } ( k ^ { 2 } \mathbf { 1 } - \boldsymbol { k } \boldsymbol { k } ) + \mathcal { W } ^ { M } ( \omega , \boldsymbol { k } ) ,
$$

where we used the explicit transverse projector for a plane wave with wavevector $\boldsymbol { k }$ .

In order to compute the macroscopic dielectric tensor, we begin by calculating $( u - B \hat { g } ) _ { a a } ^ { - 1 }$ in Eq. (11). To that end, we first notice that the operator $\hat { B }$ is Hermitian in the usual sense. The operator $\hat { g }$ would also be Hermitian if the response of medium A is dissipationless, i.e., if $\epsilon _ { A }$ is real. Nevertheless, the product $\hat { B g }$ is not Hermitian. We notice, however, that the product ${ \cal B } \hat { g }$ becomes Hermitian by redefining the internal product between two states using $\hat { g }$ as a metric tensor. Thus, we define the $g$ -product of two states $| \phi \rangle$ and $| \psi \rangle$ as $\scriptstyle ( \phi | \psi )$ , where

$$
\left( \phi | \psi \right) \equiv \left. \phi | \hat { g } | \psi \right. ,
$$

and $\langle \dots | \dots \rangle$ is the usual Hermitian scalar product. With this definition it is clear that

$$
( \phi | ( \hat { B } \hat { g } | \psi ) = \langle \phi | \hat { g } \hat { B } \hat { g } | \psi \rangle = \langle \psi | \hat { g } \hat { B } \hat { g } | \phi \rangle ^ { \ast } = ( \psi | \hat { B } \hat { g } | \phi ) ^ { \ast } ,
$$

so that $\hat { B g }$ is indeed Hermitian under the product $\left( \ldots | \ldots \right)$ and we may borrow computational methods developed for quantum mechanical calculations.

We choose an initial state $\left| 0 \right. = b _ { 0 } ^ { - 1 } \left| p \right.$ where $| p \rangle$ corresponds to a plane wave of frequency $\omega$ , wavevector $\boldsymbol { k }$ and polarization $\hat { e }$ , normalized as $\langle p | p \rangle = 1$ under the conventional internal product, and $b _ { 0 }$ is chosen such that the state $| 0 \rangle$ is $g$ -normalized, $( 0 | 0 ) = g _ { 0 } = \pm 1 \ /$ Notice that as $\hat { g }$ is not positive definite, we should allow for negative norms. We also define $| - 1 \rangle = 0$ .Following the Haydock recursive scheme [21], new states can be generated by repeatedly applying the Hermitian operator

$$
\hat { B } \hat { g } \left| n \right. = b _ { n + 1 } \left| n + 1 \right. + a _ { n } \left| n \right. + b _ { n } g _ { n } g _ { n - 1 } \left| n - 1 \right.
$$

where the real coefficients $a _ { n }$ , $b _ { n }$ and $g _ { n }$ are obtained by imposing the orthonormality condition

$$
( n | m ) = \langle n | \hat { g } | m \rangle = g _ { n } \delta _ { n m }
$$

and $g _ { n } = \pm 1$ Thus the operator $\hat { B g }$ has a tridiagonal representation in the basis $| n \rangle$ , which allows us to express the

operator $\left( u - \hat { B } \hat { g } \right)$ as

$$
\left( u - \hat { B } \hat { g } \right) =
$$

$$
\left( { \begin{array} { c c c c } { u - a _ { 0 } - b _ { 1 } g _ { 1 } g _ { 0 } } & { 0 } & { 0 } & { \ldots } \\ { - b _ { 1 } } & { u - a _ { 1 } } & { - b _ { 2 } g _ { 2 } g _ { 1 } } & { 0 } \\ { 0 } & { - b _ { 2 } } & { u - a _ { 2 } } & { - b _ { 3 } g _ { 3 } g _ { 2 } } \\ { 0 } & { 0 } & { \ddots } & { \ddots } \\ { \vdots } & { } & { } & { } \\ { \vdots } & { } & { } & { } \end{array} } \right)
$$

Finally, we have to invert and average the operator in Eq. (17). We recall that the average (Eq. (2)) is given in terms of a projection into our starting state $| p \rangle$ , i.e., the zeroth row and column element of the inverse operator which may be found for Eq. (17) in the form of a continued fraction

$$
\hat { e } \cdot ( { \pmb { \mathscr W } } _ { M } ( \omega , \pmb k ) ) ^ { - 1 } \cdot \hat { e } =
$$

$$
\frac { u } { \epsilon _ { A } } \frac { g _ { 0 } b _ { 0 } ^ { 2 } } { u - a _ { 0 } - \frac { g _ { 0 } g _ { 1 } b _ { 1 } ^ { 2 } } { u - a _ { 1 } - \frac { g _ { 1 } g _ { 2 } b _ { 2 } ^ { 2 } } { u - a _ { 2 } - \frac { g _ { 2 } g _ { 3 } b _ { 3 } ^ { 2 } } { \ddots } } } . }
$$

Choosing different independent polarizations $\hat { e }$ for the initial state $| p \rangle$ , one can compute all the independent projections of the inverse of the wave tensor. The result is then substituted in Eq. (12) to obtain the fully retarded macroscopic dielectric tensor. Further details on the method and its implementation can be found in Refs. [19,22]. We remark that the procedure above may be performed for two phase systems of arbitrary geometry and composition as long as one of them is dissipationless.

2.2 Scattering approach In this section we follow Ref. [23] to compute the solution of the multiple scattering problem of a finite array of dielectric cylinders. Consider first a single infinitely long dielectric cylinder of radius $R$ and refractive index $n$ , standing in vacuum with its axis aligned to the $\hat { z }$ axis. An incident field polarized on the $x - y$ plane is applied. The magnetic field $\pmb { H } = ( 0 , 0 , \varPhi )$ is taken parallel to the axis of the cylinder, and satisfies scalar Helmholtz equations of the form

$$
\begin{array} { r l r } { \displaystyle \frac { 1 } { r } \frac { \partial } { \partial r } \left( r \frac { \partial \phi ^ { \beta } ( r , \theta ) } { \partial r } \right) + \frac { 1 } { r ^ { 2 } } \frac { \partial ^ { 2 } \phi ^ { \beta } ( r , \theta ) } { \partial ^ { 2 } \theta } } & { } & \\ { + \left( n ^ { \beta } ( \omega ) \frac { \omega } { c } \right) ^ { 2 } \phi ^ { \beta } ( r , \theta ) = 0 } \end{array}
$$

for each frequency $\omega$ , where $n ^ { \beta } = n$ or 1 is the refractive index of the region $\beta = I , O$ , inside or outside of the cylinder, respectively. Solutions of Eq. (19) can be obtained as

the products $\varPhi _ { l } ^ { \beta } ( \kappa ^ { \beta } ) \exp ( i l \theta )$ , where $l = 0 , \pm 1 , \pm 2 , \ldots$ and $\boldsymbol { \varPhi } _ { l } ^ { \beta }$ solves the Bessel differential equation

$$
\kappa ^ { \beta } \frac { d } { d \kappa ^ { \beta } } \left( \kappa ^ { \beta } \frac { d \varPhi _ { l } ^ { \beta } } { d \kappa ^ { \beta } } \right) + \left( ( \kappa ^ { \beta } ) ^ { 2 } - l ^ { 2 } \right) \varPhi _ { l } ^ { \beta } = 0 ,
$$

where $\kappa ^ { \beta } = n ^ { \beta } q r$ .The general solution of Eq. (20) inside (I) and outside (O) the cylinder can be written in terms of the Bessel functions of the first and second kind $J _ { l }$ and $Y _ { l }$ as

$$
\begin{array} { l } { { \phi ^ { I } ( r , \theta ) = \displaystyle \sum _ { l } c _ { l } J _ { l } ( n q r ) \exp ( i l \theta ) , \ ~ } } \\ { { \phi ^ { 0 } ( r , \theta ) = \displaystyle \sum _ { l } \left[ a _ { l } J _ { l } ( q r ) + b _ { l } H _ { l } ( q r ) \right] \exp ( i l \theta ) , \ ~ } } \end{array}
$$

where we have chosen the outgoing Hankel functions $H _ { l } =$ $J _ { l } + i Y _ { l }$ as the scattered field. The coefficents $a _ { l }$ describe the incident field and $b _ { l }$ and $c _ { l }$ are to be determined by the boundary conditions. These are the continuity of $H$ and the continuity of the component of the electric field $\pmb { { \cal E } } =$ $( i / \epsilon q ) \nabla \times \pmb { H }$ parallel to the interface, $E _ { \theta }$ , at $r = R$ . As

$$
\begin{array} { l } { { \displaystyle E _ { \theta } ^ { I } ( r , \theta ) = \frac { i } { n } \sum _ { l } c _ { l } J _ { l } ^ { \prime } ( n q r ) \exp ( i l \theta ) } } \\ { { \displaystyle E _ { \theta } ^ { 0 } ( r , \theta ) = i \sum _ { l } [ a _ { l } J _ { l } ^ { \prime } ( q r ) + b _ { l } H _ { l } ^ { \prime } ( q r ) ] \exp ( i l \theta ) } } \end{array}
$$

where $J _ { l } ^ { \prime }$ and $H _ { l } ^ { \prime }$ denote the derivatives of the Bessel and Hankel functions with respect to their arguments, then

$$
\begin{array} { r } { c _ { l } J _ { l } ( n q R ) = a _ { l } J _ { l } ( q R ) + b _ { l } H _ { l } ( q R ) , } \\ { ( 1 / n ) c _ { l } J _ { l } ^ { \prime } ( n q R ) = a _ { l } J _ { l } ^ { \prime } ( q R ) + b _ { l } H _ { l } ^ { \prime } ( q R ) , } \end{array}
$$

which we solve for the scattering coefficients

$$
s _ { l } \equiv \frac { b _ { l } } { a _ { l } } = \frac { J _ { l } ^ { \prime } ( n q R ) J _ { l } ( q R ) - n J _ { l } ^ { \prime } ( q R ) J _ { l } ( n q R ) } { n J _ { l } ( n q R ) H _ { l } ^ { \prime } ( q R ) - H _ { l } ( q R ) J _ { l } ^ { \prime } ( n q R ) } .
$$

We consider now an array of $N \times N$ cilinders located at positions $\{ R _ { n } \}$ . An incident plane wave traveling along the $x$ -axis can be expressed in the frame of reference of the $n _ { t h }$ cylinder as

$$
\varPhi _ { i n } ( \pmb { r } ) = \exp ( i k X _ { n } ) \sum _ { l } i ^ { l } J _ { l } ( q r _ { n } ) \exp ( i l \theta _ { n } ) ,
$$

where $\pmb { R _ { n } } = ( X _ { n } , Y _ { n } , Z _ { n } )$ and $\mathbf { r } - \mathbf { R } _ { n }$ is described by the polar coordinates $r _ { n }$ and $\theta _ { n }$ . The Graf's addition theorem allows us to rewrite a cylindrical function centered at $\scriptstyle { R _ { n ^ { \prime } } }$ in a frame of reference centered at $\scriptstyle { \pmb R } _ { n }$ . The wave scattered by cylinder $n ^ { \prime }$ can be rewritten in the frame of reference of cylinder $n$ as

$$
\begin{array} { r l r } {  { H _ { l ^ { \prime } } ( q r _ { n ^ { \prime } } ) \exp ( i l ^ { \prime } \theta _ { n ^ { \prime } } ) = \sum _ { l } \exp ( i ( l - l ^ { \prime } ) \phi _ { n n ^ { \prime } } ) \times } } \\ & { } & { H _ { l - l ^ { \prime } } ( q R _ { n n ^ { \prime } } ) J _ { l } ( q r _ { n } ) \exp ( i l \theta _ { n } ) } \end{array}
$$

where $\pmb { R } _ { n n ^ { \prime } } = \pmb { R } _ { n } - \pmb { R } _ { n ^ { \prime } }$ is described by the polar coordinates ${ \boldsymbol { R } _ { n n ^ { \prime } } }$ and $\phi _ { n n ^ { \prime } }$ . Thus, the magnetic field in the interstices may be described in coordinates centered at the $n _ { t h }$ cylinder as the sum of the incident field, the $n$ -th scattered field and the field scattered by all the other cylinders,

$$
\begin{array} { l } { { \displaystyle \phi ^ { O } ( r ) = \exp ( i k X _ { n } ) \sum _ { l } i ^ { l } J _ { l } ( q r _ { n } ) \exp ( i l \theta _ { n } ) } } \\ { { ~ + \sum _ { l } b _ { n l } H _ { l } ( q r _ { n } ) \exp ( i l \theta _ { n } ) } } \\ { { ~ + \sum _ { n ^ { \prime } \ne n } \sum _ { l l ^ { \prime } } b _ { n ^ { \prime } l ^ { \prime } } \exp ( i ( l - l ^ { \prime } ) \phi _ { n n ^ { \prime } } ) } } \\ { { ~ \times H _ { l - l ^ { \prime } } ( q R _ { n n ^ { \prime } } ) J _ { l } ( q r _ { n } ) \exp ( i l \theta _ { n } ) } } \end{array}
$$

From Eqs. (22) and (29) we identify the coefficients $a _ { n l }$ as

$$
\begin{array} { r } { a _ { n l } = \exp ( i k X _ { n } ) i ^ { l } + \displaystyle \sum _ { n ^ { \prime } \neq n } \sum _ { l ^ { \prime } } \exp ( i ( l - l ^ { \prime } ) \phi _ { n n ^ { \prime } } ) } \\ { \times H _ { l - l ^ { \prime } } ( q R _ { n n ^ { \prime } } ) b _ { n ^ { \prime } l ^ { \prime } } . } \end{array}
$$

Introducing the scattering coefficients $s _ { n l }$ of each cylinder as in (26) yields

$$
\begin{array} { l } { { b _ { n l } - s _ { n l } \displaystyle \sum _ { n ^ { \prime } \neq n } \sum _ { l ^ { \prime } } \exp ( i ( l - l ^ { \prime } ) \phi _ { n n ^ { \prime } } ) } } \\ { { \mathrm { } \times H _ { l - l ^ { \prime } } ( q R _ { n n ^ { \prime } } ) b _ { n ^ { \prime } l ^ { \prime } } = s _ { n l } \exp ( i k X _ { n } ) i ^ { l } , } } \end{array}
$$

which can be summarized into a system of coupled equations

$$
T b = a ,
$$

where $\pmb { a } = \{ s _ { n l } \exp ( i k X _ { n } ) i ^ { l } \}$ describes the incident field, $\mathbf { \delta } _ { b } = \{ b _ { n l } \}$ describes the scattered field in the interstitial region, to be obtained, and

$$
\begin{array} { c } { { T _ { n n ^ { \prime } } ^ { l l ^ { \prime } } = \delta _ { n n ^ { \prime } } \delta _ { l l ^ { \prime } } - ( 1 - \delta _ { n n ^ { \prime } } ) \times } } \\ { { \exp ( i ( l - l ^ { \prime } ) \phi _ { n n ^ { \prime } } ) H _ { l - l ^ { \prime } } ( q R _ { n n ^ { \prime } } ) s _ { n l } . } } \end{array}
$$

3 Results We consider a metamaterial made of a square lattice of identical infinitely long dielectric cylinders of radius $R$ and refractive index $n$ set in vacuum, with a lattice constant $a$ . Using an efficient computational implementation [24] of the numerical approach presented in Sec. 2.1, we calculate the macroscopic dielectric function of the metamaterial as a function of frecuency $\omega$ and wavevector $k$ . We compare our results with the analytical solution obtained in Sec. 2.2 for a finite array of $N \times N$ cylinders.

We first consider thin weakly interacting cylinders with radius $R = 0 . 1 a$ and refractive index $n = 1 0$ Numerical calculations were performed using a 2D $6 0 1 \times 6 0 1$ grid to discretize the unit cell and performed the recursive calculation using 450 Haydock coefficients. For the analytical case we considered a large finite array of $5 1 1 \times 5 1 1$ cylinders and a maximum value of the orbital number $l = 1$ . We checked the convergence of the results by repeating the calculations with larger grids, more Haydock coefficients and larger angular momenta. Fig. 1a shows the results for the of thv et cmterial $\epsilon _ { M } ^ { T } = \epsilon _ { M } ^ { y y }$ obtaineicdietric or sive numerical approach. The results are shown as a function of the frequency, characterized by the free wavevector within the dielectric normalized to the radius nqR. Several resonances are clearly visible. In order to analyze their origin, in Fig. 1b we show the sum of the imaginary part of the scattered field coefficients $\sum _ { l } b _ { 0 l } ^ { \prime \prime }$ obtained from the analytical method. In this calculation we assumed the response of all cylinders was identical, except for the phase factor $\exp ( i k X _ { n } )$ and we included a damping of the cylindercylinder interaction at large distances to eliminate the oscillations due to reflections at the edge of the finite array, and thus mimic an infinite array. We checked convergence of this procedure by increasing the number of cylinders.

![](images/c52267b892931eff2edb6fdd615e157133e11986ee43310c76d92f490c9b62fa.jpg)  
$\epsilon _ { M } ^ { T } ~ = ~ \epsilon _ { M } ^ { y y } ( \omega , \pmb { k } )$ of the macroscopic dielectric response of a metamaterial composed of a square lattice of dielectric cylinders of refraction index $n = 1 0$ and radius $R = 0 . 1 a$ as a function of the frequency, characterized by $n q R$ , for a wavevector along the $x$ direction slightly larger than the vacuum wavevector $k ~ = ~ 1 . 0 1 q$ . (b) Sum of the imaginary part of coefficients of the scattered field $\sum _ { l } b _ { 0 l } ^ { \prime \prime }$ of an array of $5 1 1 \times 5 1 1$ dielectic cylinders, obtained from Eq. (32) using a maximum value of $l = 1$ .

Three prominent resonant features are observed in both panels of Fig. 1 at low energies. The lowest energy resonance corresponds to a magnetic dipole arising from the term $l ~ = ~ 0$ , as it lies close to that of an isolated cylinder ocurring at the first zero of the Bessel function around $n q R \ \approx \ 2 . 4$ . The second resonance around $n q R \approx \pi$ , emerges from the fulfillment of Bragg's diffraction condition for $2 a = 2 \pi / q$ . A third resonance close to $n q R \approx 3 . 8$ is originated by the term $l = 1$ .This resonance is strongly enhanced through the interaction between cylinders. The additional peaks in the macroscopic response are due to resonances caused by multiple reflections in the intesrstitial regions. We have verified that they are not due to numerical noise, but they dissapear when a very small artificial dissipation is added to the interstitial dielectric function, while the large peaks are robust. Notice however that the scattering coefficients $b _ { 0 l }$ resonate at higher energies than the macroscopic dielectric function. The reason for this discrepancy is that the transverse normal modes of the system are not actually given by the poles of the dielectric response, but by the poles of the electromagnetic Green's function.

Figure 2 a) shows the imaginary part of the electromagnetic Green's function $( \epsilon _ { T } - \check { k } ^ { 2 } / \check { q } ^ { 2 } ) ^ { - 1 }$ for a small broadening parameter $\eta = 0 . 0 0 1$ and b) the squared magnitude of the scattering amplitude $b _ { 0 } = \textstyle \sum _ { l } b _ { 0 l }$ for $k \approx q$ as a function of $n q R$ . Here, the peaks of the Green's function coincide with the peaks of the scattering coefficient as the poles of the Green's function and of the scattered coefficient correspond both to the normal modes of the system, for which one may have a finite field, and finite scattered amplitudes without an external excitation.

![](images/4ec90b9b7cdc57b107de1bcc6b2a485855592ee6753d08ed8c01eca1210a9798.jpg)  
Figure 2 a) Normalized imaginary part of the electromagnetic Green's function $\eta ^ { 2 } / ( \eta ^ { \overleftarrow { 2 } } + ( \overleftarrow { \epsilon _ { T } } - k ^ { 2 } / q ^ { 2 } ) ^ { 2 } )$ for a small dissipation parameter $\eta$ as a function of $n q R$ for $k \approx q$ for the same system as in Fig. 1. b) Squared magnitude of the scattering amplitude $b _ { 0 } = \sum _ { l } b _ { 0 l }$ for the same system as in Fig. 1.

Having verified the consistency of both computational approaches when applied to the calculation of the normal modes of a system, we now consider a more interesting and realistic case. We consider a metamaterial composed of strongly interacting cylinders with a larger radius $R =$ $0 . 3 5 a$ and a large but realistic refraction index $n = 4$ (similar to that of Si). In order to identify exotic behavior such as negative refraction, we have to examine the dispersion relation of the normal modes of the system and explore their group velocity [20,25]. Thus, we calculated the macroscopic dielectric function and obtained the Green's function of the metamaterial as a function of both frecuency $\omega$ and wavevector $k$ . For the numerical calculations we used a two dimensional $2 0 1 \times 2 0 1$ grid and 300 sets of Haydock coefficients. A very small artificial dissipative term $0 . 0 0 1 i$ has been added to the vaccum dielectric constant in order to improve convergence. The results are displayed in Fig. 3. The scattering amplitude obtained from the scattering approach, calculated for an array of $2 0 1 \times 2 0 1$ cylinders considering a maximum value of the orbital number $l = 2$ is shown for comparisson.

The upper panel of Fig. 3 shows the imaginary part of the of the Green's function $\eta / ( \eta ^ { 2 } + ( \epsilon _ { M } ^ { T } - k ^ { \overline { { 2 } } } / q ^ { 2 } ) )$ where $\epsilon _ { M } ^ { T }$ factor $\eta \ : = \ : 0 . 2$ for better visualization. The lower panel shows the absolute value squared of the scattering amplitude $b _ { 0 } = \textstyle \sum _ { l } b _ { 0 l }$ , smoothed as $\eta / ( \eta ^ { 2 } + 1 / | b _ { 0 } | ^ { 2 } )$ with a dissipation factor $\eta = 0 . 1$ .

![](images/2f9746e613d3956766104e24f077c0d2d605c9dc3873c6d7b76f43b71c61af41.jpg)  
Figure 3 Imaginary part of the Green's function $( \epsilon _ { T } -$ $k ^ { 2 } \bar { / } q ^ { 2 } ) ^ { - 1 }$ for a system composed by dielectic cylinders of radius $r \ : = \ : 0 . 3 5 a$ and refractive index $n = 4$ , obtained numerically through the macroscopic response using the package Photonic (upper panel). Magnitude of the scattered amplitude $b _ { 0 }$ normalized as $\eta / \big ( \eta ^ { 2 } + 1 / b _ { 0 } ^ { 2 } \big )$ using a dissipation factor $\eta$ for an array of $2 0 1 \times 2 0 1$ cylinders, obtained through the scattering approach (Eq. (32)) using a maximum value of $l = 2$ (lower panel). Color map is logarithmic.

Several bands may be identified in Fig. 3 showing a very good agreement between the dispersion relations obtained through our two approaches. Regions of negative dispersion, i.e., for which the frequency of the resonances decreases as the wavevector increases, yielding a negative group velocity, are clearly observed for the second, fourth and fifth bands and for large wavevectors before the first Brillouin zone boundary at $k a = \pi$ .The bands originate from the combination of isolated cylinder resonances of different values of $l$ ,while different values of $l$ participate in each band due to the relatively strong interaction between neighboring thick cylinders. This can be confirmed through Figure 4, which shows the absolute value of the individual scattering coefficients $b _ { 0 l }$ for $l = 0 , 1 , 2$ as a function of $n q R$ for a) a single cylinder, b) an array of $1 1 \times 1 1$ cylinders and c) an array of $2 0 1 \times 2 0 1$ , for a wavevector $k = q$ . The coefficients corresponding to an isolated cylinder display one broad peak each around $q = 2 . 4$ , 4.8 and 5.1 for $l = 0 , 1 , 2$ respectively, which are close to the first zeroes of the corresponding Bessel function. For an array of $1 1 \times 1 1$ cylinders these peaks appear distorted and shifted to be finally merged for a larger array of $2 0 1 \times 2 0 1$ in which case all $l$ -contributions resonate at the same frequencies.

4 Conclusions We presented two schemes to calculate the electromagnetic properties of a metamaterial made of a simple lattice of cylinders: a numerical method based on a recursive calculation of the macroscopic dielectric tensor which may be easily generalizable to arbitrary geometries and materials, and a multiple scattering approach for cylindrical geometries which allowed us a simple interpretation of the results in terms of the excitation of Mie resonances. We applied these methods to investigate the response of a system made up of cylinders of high index of refraction. The comparison between the results of both methods is not direct, as in one case we obtain a macroscopic response and in the other we obtain scattering amplitudes. Nevertheless, we showed that the poles of the macroscopic Green's function obtained by our numerical method coincide with those of the scattering coefficients, and both can be interpreted as the excitation of the normal modes of the system. For a system of thin cylinders with very high index of refraction and a relatively small coupling, we identified the nature of each mode. We found modes arising from the Mie resonances of individual cylinders and a mode arising from Bragg coherent multiple scattering. For larger cylinders the interaction yields a coupling between resonances with different angular momenta. By varying the frequency and wavevector independently we computed the dispersion relation of the normal modes. The photonic band structure obtained using both methods is in very good agreement, and reveal regions of negative dispersion. Thus, through comparison with an ad-hoc model we showed that our macroscopic approach based on Haydock's recursion and its implementation in the Photonic package is an efficient procedure for obtaining the optical properties of metamaterials made up of high index of refraction cylinders incorporating resonances which cannot be explored within the long wavelength limit. Furthermore, as it can be readily generalized to arbitrary geometries and materials, we believe it will prove to be a useful tool for the design of artificial materials with a richer geometry that might yield novel properties.

![](images/19f03adc7cb7c3ec7e8e0b06445b576a34245e4eb39f08447480140afe133655.jpg)  
Figure 4 Absolute value of the scattering coefficients $| b _ { 0 l } |$ for an a) isolated cylinder, b) array of $1 1 \times 1 1$ cylinders and array of $2 0 1 \times 2 0 1$ cylinders, as a function of $n q R$ .

Acknowledgements This work has been supported by CONACyT through a postdoctoral research fellowship. WLM acknowledges the support from DGAPA-UNAM through grant IN111119. BSM acknowledges the support from CONACYT through grant A1-S-9410.

# References

[1] R.A. Shelby, D.R. Smith, and S. Schultz, Science 292(5514), 7779 (2001).   
[2] D.R. Smith, J.B. Pendry, and M.C. Wiltshire, Science 305(5685), 788792 (2004).   
[3] K. Aydin, I. Bulu, K. Guven, M. Kafesaki, C.M. Soukoulis, and E. Ozbay, New Journal of Physics 7(1), 168 (2005).   
[4] A.J. Hoffman, L. Alekseyev, S.S. Howard, K. J. Franz, D. Wasserman, V. A. Podolskiy, E. E. Narimanov, D.L. Sivco, and C. Gmachl, Nature materials 6(12), 946 (2007).   
[5] L. Peng, L. Ran, H. Chen, H. Zhang, J. A. Kong, and T. M. Grzegorczyk, Physical Review Letters 98(15), 157403 (2007).   
[6] A. I. Kuznetsov, A. E. Miroshnichenko, M. L. Brongersma, Y.S. Kivshar, and B. Lukyanchuk, Science 354(6314), aag2472 (2016).   
[7] S. Linden, C. Enkrich, M. Wegener, J. Zhou, T. Koschny, and C.M. Soukoulis, Science 306(5700), 13511353 (2004).   
[8] T.J. Yen, W.J. Padilla, N. Fang, D.C. Vier, D.R. Smith, J. B. Pendry, D.N. Basov, and X. Zhang, Science 303(5663), 14941496 (2004).   
[9] Q. Zhao, J. Zhou, F. Zhang, and D. Lippens, Materials Today 12(12), 60 69 (2009).   
[10] Y. Kivshar and A. Miroshnichenko, Opt. Photon. News 28(1), 2431 (2017).   
[11] S. Kruk and Y. Kivshar, ACS Photonics 4(11), 26382649 (2017).   
[12] K. Vynck, D. Felbacq, E. Centeno, A. I. Cbuz, D. Cassagne, and B. Guizal, Physical Review Letters 102(Mar), 133901 (2009).   
[1] J.A. Schuller, R. Zia, T. Taubner, and M.L. Brongersma, Physical Review Letters 99(10), 107401 (2007).   
[14] S. Jahani and Z. Jacob, Nature nanotechnology 11(1), 23 (2016).   
[15] C. Zhang, Y. Xu, J. Liu, J. Li, J. Xiang, H. Li, J. Li, Q. Dai, S. Lan, and A. Miroshnichenko, Lighting up silicon nanoparticles with mie resonances, 2018.   
[1] D. Felbacq, G. Tayeb, and D. Maystre, Journal of the Optical Society of America A 11(9), 25262538 (1994).   
[17] W. L. Mochán, G. P. Ortiz, and B. S. Mendoza, Optics express 18(21), 2211922127 (2010).   
[18] U. R. Meza, B.S. Mendoza, and W. L. Mochán, Physical Review B 99(12), 125408 (2019).   
[19] J. S. Pérez-Huerta, G. P. Ortiz, B. S. Mendoza, and W.L. Mochán, New Journal of Physics 15(4), 043037 (2013).   
[20] V. M. Agranovich and Y. N. Gartstein, Physics-Uspekhi 49(10), 1029 (2006).   
[21] R. Haydock, Computer Physics Communications 20(1), 1116 (1980).   
[22] L. Juárez-Reyes and W. L. Mochán, Physica Status Solidi (b) 255(4), 1700495 (2018).   
[23] D. Gagnon and L. J. Dubé, Journal of Optics 17(10), 103501 (2015).   
[24] W.L. Mochán, G. Ortiz, B. S. Mendoza, and J. S. PérezHuerta, Photonic, Comprehensive Perl Archive Network (CPAN), 2016, Perl package for calculations on metamaterials and photonic structures.   
[25] V. Agranovich and Y. N. Gartstein, Metamaterials 3(1), 1− 9 (2009).