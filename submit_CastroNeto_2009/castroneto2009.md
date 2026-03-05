# The electronic properties of graphene

A. H. Castro Neto1, F. Guinea2, N. M. R. Peres3, K. S. Novoselov4, and A. K. Geim4

1Department of Physics, Boston University, 590 Commonwealth Avenue, Boston, MA 02215,USA Instituto de Ciencia de Materiales de Madrid. CSIC. Cantoblanco. E-28049 Madrid, Spain Center of Physics and Department of Physics, Universidade do Minho, P-4710-057, Braga, Portugal and 4Department of Physics and Astronomy, University of Manchester, Manchester, M13 9PL, UK

(Dated: November 26, 2024)

This article reviews the basic theoretical aspects of graphene, a one atom thick allotrope of carbon, with unusual two-dimensional Dirac-like electronic excitations. The Dirac electrons can be controlled by application of external electric and magnetic fields, or by altering sample geometry and/or topology. We show that the Dirac electrons behave in unusual ways in tunneling, confinement, and integer quantum Hall effect. We discuss the electronic properties of graphene stacks and show that they vary with stacking order and number of layers. Edge (surface) states in graphene are strongly dependent on the edge termination (zigzag or armchair) and affect the physical properties of nanoribbons. We also discuss how different types of disorder modify the Dirac equation leading to unusual spectroscopic and transport properties. The effects of electron-electron and electron-phonon interactions in single layer and multilayer graphene are also presented.

# Contents

I. Introduction

II. Elementary electronic properties of graphene

A. Single layer: tight-binding approach 1. Cyclotron mass 2. Density of states   
B. Dirac fermions 1. Chiral Tunneling and Klein paradox 2. Confinement and zitterbewegung   
C. Bilayer graphene: tight-binding approach   
D. Epitaxial graphene   
E. Graphene stacks 1. Electronic structure of bulk graphite   
F. Surface states in graphene   
G. Surface states in graphene stacks   
H. The spectrum of graphene nanoribbons 1. Zigzag nanoribbons 2. Armchair nanoribbons I. Dirac fermions in a magnetic field   
J. The anomalous integer quantum Hall effec   
K. Tight-binding model in a magnetic field   
L. Landau levels in graphene stacks   
M. Diamagnetism   
N. Spin orbit coupling

III. Flexural phonons, elasticity, and crumpling

IV. Disorder in graphene

A. Ripples   
B. Topological lattice defects   
C. Impurity states   
D. Localized states near edges, cracks, and voids   
E. Self-doping   
F. Vector potential and gauge field disorder 1. Gauge field induced by curvature 2. Elastic strain 3. Random gauge fields   
G. Coupling to magnetic impurities   
H. Weak and strong localization   
I. Transport near the Dirac point   
J. Boltzmann Equation description of DC transport in doped graphene

26   
28   
29   
30   
30   
31   
32   
32   
33   
34   
34   
34   
36   
37

K. Magnetotransport and universal conductivity 38 1. The full self-consistent Born approximation (FSBA) 38

V. Many-body effects 40

A. Electron-phonon interactions 40   
B. Electron-electron interactions 42   
1. Screening in graphene stacks 45   
C. Short range interactions 45   
1. Bilayer graphene: exchange 46   
2. Bilayer graphene: short range interactions 46   
D. Interactions in high magnetic fields 47

VI. Conclusions 47

VII. Acknowledgments 48   
References 48

# I. INTRODUCTION

Carbon is the materia prima for life on the planet and the basis of all organic chemistry. Because of the flexibility of its bonding, carbon-based systems show an unlimited number of different structures with an equally large variety of physical properties. These physical properties are, in great part, the result of the dimensionality of these structures. Among systems with only carbon atoms, graphene - a two-dimensional (2D) allotrope of carbon - plays an important role since it is the basis for the understanding of the electronic properties in other allotropes. Graphene is made out of carbon atoms arranged on a honeycomb structure made out of hexagons (see Fig. 1), and can be thought as composed of benzene rings stripped out from their hydrogen atoms (Pauling, 1972). Fullerenes (Andreoni, 2000) are molecules where carbon atoms are arranged spherically, and hence, from the physical point of view, are zero-dimensional objects with discrete energy states. Fullerenes can be obtained from graphene with the introduction of pentagons (that create positive curvature defects), and hence, fullerenes can be thought as wrapped up graphene. Carbon nanotubes (Charlier et al., 2007; Saito et al., 1998) are obtained by rolling graphene along a given direction and reconnecting the carbon bonds. Hence, carbon nanotubes have only hexagons and can be thought as onedimensional (1D) objects. Graphite, a three dimensional (3D) allotrope of carbon, became widely known to mankind after the invention of the pencil in 1564 (Petroski, 1989) and its usefulness as an instrument for writing comes from the fact that graphite is made out of stacks of graphene layers that are weakly coupled by van der Waals forces. Hence, when one presses a pencil against a sheet of paper one is actually producing graphene stacks and, somewhere among them, there could be individual graphene layers. Although graphene is the mother for all these different allotropes and has been presumably produced every time someone writes with a pencil, it was only isolated 440 years after its invention (Novoselov et al., 2004). The reason is that, first, no one actually expected graphene to exist in the free state and, second, even with the benefit of hindsight, no experimental tools existed to search for one-atom-thickfakes among the pencil debris covering macroscopic areas (Geim and MacDonald, 2007). Graphene was eventually spotted due to the subtle optical effect it creates on top of a cleverly chosen SiO $^ 2$ substrate (Novoselov et al., 2004) that allows its observation with an ordinary optical microscope (Abergel et al., 2007; Blake et al., 2007; Casiraghi et al., 2007). Hence, graphene is relatively straightforward to make, but not so easy to find.

The structural flexibility of graphene is reflected in its electronic properties. The sp $^ 2$ hybridization between one s-orbital and two p-orbitals leads to a trigonal planar structure with a formation of a $\sigma$ -bond between carbon atoms which are separated by 1.42 Å. The $\sigma$ -band is responsible for the robustness of the lattice structure in all allotropes. Due to the Pauli principle these bands have a filled shell and hence, form a deep valence band. The unaffected p-orbital, which is perpendicular to the planar structure, can bind covalently with neighboring carbon atoms leading to the formation of a $\pi$ -band. Since each p-orbital has one extra electron, the $\boldsymbol { \mathscr { u } }$ -band is half-filled.

Half-filled bands in transition elements have played an important role in the physics of strongly correlated systems since, due to its strong tight binding character, the Coulomb energies are very large, leading to strong collective effects, magnetism, and insulating behavior due to correlation gaps or Mottness (Phillips, 2006). In fact, Linus Pauling proposed in the 1950's that, on the basis of the electronic properties of benzene, graphene should be a resonant valence bond structure (RVB) (Pauling, 1972). RVB states have become very popular in the literature of transition metal oxides, and particularly in studies of cuprate oxides superconductors (Maple, 1998). This point of view, should be contrasted with contemporaneous band structure studies of graphene (Wallace,

![](images/f0af188793762f2ac7e8e97a42a843226845e9ab99e0879d5866f8825b27e498.jpg)  
Figure 1 (Color online) Graphene (top left) is a honeycomb lattice of carbon atoms. Graphite (top right) can be viewed a stack of graphene layers. Carbon nanotubes are rolledup cylinders of graphene (bottom left). Fullerenes (C60) are molecules consisting of wrapped graphene by the introduction of pentagons on the hexagonal lattice (Castro Neto et al., 2006a).

1947) that found it to be a semimetal with unusual linearly dispersing electronic excitations called Dirac electrons. While most of the current experimental data in graphene supports the band structure point of view, the role of the electron-electron interactions in graphene is a subject of intense research.

It was P. R. Wallace who in 1946 wrote the first papers on the band structure of graphene and showed the unusual semimetallic behavior in this material (Wallace, 1947). At that point in time, the thought of a purely 2D structure was a mere fantasy and Wallace's studies of graphene served him as a starting point to study graphite, a very important material for nuclear reactors in the post-World War II era. During the following years, the study of graphite culminated with the Slonczewski-Weiss-McClure (SWM) band structure of graphite which provided a detailed description of the electronic properties in this material (McClure, 1957; Slonczewski and Weiss, 1958) and was very successful in describing the experimental data (Boyle and Nozières, 1958; Dillon et al., 1977; McClure, 1958, 1964; Soule et al., 1964; Spry and Scherer, 1960; Williamson et al., 1965). Interestingly enough, from 1957 to 1968, the assignment of the electron and hole states within the SWM model were the opposite to what is accepted today. In 1968, Schroeder et al. (Schroeder et al., 1968) established the currently accepted location of electron and hole pockets (McClure, 1971). The SWM model has been revisited in recent years because of its inability to describe the van der Waals-like interactions between graphene planes, a problem that requires the understanding of many-body effects that go beyond the band structure description (Rydberg et al., 2003). These issues, however, do not arise in the context of a single graphene crystal but they show up with great importance when graphene layers are stacked on top of each other, as in the case, for instance, of the bilayer graphene. Stacking can change the electronic properties considerably and the layering structure can be used in order to control the electronic properties.

One of the most interesting aspects of the graphene problem is that its low energy excitations are massless, chiral, Dirac fermions. In neutral graphene the chemical potential crosses exactly the Dirac point. This particular dispersion, that is only valid at low energies, mimics the physics of quantum electrodynamics (QED) for massless fermions except by the fact that in graphene the Dirac fermions move with a speed $v _ { F }$ which is 300 times smaller than the speed of light, c. Hence, many of the unusual properties of QED can show up in graphene but at much smaller speeds (Castro Neto et al., 2006a; Katsnelson and Novoselov, 2007; Katsnelson et al., 2006). Dirac fermions behave in very unusual ways when compared to ordinary electrons if subjected to magnetic fields, leading to new physical phenomena (Gusynin and Sharapov, 2005; Peres et al., 2006c) such as the anomalous integer quantum Hall effect (IQHE) measured experimentally (Novoselov et al., 2005a; Zhang et al., 2005). Besides being qualitatively different from the IQHE observed in Si and GaAlAs (heterostructures) devices (Stone, 1992), the IQHE in graphene can be observed at room temperature because of the large cyclotron energies for "relativistic" electrons (Novoselov et al., 2007). In fact, the anomalous IQHE is the trademark of Dirac fermion behavior.

Another particularly interesting feature of Dirac fermions is their insensitivity to external electrostatic potentials due to the so-called Klein paradox, that is, the fact that Dirac fermions can be transmitted with probability one through a classically forbidden region (Calogeracos and Dombey, 1999; Itzykson and Zuber, 2006). In fact, Dirac fermions behave in a very unusual way in the presence of confining potentials leading to the phenomenon of zitterbewegung, or jittery motion of the wavefunction (Itzykson and Zuber, 2006). In graphene these electrostatic potentials can be easily generated by disorder. Since disorder is unavoidable in any material, there has been great interest in trying to understand how disorder affects the physics of electrons in graphene and its transport properties. In fact, under certain conditions, Dirac fermions are immune to localization effects observed in ordinary electrons (Lee and Ramakrishnan, 1985) and it has been established experimentally that electrons can propagate without scattering over large distances of the order of micrometers in graphene (Novoselov et al., 2004). The sources of disorder in graphene are many and can vary from ordinary effects commonly found in semiconductors, such as ionized impurities in the Si substrate, to adatoms and various molecules adsorbed in the graphene surface, to more unusual defects such as ripples associated with the soft structure of graphene (Meyer et al., 2007a). In fact, graphene is unique in the sense that it shares properties of soft membranes (Nelson et al., 2004) and at the same time it behaves in a metallic way, so that the Dirac fermions propagate on a locally curved space. Here, analogies with problems of quantum gravity become apparent (Fauser et al., 2007). The softness of graphene is related with the fact that it has out-ofplane vibrational modes (phonons) that cannot be found in 3D solids. These flexural modes, responsible for the bending properties of graphene, also account for the lack of long range structural order in soft membranes leading the phenomenon of crumpling (Nelson et al., 2004). Nevertheless, the presence of a substrate or scaffolds that hold graphene in place, can stabilize a certain degree of order in graphene but leaves behind the so-called ripples (which can be viewed as frozen flexural modes).

It was realized very early on that graphene should also present unusual mesoscopic effects (Katsnelson, 2007a; Peres et al., 2006a). These effects have their origin in the boundary conditions required for the wavefunctions in mesoscopic samples with various types of edges graphene can have (Akhmerov and Beenakker, 2007; Nakada et al., 1996; Peres et al., 2006c; Wakabayashi et al., 1999). The most studied edges, zigzag and armchair, have drastically different electronic properties. Zigzag edges can sustain edge (surface) states and resonances that are not present in the armchair case. Moreover, when coupled to conducting leads, the boundary conditions for a graphene ribbon strongly affects its conductance and the chiral Dirac nature of the fermions in graphene can be exploited for applications where one can control the valley flavor of the electrons besides its charge, the so-called valleytronics (Rycerz et al., 2007). Furthermore, when superconducting contacts are attached to graphene, they lead to the development of supercurrent flow and Andreev processes characteristic of superconducting proximity effect (Heersche et al., 2007). The fact that Cooper pairs can propagate so well in graphene attests for the robust electronic coherence in this material. In fact, quantum interference phenomena such as weak localization, universal conductance fluctuations (Morozov et al., 2006), and the Aharonov-Bohm effect in graphene rings have already been observed experimentally (Recher et al., 2007; Russo et al., 2007). The ballistic electronic propagation in graphene can be used for field effect devices such as p-n (Cheianov et al., 2007a; Cheianov and Fal'ko, 2006; Fogler et al., 2007a; Huard et al., 2007; Lemme et al., 2007; Tworzydlo et al., 2007; Williams et al., 2007; Zhang and Fogler, 2007) and p-n-p (Ossipov et al., 2007) junctions, and as "neutrino" billiards (Berry and Modragon, 1987; Miao et al., 2007). It has also been suggested that Coulomb interactions are considerably enhanced in smaller geometries, such as graphene quantum dots (Milton Pereira Junior et al., 2007), leading to unusual Coulomb blockade effects (Geim and Novoselov, 2007) and perhaps to magnetic phenomena such as the Kondo effect. The amazing transport properties of graphene allow for their use in a plethora of applications ranging from single molecule detection (Schedin et al., 2007; Wehling et al., 2007) to spin injection (Cho et al., 2007; Hill et al., 2007; Ohishi et al., 2007; Tombros et al., 2007).

Because of its unusual structural and electronic flexibility, graphene can be tailored chemically and/or structurally in many different ways: deposition of metal atoms (Calandra and Mauri, 2007; Uchoa et al., 2007) or molecules (Leenaerts et al., 2007; Schedin et al., 2007; Wehling et al., 2007) on top; intercalation (as it is done in graphite intercalated compounds (Dresselhaus and Dresselhaus, 2002; Dresselhaus et al., 1983; Tanuma and Kamimura, 1985)); incorporation of nitrogen and/and boron in its structure (Martins et al., 2007; Peres et al., 2007a) (in analogy with what has been done in nanotubes (Stephan et al., 1994)); using different substrates that modify the electronic structure (Calizo et al., 2007; Das et al., 2007; Faugeras et al., 2007; Giovannetti et al., 2007; Varchon et al., 2007; Zhou et al., 2007). The control of graphene properties can be extended in new directions allowing for creation of graphene-based systems with magnetic and superconducting properties (Uchoa and Castro Neto, 2007) that are unique in their 2D properties. Although the graphene field is still in its infancy, the scientific and technological possibilities of this new material seem to be unlimited. The understanding and control of the properties of this material can open doors for a new frontier in electronics. As the current status of the experiment and potential applications have recently been reviewed (Geim and Novoselov, 2007), in this article we mostly concentrate on the theory and more technical aspects of electronic properties of this exciting new material.

# II. ELEMENTARY ELECTRONIC PROPERTIES OFGRAPHENE

# A. Single layer: tight-binding approach

Graphene is made out of carbon atoms arranged in hexagonal structure as shown in Fig. 2. The structure is not a Bravais lattice but can be seen as a triangular lattice with a basis of two atoms per unit cell. The lattice vectors can be written as:

$$
a _ { 1 } = { \frac { a } { 2 } } ( 3 , \sqrt { 3 } ) , \qquad a _ { 2 } = { \frac { a } { 2 } } ( 3 , - \sqrt { 3 } ) ,
$$

where $a \approx 1 . 4 2$ Å is the carbon-carbon distance. The reciprocal lattice vectors are given by:

Brillouin zone (BZ). These are named Dirac points for reasons that will become clear later. Their positions in momentum space are given by:

$$
K = \left( { \frac { 2 \pi } { 3 a } } , { \frac { 2 \pi } { 3 { \sqrt { 3 } } a } } \right) , \qquad K ^ { \prime } = \left( { \frac { 2 \pi } { 3 a } } , - { \frac { 2 \pi } { 3 { \sqrt { 3 } } a } } \right) .
$$

The three nearest neighbors vectors in real space are given by:

$$
\delta _ { 1 } = { \frac { a } { 2 } } ( 1 , \sqrt { 3 } ) \qquad \delta _ { 2 } = { \frac { a } { 2 } } ( 1 , - \sqrt { 3 } ) \qquad \delta _ { 3 } = - a ( 1 , 0 )
$$

while the six second-nearest neighbors are located at: $\delta _ { 1 } ^ { \prime } = \pm { \pmb a } _ { 1 } , \delta _ { 2 } ^ { \prime } = \pm { \pmb a } _ { 2 } , \delta _ { 3 } ^ { \prime } = \pm ( { \pmb a } _ { 2 } - { \pmb a } _ { 1 } )$ .

Of particular importance for the physics of graphene are the two points $K$ and $K ^ { \prime }$ at the corners of the graphene

The tight-binding Hamiltonian for electrons in graphene considering that electrons can hop both to nearest and next nearest neighbor atoms has the form (we use units such that $\hbar = 1$ ):

$$
\begin{array} { l } { { \displaystyle H = \mathrm { ~ - ~ } t \sum _ { \langle i , j \rangle , \sigma } \left( a _ { \sigma , i } ^ { \dagger } b _ { \sigma , j } + \mathrm { h . c . } \right) } } \\ { { \mathrm { ~ - ~ } t ^ { \prime } \sum _ { \langle \langle i , j \rangle \rangle , \sigma } \left( a _ { \sigma , i } ^ { \dagger } a _ { \sigma , j } + b _ { \sigma , i } ^ { \dagger } b _ { \sigma , j } + \mathrm { h . c . } \right) , } } \end{array}
$$

where $u _ { i , \sigma }$ $a _ { i , \sigma } ^ { \dag } )$ annihilates (creates) an electron with spin $\sigma$ $\sigma = \uparrow , \downarrow$ )on site $\mathbf { R } _ { i }$ on sublattice A (an equivalent definition is used for sublattice B), $t$ $\approx 2 . 8 \ \mathrm { e V }$ is the nearest neighbor hopping energy (hopping between different sublattices), $t ^ { \prime }$ 1is the next nearest neighbor hopping energy (hopping in the same sublattice). The energy bands derived from this Hamiltonian have the form (Wallace, 1947):

$$
\begin{array} { l } { { \displaystyle E _ { \pm } ( { \bf k } ) = \pm t \sqrt { 3 + f ( { \bf k } ) } - t ^ { \prime } f ( { \bf k } ) } , \ ~ } \\ { { \displaystyle f ( { \bf k } ) = 2 \cos \left( \sqrt { 3 } k _ { y } a \right) + 4 \cos \left( \frac { \sqrt { 3 } } { 2 } k _ { y } a \right) \cos \left( \frac { 3 } { 2 } k _ { x } a \right) } , \ ~ } \end{array}
$$

where the plus sign applies to the upper $( \pi )$ and the minus sign the lower $( \pi ^ { * } )$ band. It is clear from (6) that the spectrum is symmetric around zero energy if $t ^ { \prime } =$ 0. For finite values of $t ^ { \prime }$ the electron-hole symmetry is broken and the $\pi$ and $\pi ^ { * }$ bands become asymmetric. In Fig. 3 we show the full band structure of graphene with both $t$ and $t ^ { \prime }$ . In the same figure we also show a zoom in of the band structure close to one of the Dirac points (at the $\mathrm { K }$ or $\mathrm { K } '$ point in the BZ). This dispersion can be obtained by expanding the full band structure, eq.(6),

$$
b _ { 1 } = \frac { 2 \pi } { 3 a } ( 1 , \sqrt { 3 } ) , ~ b _ { 2 } = \frac { 2 \pi } { 3 a } ( 1 , - \sqrt { 3 } ) .
$$

![](images/c8b092694b085f527aa640596d92e87f477f891440650753b9e0a7bd6149eb6c.jpg)  
Figure 2 (Color online) Left: Lattice structure of graphene, made out of two interpenetrating triangular lattices ( $\mathbf { \delta } _ { \mathbf { u } 1 }$ and $\mathbf { \delta } _ { a 2 }$ are the lattice unit vectors, and $\delta _ { i }$ , $i = 1 , 2 , 3$ are the nearest neighbor vectors); Right: corresponding Brillouin zone. The Dirac cones are located at the K and $\mathrm { K } '$ points.

close to the $\mathbf { K }$ (or $\mathbf { K } ^ { \prime }$ ) vector, eq.(3), as: $\mathbf { k } = \mathbf { K } + \mathbf { q }$ , with $| \mathbf { q } | \ll | \mathbf { K } |$ (Wallace, 1947):

$$
E _ { \pm } ( { \bf q } ) \approx \pm v _ { F } | { \bf q } | + { \mathcal O } ( ( q / K ) ^ { 2 } ) ,
$$

where $\mathbf { \Delta } _ { q }$ is the momentum measured relatively to the Dirac points and $v _ { F }$ represents the Fermi velocity, given by $v _ { F } = 3 t a / 2$ , with a value $v _ { F } \simeq 1 \times 1 0 ^ { 6 } ~ \mathrm { m / s }$ . This result was first obtained by Wallace (Wallace, 1947).

The most striking difference between this result and the usual case, $\epsilon ( { \bf q } ) = q ^ { 2 } / ( 2 m )$ where $m$ is the electron mass, is that the Fermi velocity in (7) does not depend on the energy or momentum: in the usual case we have $v = k / m = \sqrt { 2 E / m }$ and hence the velocity changes substantially with energy. The expansion of the spectrum around the Dirac point including $t ^ { \prime }$ up to second order in $q / K$ is given by:

$$
E _ { \pm } ( { \bf q } ) \simeq 3 t ^ { \prime } \pm v _ { F } | { \bf q } | - \left( \frac { 9 t ^ { \prime } a ^ { 2 } } { 4 } \pm \frac { 3 t a ^ { 2 } } { 8 } \sin ( 3 \theta _ { \bf q } ) \right) | { \bf q } | ^ { 2 } ,
$$

where

$$
\theta _ { { \bf q } } = \arctan \left( { \frac { q _ { x } } { q _ { y } } } \right) ,
$$

is the angle in momentum space. Hence, the presence of $t ^ { \prime }$ shifts in energy the position of the Dirac point and breaks electron-hole symmetry. Notice that up to order $( q / K ) ^ { 2 }$ the dispersion depends on the direction in momentum space and has a three fold symmetry. This is the so-called trigonal warping of the electronic spectrum (Ando et al., 1998; Dresselhaus and Dresselhaus, 2002).

# 1. Cyclotron mass

The energy dispersion (7) resembles the energy of ultra-relativistic particles; these particles are quantum mechanically described by the massless Dirac equation (see section II.B for more on this analogy). An immediate consequence of this massless Dirac-like dispersion is a cyclotron mass that depends on the electronic density as its square root (Novoselov et al., 2005a; Zhang et al., 2005). The cyclotron mass is defined, within the semiclassical approximation (Ashcroft and Mermin, 1976), as

![](images/c5c9b370f66b65bdf809b05c9d1ce162be1a7a3ca1ea08e6c43e969dd0e4279b.jpg)  
Figure 3 (Color online) Left: Energy spectrum (in units of $t$ ) for finite values of $t$ and $t ^ { \prime }$ , with $t = 2 . 7 \ \mathrm { e V }$ and $t ^ { \prime } = 0 . 2 t$ . Right: zoom-in of the energy bands close to one of the Dirac points.

$$
m ^ { * } = \frac { 1 } { 2 \pi } \left[ \frac { \partial A ( E ) } { \partial E } \right] _ { E = E _ { F } } ,
$$

with $A ( E )$ the area in $k -$ space enclosed by the orbit and given by:

$$
A ( E ) = \pi q ( E ) ^ { 2 } = \pi { \frac { E ^ { 2 } } { v _ { F } ^ { 2 } } } .
$$

Using (11) in (10) one obtains:

$$
m ^ { * } = \frac { E _ { F } } { v _ { F } ^ { 2 } } = \frac { k _ { F } } { v _ { F } } .
$$

The electronic density, $n$ , is related to the Fermi momentum, $k _ { F }$ , as $k _ { F } ^ { 2 } / \pi = n$ (with contributions from the two Dirac points $\pmb { K }$ and $\pmb { K } ^ { \prime }$ and spin included) which leads to:

$$
m ^ { * } = \frac { \sqrt { \pi } } { v _ { F } } \sqrt { n } .
$$

Fitting (13) to the experimental data (see Fig.4) provides an estimation for the Fermi velocity and the hopping parameter as $v _ { F } ~ \approx ~ 1 0 ^ { 6 } \mathrm { m s ^ { - 1 } }$ and $t \ \approx \ 3 \mathrm { e V }$ respectively. The experimental observation of the $\sqrt { n }$ dependence of the cyclotron mass provides evidence for the existence of massless Dirac quasiparticles in graphene (Deacon et al., 2007; Jiang et al., 2007a; Novoselov et al., 2005a; Zhang et al., 2005) - the usual parabolic (Schrödinger) dispersion implies a constant cyclotron mass.

# 2Density of states

The density of states per unit cell, derived from (6), is given in Fig. 5 for both $t ^ { \prime } = 0$ and $t ^ { \prime } \neq 0$ , showing in both cases a semimetallic behavior (Bena and Kivelson, 2005; Wallace, 1947). For $t ^ { \prime } = 0$ it is possible to derive an analytical expression for the density of states per unit cell, which has the form (Hobson and Nierenberg, 1953):

![](images/9261f926874fb037aafa3e3ca007a595324755574c237f777aab44da5ec4eb4e.jpg)  
Figure 4 (Color online) Cyclotron mass of charge carriers in graphene as a function of their concentration $n$ . Positive and negative $n$ correspond to electrons and holes, respectively. Symbols are the experimental data extracted from temperature dependence of the SdH oscillations; solid curves is the best fit to Eq. (13). $m _ { 0 }$ is the free electron mass. Adapted from Novoselov et al., $2 0 0 5 \mathrm { a }$ .

where $A _ { c }$ is the unit cell area given by $A _ { c } = 3 \sqrt { 3 } a ^ { 2 } / 2$ . It is worth noting that the density of states for graphene is very different from the density of states of carbon nanotubes (Saito et al., 1992a,b). The latter show $1 / \sqrt { E }$ singularities due to the 1D nature of their electronic spectrum, which comes about due to the quantization of the momentum in the direction perpendicular to the tube axis. From this perspective, graphene nanoribbons, which also have momentum quantization perpendicular to the ribbon length, have properties very similar to carbon nanotubes.

![](images/ce22c751c399f0439dae9aaf4051183aea2026749b8512c8238521b2744c3896.jpg)  
Figure 5 (Color online) Density of states per unit cell as a function of energy (in units of $t$ ) computed from the energy dispersion (6), $t ^ { \prime } = 0 . 2 t$ (top) and for $t ^ { \prime } = 0$ (bottom). Also shown is a zoom in of the density of states close to the neutrality point of one electron per site. For the case $t ^ { \prime } = 0$ the electron-hole nature of the spectrum is apparent and the density of states close to the neutrality point can be approximated by $\rho ( \epsilon ) \propto | \epsilon |$ .

$$
\begin{array} { r l } & { \rho ( E ) = \frac { 4 } { \pi ^ { 2 } } \left| \frac { E } { { \xi } ^ { 2 } } \right| ^ { 2 } \frac { 1 } { \sqrt { Z _ { 0 } } } \mathbf { F } \left( \frac { \pi } { 2 } \sqrt { \frac { E _ { 1 } } { Z _ { 0 } } } \right) } \\ & { Z _ { 0 = } \left\{ \begin{array} { l l } { \left( 1 + \left| \frac { E } { t } \right| \right) ^ { 2 } - \frac { \left( \left( \frac { E } { t } \right) ^ { 2 } - 1 \right) ^ { 2 } } { 4 } } & { ; \quad - t \le E \le t } \\ { \left( 1 + \left| \frac { E } { t } \right| \right) } & { ; \quad - 3 t \le E \le - t \lor t \le E \le 3 t } \end{array} \right. } \\ & { Z _ { 1 = } \left\{ \begin{array} { l l } { 4 \left| \frac { E } { t } \right| } & { ; \quad - t \le E \le t } \\ { \left( 1 + \left| \frac { E } { t } \right| \right) ^ { 2 } - \frac { \left( \left( \frac { E } { t } \right) ^ { 2 } - 1 \right) ^ { 2 } } { 4 } } & { ; \quad - 3 t \le E \le - t \lor t \le E \le 3 t } \end{array} \right. } \end{array}
$$

where $\mathbf { F } ( \pi / 2 , x )$ is the complete elliptic integral of the first kind. Close to the Dirac point the dispersion is approximated by (7) and the expression for the density of states per unit cell is given by (with a degeneracy of four included):

# B. Dirac fermions

Let us consider Hamiltonian (5) with $t ^ { \prime } = 0$ and consider the Fourier transform of the electron operators:

$$
\rho ( E ) = \frac { 2 A _ { c } } { \pi } \frac { | E | } { v _ { F } ^ { 2 } }
$$

$$
a _ { n } = \frac { 1 } { \sqrt { N _ { c } } } \sum _ { k } e ^ { - i \mathbf { k } \cdot \mathbf { R } _ { n } } a ( \mathbf { k } ) , ,
$$

where $N _ { c }$ is the number of unit cells. Using this transformation, let us write the field $a _ { n }$ as a sum of two terms, coming from expanding the Fourier sum around $\pmb { K } ^ { \prime }$ and $\kappa$ . This produces an approximation for the representation of the field $a _ { n }$ as a sum of two new fields, written as

$$
\begin{array} { r l r } { a _ { n } } & { \simeq } & { e ^ { - i { \bf K } \cdot { \bf R } _ { n } } a _ { 1 , n } + e ^ { - i { \bf K ^ { \prime } } \cdot { \bf R } _ { n } } a _ { 2 , n } , } \\ { b _ { n } } & { \simeq } & { e ^ { - i { \bf K } \cdot { \bf R } _ { n } } b _ { 1 , n } + e ^ { - i { \bf K ^ { \prime } } \cdot { \bf R } _ { n } } b _ { 2 , n } , } \end{array}
$$

where the index $i ~ = ~ 1$ $i \ = \ 2$ ) refers to the K (K') point. These new fields, $a _ { i , n }$ and $b _ { i , n }$ are assumed to vary slowly over the unit cell. The procedure for deriving a theory that is valid close to the Dirac point consists in using this representation in the tight-binding Hamiltonian and expanding the operators up to a linear order in $\delta$ . In the derivation one uses the fact that $\begin{array} { r } { \sum _ { \delta } e ^ { \pm i { \cal K } \cdot \delta } = \sum _ { \delta } e ^ { \pm i { \cal K } ^ { \prime } \cdot \delta } = 0 } \end{array}$ After some straightforward algebra we arrive at (Semenoff, 1984):

$$
\begin{array} { l } { { \displaystyle \simeq - t \int d x d y \hat { \Psi } _ { 1 } ^ { \dagger } ( r ) [ ( \begin{array} { c c } { { 0 } } & { { 3 a ( 1 - i \sqrt { 3 } ) / 4 } } \\ { { - 3 a ( 1 + i \sqrt { 3 } ) / 4 } } & { { 0 } } \end{array} ) \partial _ { x } + ( \begin{array} { c c } { { 0 } } & { { 3 a ( - i - i \sqrt { 3 } ) / 4 } } \\ { { - 3 a ( i - \sqrt { 3 } ) / 4 } } & { { 0 } } \end{array} ) \partial _ { x } + ( \begin{array} { c c } { { 0 } } & { { 3 a ( - i - i \sqrt { 3 } ) / 4 } } \\ { { - 3 a ( i - \sqrt { 3 } ) / 4 } } & { { 0 } } \end{array} ) \partial _ { x } + ( \begin{array} { c c } { { 0 } } & { { 3 a ( - i - i \sqrt { 3 } ) / 4 } } \\ { { - 3 a ( i - i \sqrt { 3 } ) / 4 } } & { { 0 } } \end{array} ) \partial _ { x } } } \\ { { + } } & { { 0 } } \\ { { = - i v _ { F } \int d x d y ( \hat { \Psi } _ { 1 } ^ { \dagger } ( r ) \sigma \cdot \nabla \hat { \Psi } _ { 1 } ( r ) + \hat { \Psi } _ { 2 } ^ { \dagger } ( r ) \sigma ^ { * } \cdot \nabla \hat { \Psi } _ { 2 } ( r ) ) , } } \end{array}
$$

with Pauli matrices $\boldsymbol { \sigma } = ( \sigma _ { x } , \sigma _ { y } )$ , $\pmb { \sigma } ^ { * } = ( \sigma _ { x } , - \sigma _ { y } )$ , and $\hat { \Psi } _ { i } ^ { \dagger } ~ = ~ ( a _ { i } ^ { \dagger } , b _ { i } ^ { \dagger } )$ b $i = 1 , 2$ ). It is clear that the effective Hamiltonian (18) is made of two copies of the massless Dirac-like Hamiltonian, one holding for $\mathbfcal { p }$ around $\kappa$ and other for $\mathbfcal { p }$ around $\pmb { K } ^ { \prime }$ . Notice that, in first quantized language, the two-component electron wavefunction, $\psi ( \mathbf { r } )$ , close to the K point, obeys the 2D Dirac equation:

$$
- i v _ { F } { \pmb \sigma } \cdot \nabla \psi ( { \bf r } ) = E \psi ( { \bf r } ) .
$$

The wavefunction, in momentum space, for the momentum around $\kappa$ has the form:

$$
\psi _ { \pm , { \bf K } } ( { \bf k } ) = \frac { 1 } { \sqrt { 2 } } \left( \begin{array} { c } { { e ^ { - i \theta _ { \bf k } / 2 } } } \\ { { \pm e ^ { i \theta _ { \bf k } / 2 } } } \end{array} \right) ,
$$

for ${ \cal H } _ { K } = v _ { F } \pmb { \sigma } \cdot { \bf k }$ , where the $\pm$ signs correspond to the eigenenergies $E = \pm v _ { F } k$ , that is, for the $\pi$ and $\pi ^ { * }$ band, respectively, and $\theta _ { \mathbf { k } }$ is given by (9). The wavefunction for the momentum around $\pmb { K } ^ { \prime }$ has the form:

$$
\psi _ { \pm , { \bf K } ^ { \prime } } ( { \bf k } ) = \frac { 1 } { \sqrt { 2 } } \left( \begin{array} { c } { { e ^ { i \theta _ { \bf k } / 2 } } } \\ { { \pm e ^ { - i \theta _ { \bf k } / 2 } } } \end{array} \right) ,
$$

for ${ \cal H } _ { K ^ { \prime } } \ = \ v _ { { \cal F } } \pmb { \sigma } ^ { * } \cdot { \bf k }$ . Notice that the wavefunctions at $\mathbf { K }$ and $\mathbf { K } ^ { \prime }$ are related by time reversal symmetry: if we set the origin of coordinates in momentum space in the M-point of the BZ (see Fig.2), time reversal becomes equivalent to a reflection along the $k _ { x }$ axis, that is, $( k _ { x } , k _ { y } )  ( k _ { x } , - k _ { y } )$ . Also note that if the phase $\theta$ is rotated by $2 \pi$ the wavefunction changes sign indicating a phase of $\pi$ (in the literature this is commonly called a Berry's phase). This change of phase by $\pi$ under rotation is characteristic of spinors. In fact, the wavefunction is a two component spinor.

A relevant quantity used to characterize the eigenfunctions is their helicity defined as the projection of the momentum operator along the (pseudo)spin direction. The quantum mechanical operator for the helicity has the form:

It is clear from the definition of $\hat { h }$ that the states $\psi _ { K } ( r )$ and $\psi _ { K ^ { \prime } } ( \pmb { r } )$ are also eigenstates of $\hat { h }$ :

$$
\hat { h } \psi _ { K } ( \pmb { r } ) = \pm \frac { 1 } { 2 } \psi _ { K } ( \pmb { r } ) ,
$$

and an equivalent equation for $\psi _ { K ^ { \prime } } ( \boldsymbol { r } )$ with inverted sign. Therefore electrons (holes) have a positive (negative) helicity. Equation (23) implies that $_ { \pmb { \sigma } }$ has its two eigenvalues either in the direction of (↑) or against $( \Downarrow )$ the momentum $\mathbfcal { p }$ . This property says that the states of the system close to the Dirac point have well defined chirality or helicity. Notice that chirality is not defined in regards to the real spin of the electron (that has not yet appeared in the problem) but to a pseudo-spin variable associated with the two components of the wavefunction. The helicity values are good quantum numbers as long as the Hamiltonian (18) is valid. Therefore the existence of helicity quantum numbers holds only as an asymptotic property, which is well defined close to the Dirac points $\kappa$ and $\pmb { K } ^ { \prime }$ . Either at larger energies or due to the presence of a finite $t ^ { \prime }$ the helicity stops being a good quantum number.

# Chiral Tunneling and Klein paradox

In this section we want to address the scattering of chiral electrons in two dimensions by a square barrier (Katsnelson, 2007b; Katsnelson et al., 2006). The one dimensional scattering of chiral electrons was discussed earlier in the context of carbon nanotubes (Ando et al., 1998; McEuen et al., 1999)

We start by noticing that by a gauge transformation the wavefunction (20) can be written as:

$$
\psi _ { K } ( k ) = { \frac { 1 } { \sqrt { 2 } } } \binom { 1 } { \pm e ^ { i \theta _ { k } } } \ .
$$

$$
{ \hat { h } } = { \frac { 1 } { 2 } } { \boldsymbol { \sigma } } \cdot { \frac { \boldsymbol { p } } { | { \boldsymbol { p } } | } } .
$$

We further assume that the scattering does not mix the momenta around $\kappa$ and $\pmb { K } ^ { \prime }$ points. In Fig. 6 we depict

![](images/a0db4fa6c3f90387ffd51a57edbeae489eca17ce4c87653010dfe39aae46ee5c.jpg)  
Figure 6 (Color online) Top: Schematic picture of the scattering of Dirac electrons by a square potential. Bottom: definition of the angles $\phi$ and $\theta$ used in the scattering formalism in the three regions I, II, and III.

the scattering process due to the square barrier of width $D$ .

The wavefunction in the different regions can be written in terms of incident and reflected waves. In region I we have:

$$
\begin{array} { l } { { \psi _ { I } ( r ) ~ = ~ { \frac { 1 } { \sqrt { 2 } } } \left( \begin{array} { c } { { 1 } } \\ { { s e ^ { i \phi } } } \end{array} \right) e ^ { i ( k _ { x } x + k _ { y } y ) } } } \\ { { ~ + ~ { \frac { r } { \sqrt { 2 } } } \left( \begin{array} { c } { { 1 } } \\ { { s e ^ { i ( \pi - \phi ) } } } \end{array} \right) e ^ { i ( - k _ { x } x + k _ { y } y ) } , } } \end{array}
$$

with $\phi = \arctan ( k _ { y } / k _ { x } )$ , $k _ { x } = k _ { F } \cos \phi$ , $k _ { y } = k _ { F } \sin \phi $ and $k _ { F }$ the Fermi momentum. In region II we have:

$$
\begin{array} { l } { { \psi _ { I I } ( { \bf r } ) ~ = ~ { \frac { a } { \sqrt { 2 } } } \left( \begin{array} { c } { { 1 } } \\ { { s ^ { \prime } e ^ { i \theta } } } \end{array} \right) e ^ { i ( q _ { x } x + k _ { y } y ) } } } \\ { { ~ + ~ { \frac { b } { \sqrt { 2 } } } \left( \begin{array} { c } { { 1 } } \\ { { s ^ { \prime } e ^ { i ( \pi - \theta ) } } } \end{array} \right) e ^ { i ( - q _ { x } x + k _ { y } y ) } , } } \end{array}
$$

with $\theta = \arctan ( k _ { y } / q _ { x } )$ and

$$
q _ { x } = \sqrt { ( V _ { 0 } - E ) ^ { 2 } / ( v _ { F } ^ { 2 } ) - k _ { y } ^ { 2 } } ,
$$

with $s = \mathrm { s g n } ( E )$ and $s ^ { \prime } = \mathrm { s g n } ( E - V _ { 0 } )$ . The coefficients $r , \ a$ , $b$ and $t$ are determined from the continuity of the wavefunction, which implies that the wavefunction has to obey the conditions $\psi _ { I } ( x = 0 , y ) = \psi _ { I I } ( x = 0 , y )$ and $\psi _ { I I } ( x = D , y ) = \psi _ { I I I } ( x = D , y )$ . Unlike the Schödinger equation we only need to match the wavefunction but not its derivative. The transmission through the barrier is obtained from $T ( \phi ) = t t ^ { * }$ and has the form:

![](images/4a464287962f2616874f4e47e109d4fe134dfeef95a98cdfbe1b7538f71d3e94.jpg)  
Figure 7 (Color online) Angular behavior of $T ( \phi )$ for two different values of $V _ { \mathrm { 0 } }$ . $V _ { 0 } = 2 0 0$ meV dashed line, $V _ { 0 } = 2 8 5$ meV solid line. The remaining parameters are $D = 1 1 0 \ \mathrm { n m }$ (top), $D = 5 0$ nm (bottom) $E = 8 0 ~ \mathrm { m e V }$ , $k _ { F } = 2 \pi / \lambda$ , $\lambda = 5 0$ nm.

and finally in region III we have a transmitted wave only:

$$
\psi _ { I I I } ( \pmb { r } ) = \frac { t } { \sqrt { 2 } } \left( \begin{array} { c } { 1 } \\ { s e ^ { i \phi } } \end{array} \right) e ^ { i ( k _ { x } x + k _ { y } y ) } ,
$$

$$
T ( \phi ) = \frac { \cos ^ { 2 } \theta \cos ^ { 2 } \phi } { [ \cos ( D q _ { x } ) \cos \phi \cos \theta ] ^ { 2 } + \sin ^ { 2 } ( D q _ { x } ) ( 1 - s s ^ { \prime } \sin \phi \sin \theta ) ^ { 2 } } .
$$

This expression does not take into account a contribution from evanescent waves in region II, which is usually negligible, unless the chemical potential in region $\mathrm { I I }$ is at the Dirac energy (see section IV.I).

Notice that $T ( \phi ) ~ = ~ T ( - \phi )$ and for values of $D q _ { x }$ satisfying the relation $\begin{array} { l l l } { \Delta q _ { x } } & { = } & { n \pi } \end{array}$ , with $n$ an integer, the barrier becomes completely transparent since $T ( \phi ) = 1$ , independently of the value of $\phi$ . Also, for normal incidence ( $\phi ~  ~ 0$ and $\theta ~  ~ 0$ ) and for any value of $D q _ { x }$ one obtains $T ( 0 ) ~ = ~ 1$ , and the barrier is again totally transparent. This result is a manifestation of the Klein paradox (Calogeracos and Dombey, 1999; Itzykson and Zuber, 2006) and does not occur for non-relativistic electrons. In this latter case and for normal incidence, the transmission is always smaller than one. In the limit $| V _ { 0 } | \gg | E |$ , eq. (29) has the following asymptotic form

$$
T ( \phi ) \simeq \frac { \cos ^ { 2 } \phi } { 1 - \cos ^ { 2 } ( D q _ { x } ) \sin ^ { 2 } \phi } .
$$

In Fig. 7 we show the angular dependence of $T ( \phi )$ for two different values of the potential $V _ { 0 }$ ; it is clear that there are several directions for which the transmission is one. Similar calculations were done for a graphene bilayer (Katsnelson et al., 2006) with its most distinctive behavior being the absence of tunneling in the forward ( $k _ { y } = 0$ ) direction.

The simplest example of a potential barrier is a square potential discussed previously. When intervalley scattering and the lack of symmetry between sublattices are neglected, a potential barrier shows no reflection for electrons incident in the normal direction (Katsnelson et al., 2006). Even when the barrier separates regions where the Fermi surface is electron like on one side and hole like on the other, a normally incident electron continues propagating as a hole with $1 0 0 \%$ efficiency. This phenomenon is another manifestation of the chirality of the Dirac electrons within each valley, which prevents backscattering in general. The transmission and reflection probabilities of electrons at different angles depend on the potential profile along the barrier. A slowly varying barrier is more efficient in reflecting electrons at non-zero incident angles (Cheianov and Fal'ko, 2006).

Electrons moving through a barrier separating p- and n-doped graphene, a p-n junction, are transmitted as holes. The relation between the velocity and the momentum for a hole is the inverse of that for an electron. This implies that, if the momentum parallel to the barrier is conserved, the velocity of the quasiparticle is inverted. When the incident electrons emerge from a source, the transmitting holes are focused into an image of the source. This behavior is the same as that of photons moving in a medium with negative reflection index (Cheianov et al., 2007a). Similar effects can occur in graphene quantum dots, where the inner and outer regions contain electrons and holes, respectively (Cserti et al., 2007b). Note that the fact that barriers do not impede the transmission of normally incident electrons does not preclude the existence of sharp resonances, due to the confinement of electrons with a finite parallel momentum. This leads to the possibility of fabricating quantum dots with potential barriers (Silvestrov and Efetov, 2007). Finally, at half-filling, due to disorder graphene can be divided in electron and hole charge puddles (Katsnelson et al., 2006; Martin et al., 2007). Transport is determined by the transmission across the p-n junctions between these puddles (Cheianov et al., 2007b; Shklovskii, 2007). There is a rapid progress in the measurement of transport properties of graphene ribbons with additional top gates that play the role of tunable potential barriers (Han et al., 2007; Huard et al., 2007; Lemme et al., 2007; Özyilmaz et al., 2007; Williams et al., 2007).

A magnetic field and potential fluctuations break both inversion symmetry of the lattice and time reversal symmetry. The combination of these effects break also the symmetry between the two valleys. The transmission coefficient becomes valley dependent, and, in general, electrons from different valleys propagate along different paths. This opens the possibility of manipulating the valley index (Tworzydlo et al., 2007) (valleytronics) in a way similar to the control of the spin in mesoscopic devices (spintronics). For large magnetic fields, a p-n junction separates regions with different quantized Hall conductivities. At the junction, chiral currents can flow at both edges (Abanin and Levitov, 2007), inducing backscattering between the Hall currents at the edges of the sample.

The scattering of electrons near the Dirac point by graphene-superconductor junctions differs from Andreev scattering process in normal metals (Titov and Beenakker, 2006). When the distance between the Fermi energy and the Dirac energy is smaller than the superconducting gap, the superconducting interaction hybridizes quasiparticles from one band with quasiholes in the other. As in the case of scattering at a p-n junction, the trajectories of the incoming electron and reflected hole (note that hole here is meant as in the BCS theory of superconductivity) are different from those in similar processes in metals with only one type of carrier (Bhattacharjee and Sengupta, 2006; Maiti and Sengupta, 2007).

# Confinement and zitterbewegung

Zitterbewegung, or jittery motion of the wavefunction of the Dirac problem, occurs when one tries to confine the Dirac electrons (Itzykson and Zuber, 2006). Localization of a wavepacket leads, due to the Heisenberg principle, to uncertainty in the momentum. For a Dirac particle with zero rest mass, uncertainty in the momentum translates into uncertainty in the energy of the particle as well (this should be contrasted with the non-relativistic case where the position-momentum uncertainty relation is independent of the energy-time uncertainty relation). Thus, for a ultra-relativistic particle, a particle-like state can have hole-like states in its time evolution. Consider, for instance, if one tries to construct a wave packet at some time $t = 0$ , and let us assume, for simplicity, that this packet has a Gaussian shape of width $w$ with momentum close to $\mathbf { K }$ :

$$
\psi _ { 0 } ( \mathbf { r } ) = { \frac { e ^ { - r ^ { 2 } / ( 2 w ^ { 2 } ) } } { \sqrt { \pi } w } } e ^ { i \mathbf { K } \cdot \mathbf { r } } \phi ,
$$

where $\phi$ is spinor composed of positive energy states (associated with $\psi _ { + , \mathbf { K } }$ of (20)). The eigenfunction of the Dirac equation can be written in terms of the solution (20) as:

$$
\psi ( \mathbf { r } , t ) = \int { \frac { d ^ { 2 } k } { ( 2 \pi ) ^ { 2 } } } \sum _ { a = \pm 1 } \alpha _ { a , \mathbf { k } } \psi _ { a , \mathbf { K } } ( \mathbf { k } ) e ^ { - i a ( \mathbf { k } \cdot \mathbf { r } + v _ { F } k t ) } ,
$$

where $\alpha _ { \pm , \mathbf { k } }$ are Fourier coefficients. We can rewrite (31) in terms of (32) by inverse Fourier transform and find that:

$$
\alpha _ { \pm , { \bf k } } = \sqrt { \pi } w e ^ { - k ^ { 2 } w ^ { 2 } / 2 } \psi _ { \pm , { \bf K } } ^ { \dagger } ( { \bf k } ) \phi .
$$

Notice that the relative weight of positive energy states with respect to negative energy states, $| \alpha _ { + } / \alpha _ { - } |$ , given by (20) is one, that is, there are as many positive energy states as negative energy states in a wavepacket. Hence, these will cause the wavefunction to be delocalized at any time $t \neq 0$ . Thus, a wave packet of electron-like states has hole-like components, a result that puzzled many researchers in the early days of QED (Itzykson and Zuber, 2006).

![](images/df5b3ccc3f118e3ae04201334124ccc61f5a731c5a1ed3a45324609a93ad0406.jpg)  
Figure 8 (Color on line) Energy spectrum (in units of $t$ )for a graphene ribbon $6 0 0 a$ wide, as a function of the momentum $k$ along the ribbon (in units of $1 / ( \sqrt { 3 } a ) )$ , in the presence of confining potential with $V _ { 0 } = 1$ eV, $\lambda = 1 8 0 a$ .

Consider the tight-binding description (Chen et al., 2007a; Peres et al., 2006b) of Sec. II.A when a potential $V _ { i }$ on site $\mathbf { R } _ { i }$ is added to the problem:

$$
H _ { e } = \sum _ { i } V _ { i } n _ { i } ,
$$

where $n _ { i }$ is the local electronic density. For simplicity, we assume that the confining potential is 1D, that is, that $V _ { i }$ vanishes in the bulk but becomes large at the edge of the sample. Let us assume a potential that decays exponentially away from the edges into the bulk with a penetration depth, $\lambda$ . In Fig. 8 we show the electronic spectrum for a graphene ribbon of width $L = 6 0 0 a$ , in the presence of a confining potential,

$$
V ( x ) = V _ { 0 } \left[ e ^ { - ( x - L / 2 ) / \lambda } + e ^ { - ( L / 2 - x ) / \lambda } \right] ,
$$

where $x$ is the direction of confinement and $V _ { 0 }$ the strength of the potential. One can clearly see that in the presence of the confining potential the electron-hole symmetry is broken and, for $V _ { 0 } ~ > ~ 0$ , the hole part of the spectrum is strongly distorted. In particular, for $k$ close to the Dirac point, we see that the hole dispersion is given by: $E _ { n , \sigma = - 1 } ( k ) \approx - \gamma _ { n } k ^ { 2 } - \zeta _ { n } k ^ { 4 }$ where $n$ is a positive integer, and $\gamma _ { n } ~ < ~ 0$ b $\gamma _ { n } > 0 )$ for $n < N ^ { * }$ $( n > N ^ { * } )$ . Hence, at $n = N ^ { * }$ the hole effective mass $d i$ . verges ( $\gamma _ { N ^ { * } } = 0 \mathrm { . }$ ) and, by tuning the chemical potential, $\mu$ , via a back gate, to the hole region of the spectrum ( $\mu < 0 )$ one should be able to observe an anomaly in the Shubnikov-de Haas (SdH) oscillations. This is how zitterbewegung could manifest itself in magnetotransport.

# C. Bilayer graphene: tight-binding approach

The tight-binding model developed for graphite can be easily extended to stacks with a finite number of graphene layers. The simplest generalization is a bilayer (McCann and Fal'ko, 2006). A bilayer is interesting because the IQHE shows anomalies, although different from those observed in a single layer (Novoselov et al., 2006), and also a gap can open between the conduction and valence band (McCann and Fal'ko, 2006). The bilayer structure, with the AB stacking of 3D graphite, is shown in Fig.9.

![](images/6d8538b7594df1caff4262c37c4d35f06a456fc8bd02d897899efbbd20a467c2.jpg)  
Figure 9 (Color online)(a) Lattice structure of the bilayer with the various hopping parameters according to the SWM model. The A-sublattices are indicated by the darker spheres. (b) Brillouin zone. Adapted from Malard et al., 2007.

The tight-binding Hamiltonian for this problem can be written as:

$$
\begin{array} { l } { \mathcal { H } _ { \mathrm { t , b . } } ~ = ~ - \gamma _ { 0 } \displaystyle \sum _ { \omega , j , \sigma } \big ( a _ { m , i , \sigma } ^ { \dagger } b _ { m , j , \sigma } + \mathrm { h . c . } \big ) } \\ { ~ - ~ \gamma _ { 1 } \displaystyle \sum _ { j , \sigma } \big ( a _ { 1 , j , \sigma } ^ { \dagger } a _ { 2 , j , \sigma } + \mathrm { h . c . } \big ) , } \\ { ~ - ~ \gamma _ { 3 } \displaystyle \sum _ { j , \sigma } \big ( a _ { 1 , j , \sigma } ^ { \dagger } b _ { 2 , j , \sigma } + a _ { 2 , j , \sigma } ^ { \dagger } b _ { 1 , j , \sigma } + \mathrm { h . c . } \big ) } \\ { ~ - ~ \gamma _ { 4 } \displaystyle \sum _ { j , \sigma } \big ( b _ { 1 , j , \sigma } ^ { \dagger } b _ { 2 , j , \sigma } + \mathrm { h . c . } \big ) , } \end{array}
$$

where $a _ { m , i , \sigma }$ $\left( b _ { m , i \sigma } \right)$ annihilates an electron with spin $\sigma$ , on sublattice A (B), in plane $m = 1 , 2$ , at site $\mathbf { R } _ { i }$ . Here we use the graphite nomenclature for the hopping parameters: $\gamma _ { 0 } = t$ is the in-plane hopping energy and $\gamma _ { 1 }$ b $\gamma _ { 1 } = t _ { \perp } \approx 0 . 4 ~ \mathrm { e V }$ in graphite (Brandt et al., 1988; Dresselhaus and Dresselhaus, 2002)) is the hopping energy between atom $\mathrm { A } _ { 1 }$ and atom $\mathrm { A _ { 2 } }$ (see Fig. 9), and $\gamma _ { 3 }$ b $\mathrm { \Omega \gamma _ { 3 } \approx 0 . 3 ~ e V }$ in graphite (Brandt et al., 1988; Dresselhaus and Dresselhaus, 2002)) is the hopping energy between atom $\mathrm { A } _ { 1 }$ ′ $\left( \mathrm { A } _ { 2 } \right)$ and atom $\mathrm { B _ { 2 } }$ ′ $\left( \operatorname { B } _ { 1 } \right)$ , and $\gamma _ { 4 }$ b $\gamma _ { 4 } ~ \approx ~ - 0 . 0 4 ~ \mathrm { e V }$ in graphite (Brandt et al., 1988; Dresselhaus and Dresselhaus, 2002)) that connects $\mathrm { B _ { 1 } }$ and $\mathrm { B _ { 2 } }$ .

In the continuum limit, by expanding the momentum close to the $\mathrm { K }$ point in the BZ, the Hamiltonian reads,

$$
\mathcal { H } = \sum _ { \mathbf { k } } \Psi _ { \mathbf { k } } ^ { \dagger } \cdot \mathcal { H } _ { K } \cdot \Psi _ { \mathbf { k } }
$$

![](images/ae8818fda78919150d100bf67531766e866c97b0a6f08c453c999b3e5e378a1e.jpg)  
Figure 10 (Color online) Band structure for bilayer graphene for $V = 0$ and $\gamma _ { 3 } = 0$ .

where (ignoring $\gamma _ { 4 }$ for the time being):

$$
\mathcal { H } _ { K } \equiv \left( \begin{array} { c c c c } { - V } & { v _ { F } k } & { 0 } & { 3 \gamma _ { 3 } a k ^ { * } } \\ { v _ { F } k ^ { * } } & { - V } & { \gamma _ { 1 } } & { 0 } \\ { 0 } & { \gamma _ { 1 } } & { V } & { v _ { F } k } \\ { 3 \gamma _ { 3 } a k } & { 0 } & { v _ { F } k ^ { * } } & { V } \end{array} \right) ~ ,
$$

where $k = k _ { x } + i k _ { y }$ is a complex number, and we have added $V$ which is here half the shift in electro-chemical potential between the two layers (this term will appear if a potential bias is applied between the layers), and

$$
\Psi _ { \mathbf { k } } ^ { \dagger } = \left( a _ { 1 } ^ { \dagger } ( \mathbf { k } ) , a _ { 2 } ^ { \dagger } ( \mathbf { k } ) , b _ { 1 } ^ { \dagger } ( \mathbf { k } ) , b _ { 2 } ^ { \dagger } ( \mathbf { k } ) \right)
$$

is a four component spinor.

If $V ~ = ~ 0$ and $\gamma _ { 3 } , v _ { F } k \ll \gamma _ { 1 }$ , one can eliminate the high energy states perturbatively and write an effective Hamiltonian:

$$
\mathcal { H } _ { K } \equiv \left( \begin{array} { c c c } { { 0 } } & { { \frac { v _ { F } ^ { 2 } k ^ { 2 } } { \gamma _ { 1 } } + 3 \gamma _ { 3 } a k ^ { * } } } \\ { { \frac { v _ { F } ^ { 2 } ( k ^ { * } ) ^ { 2 } } { \gamma _ { 1 } } + 3 \gamma _ { 3 } a k } } & { { 0 } } \end{array} \right) .
$$

![](images/1d96faea4ba489e3f12ba10ff8b9460fe9a83281464dd2ebc755e511e9a4d36e.jpg)  
Figure 11 (Color online) Band structure for bilayer graphene for $V \neq 0$ and $\gamma _ { 3 } = 0$ .

The hopping $\gamma _ { 4 }$ leads to a $k$ dependent coupling between the sublattices or a small renormalization of $\gamma _ { 1 }$ . The same role is played by the inequivalence between sublattices within a layer.

For $\gamma _ { 3 } ~ = ~ 0$ , (40) gives two parabolic bands, $\epsilon _ { k , \pm } \approx$ $\pm v _ { F } ^ { 2 } k ^ { 2 } / t _ { \perp }$ which touch at $\epsilon = 0$ (as shown in Fig.10). The spectrum is electron-hole symmetric. There are two additional bands which start at $\pm t _ { \perp }$ . Within this approximation, the bilayer is metallic, with a constant density of states. The term $\gamma _ { 3 }$ changes qualitatively the spectrum at low energies since it introduces a trigonal distortion, or warping, of the bands (notice that this trigonal distortion, unlike the one introduced by large momentum in (8), occurs at low energies). The electron-hole symmetry is preserved but, instead of two bands touching at $k = 0$ , we obtain three sets of Dirac-like linear bands. One Dirac point is at $\epsilon = 0$ and $k = 0$ , while the three other Dirac points, also at $\epsilon = 0$ , lie at three equivalent points with a finite momentum. The stability of points where bands touch can be understood using topological arguments (Mañes et al., 2007). The winding number of a closed curve in the plane around a given point is an integer representing the total number of times that the curve travels counterclockwise around the point so that the wavefunction remains unaltered. The winding number of the point where the two parabolic bands come together for $\gamma _ { 3 } = 0$ has winding number $+ 2$ . The trigonal warping term, $\gamma _ { 3 }$ , splits it into a Dirac point at $k = 0$ and winding number $^ { - 1 }$ , and three Dirac points at $k \neq 0$ and winding numbers $+ 1$ . An in-plane magnetic field, or a small rotation of one layer with respect to the other splits the $\gamma _ { 3 } = 0$ degeneracy into two Dirac points with winding number $+ 1$ .

The term $V$ in (38) breaks the equivalence of the two layers, or, alternatively, inversion symmetry. In this case, the dispersion relation becomes:

$$
\begin{array} { l } { { \epsilon _ { \pm , \mathbf k } ^ { 2 } = V ^ { 2 } + v _ { F } ^ { 2 } k ^ { 2 } + t _ { \perp } ^ { 2 } / 2 } } \\ { { \pm \sqrt { 4 V ^ { 2 } v _ { F } ^ { 2 } k ^ { 2 } + t ^ { 2 } v _ { F } ^ { 2 } k ^ { 2 } + t _ { \perp } ^ { 4 } / 4 } , } } \end{array}
$$

given rise to the dispersion shown in Fig. 11, and to the opening of a gap close, but not directly at, the K point. For small momenta, and $V \ll t$ , the energy of the conduction band can be expanded:

$$
\epsilon _ { k } \approx V - ( 2 V v _ { F } ^ { 2 } k ^ { 2 } ) / t _ { \perp } + ( v _ { F } ^ { 4 } k ^ { 4 } ) / ( 2 t _ { \perp } ^ { 2 } V ) .
$$

The dispersion for the valence band can be obtained by replacing $\epsilon _ { k }$ by $- \epsilon _ { k }$ . The bilayer has a gap at $k ^ { 2 } \approx$ $( 2 V ^ { 2 } ) / v _ { F } ^ { 2 }$ . Notice, therefore, that the gap in the biased bilayer system depends on the applied bias and hence can be measured experimentally (Castro et al., 2007a; McCann, 2006; McCann and Fal'ko, 2006). The ability to open a gap makes bilayer graphene most interesting for technological applications.

# D. Epitaxial graphene

It has been known for a long time that monolayers of graphene could be grown epitaxially on metal surfaces by using catalytic decomposition of hydrocarbons or carbon oxide (Campagnoli and Tosatti, 1989;

Eizenberg and Blakely, 1979; Oshima and Nagashima, 1997; Shelton et al., 1974; Sinitsyna and Yaminsky, 2006). When such surfaces are heated, oxygen or hydrogen desorbs, and the carbon atoms form a graphene monolayer. The resulting graphene structures could reach sizes up to a micrometer, with few defects and were characterized by different surface-science techniques and local scanning probes (Himpsel et al., 1982). For example, graphene grown on ruthenium has zigzag edges and also ripples associated with a ( $1 0 \times 1 0$ )reconstruction (Vázquez de Parga et al., 2007).

Graphene can also be formed on the surface of SiC. Upon heating, the silicon from the top layers desorbs, and a few layers of graphene are left on the surface (Berger et al., 2004; Bommel et al., 1975; Coey et al., 2002; Forbeaux et al., 1998; Hass et al., 2007a; de Heer et al., 2007; Rollings et al., 2005). The number of layers can be controlled by limiting time or temperature of the heating treatment. The quality and the number of layers in the samples depends on the SiC face used for their growth (de Heer et al., 2007) (the carbon terminated surface produces few layers but with a low mobility, whereas the silicon terminated surface produces several layers but with higher mobility). Epitaxially grown multilayers exhibit SdH oscillations with a Berry phase shift of $\pi$ (Berger et al., 2006), which is the same as the phase shift for Dirac fermions observed in a single layer as well as for some subbands present in multilayer graphene (see further) and graphite (Luk'yanchuk and Kopelevich, 2004). The carbon layer directly on top of the substrate is expected to be strongly bonded to it, and it shows no $\pi$ bands (Varchon et al., 2007). The next layer shows a $( 6 { \sqrt { 3 } } \times 6 { \sqrt { 3 } } )$ reconstruction due to the substrate, and has graphene properties. An alternate route to produce few layers graphene is based on synthesis from nanodiamonds (Affoune et al., 2001).

Angle resolved photo-emission experiments (ARPES) show that epitaxial graphene grown on SiC has linearly dispersing quasiparticles (Dirac fermions) (Bostwick et al., 2007b; Ohta et al., 2007; Zhou et al., 2006b), in agreement with the theoretical expectation. Nevertheless, these experiments show that the electronic properties can change locally in space indicating a certain degree of inhomogeneity due to the growth method (Zhou et al., 2007). Similar inhomogeneities due to disorder in the c-axis orientation of graphene planes is observed in graphite (Zhou et al., 2006a). Moreover, graphene grown this way is heavily doped due to the charge transfer from the substrate to the graphene layer (with the chemical potential well above the Dirac point) and therefore all samples have strong metallic character with large electronic mobilities (Berger et al., 2006; de Heer et al., 2007). There is also evidence for strong interaction between a substrate and the graphene layer leading to the appearance of gaps at the Dirac point (Zhou et al., 2007). Indeed, gaps can be generated by the breaking of the sublattice symmetry and, as in the case of other carbon based systems such as polyacethylene (Su et al., 1979, 1980), it can lead to soliton-like excitations (Hou et al., 2007; Jackiw and Rebbi, 1976). Multilayer graphene grown on SiC have also been studied with ARPES (Bostwick et al., 2007a; Ohta et al., 2007, 2006) and the results seem to agree quite well with band structure calculations (Mattausch and Pankratov, 2007). Spectroscopy measurements also show the transitions associated with Landau levels (Sadowski et al., 2006), and weak localization effects at low magnetic fields, also expected for Dirac fermions (Wu et al., 2007). Local probes reveal a rich structure of terraces (Mallet et al., 2007) and interference patterns due to defects at or below the graphene layers (Rutter et al., 2007).

![](images/e8fcdb43815863d0a36e2369f7cc06be76cdb9ac95bd936ddf187b5b5c3e1d95.jpg)  
Figure 12 (Color online) Sketch of the three inequivalent orientations of graphene layers with respect to each other.

# E. Graphene stacks

In stacks with more than one graphene layer, two consecutive layers are normally oriented in such a way that the atoms in one of the two sublattices, $A _ { n }$ , of the honeycomb structure of one layer are directly above one half of the atoms in the neighboring layer, sublattice $A _ { n \pm 1 }$ . The second set of atoms in one layer sits on top of the (empty) center of an hexagon in the other layer. The shortest distance between carbon atoms in different layers is $d _ { A _ { n } A _ { n \pm 1 } } = c = 3 . 4 \mathrm { \AA }$ . The next distance is $d _ { A _ { n } B _ { n \pm 1 } } =$ $\sqrt { c ^ { 2 } + a ^ { 2 } }$ . This is the most common arrangement of nearest neighbor layers observed in nature, although a stacking order in which all atoms in one layer occupy positions directly above the atoms in the neighboring layers (hexagonal stacking) has been considered theoretically (Charlier et al., 1991) and appears in graphite intercalated compounds (Dresselhaus and Dresselhaus, 2002).

The relative position of two neighboring layers allows for two different orientations of the third layer. If we label the positions of the two first atoms as 1 and 2, the third layer can be of type 1, leading to the sequence 121, or it can fill a third position different from 1 and 2 (see Fig. 12), labeled 3. There are no more inequivalent positions where a new layer can be placed, so that thicker stacks can be described in terms of these three orientations. In the most common version of bulk graphite the stacking order is 1212... (Bernal stacking). Regions with the stacking $1 2 3 1 2 3 \cdots$ (rhombohedral stacking) have also been observed in different types of graphite (Bacon, 1950; Gasparoux, 1967). Finally, samples with no discernible stacking order (turbostratic graphite) are also commonly reported.

![](images/9d9e467b6137a863b99ca6329b31cca9a4845faa9247d6cead48be0731f5a17b.jpg)  
Figure 13 (Color online) Electronic bands of graphene multi layers: top left: biased bilayer; top right: trilayer with Bernal stacking; bottom left: trilayer with orthorhombic stacking; bottom right: stack with four layers where the top and bottom layers are shifted in energy with respect to the two middle layers by $+ 0 . 1$ eV.

Beyond two layers, the stack ordering can be arbitrarily complex. Simple analytical expressions for the electronic bands can be obtained for perfect Bernal ( $1 2 1 2 \cdots { } )$ and rhombohedral ( $1 2 3 1 2 3 \cdots )$ stacking (Guinea et al., 2006). Even if we consider one interlayer hopping, $t _ { \perp } = \gamma _ { 1 }$ , the two stacking orders show rather different band structures near $\epsilon = 0$ . A Bernal stack with $N$ layers, $N$ even, has $N / 2$ electron like and $N / 2$ hole like parabolic subbands touching at $\epsilon = 0$ . When $N$ is odd, an additional subband with linear (Dirac) dispersion emerges. Rhombohedral systems have only two subbands that touch at $\epsilon = 0$ . These subbands disperse as $k ^ { N }$ , and become surface states localized at the top and bottom layer when $N  \infty$ . In this limit, the remaining $2 N - 2$ subbands of a rhombohedral stack become Dirac like, with the same Fermi velocity as a single graphene layer. The subband structure of a tri-layer with the Bernal stacking includes two touching parabolic bands, and one with Dirac dispersion, combining the features of bilayer and monolayer graphene.

The low energy bands have different weights on the two sublattices of each graphene layer. The states at a site directly coupled to the neighboring planes are pushed to energies $\epsilon \approx \pm t _ { \perp }$ . The bands near $\epsilon = 0$ are localized mostly at the sites without neighbors in the next layers. For the Bernal stacking, this feature implies that the density of states at $\epsilon = 0$ at sites without nearest neighbors in the contiguous layers is finite, while it vanishes linearly at the other sites. In stacks with rhombohedral stacking, all sites have one neighbor in another plane, and the density of states vanishes at $\epsilon = 0$ (Guinea et al., 2006). This result is consistent with the well known fact that only one of the two sublattices at a graphite surface can be resolved by scanning tunneling microscopy (STM) (Tománek et al., 1987).

As in the case of a bilayer, an inhomogeneous charge distribution can change the electrostatic potential in the different layers. For more than two layers, this breaking of the equivalence between layers can take place even in the absence of an applied electric field. It is interesting to note that a gap can open in a stack with Bernal ordering and four layers, if the electronic charge at the two surface layers is different from that at the two inner ones. Systems with a higher number of layers do not show a gap, even in the presence of charge inhomogeneity. Four representative examples are shown in Fig. 13. The band structure analyzed here will be modified by the inclusion of the trigonal warping term, $\gamma _ { 3 }$ . Experimental studies of graphene stacks have showed that, with increasing number of layers, the system becomes increasingly metallic (concentration of charge carriers at zero energy gradually increases), and there appear several types of electron-and-hole-like carries (Morozov et al., 2005; Novoselov et al., 2004). An inhomogeneous charge distribution between layers becomes very important in this case, leading to 2D electron and hole systems that occupy only a few graphene layers near the surface and can completely dominate transport properties of graphene stacks (Morozov et al., 2005).

The degeneracies of the bands at $\epsilon = 0$ can be studied using topological arguments (Mañes et al., 2007). Multilayers with an even number of layers and Bernal stacking have inversion symmetry, leading to degeneracies with winding number $+ 2$ , as in the case of a bilayer. The trigonal lattice symmetry implies that these points can lead, at most, to four Dirac points. In stacks with an odd number of layers, these degeneracies can be completely removed. The winding number of the degeneracies found in stacks with $N$ layers and orthorhombic ordering is $\pm N$ The inclusion of trigonal warping terms will lead to the existence of many weaker degeneracies near $\epsilon = 0$ .

Furthermore, it is well known that in graphite the planes can be rotated relative each other giving rise to Moiré patterns that are observed in STM of graphite surfaces (Rong and Kuiper, 1993). The graphene layers can be rotated relative to each other due to the weak coupling between planes that allows for the presence of many different orientational states that are quasidegenerate in energy. For certain angles the graphene layers become commensurate with each other leading to a lowering of the electronic energy. Such phenomenon is quite similar to the commensurate-incommensurate transitions observed in certain charge density wave systems or adsorption of gases on graphite (Bak, 1982). This kind of dependence of the electronic structure on the relative rotation angle between graphene layers leads to what is called superlubricity in graphite (Dienwiebel et al., 2004), namely, the vanishing of the friction between layers as a function of the angle of rotation. In the case of bilayer graphene, a rotation by a small commensurate angle leads to the effective decoupling between layers and the recovery of the linear Dirac spectrum of the single layer albeit with a modification on the value of the Fermi velocity (Lopes dos Santos et al., 2007).

<table><tr><td rowspan=1 colspan=1>γ0</td><td rowspan=1 colspan=1>3.16 eV</td></tr><tr><td rowspan=1 colspan=1>γ1</td><td rowspan=1 colspan=1>0.39 eV</td></tr><tr><td rowspan=1 colspan=1>γ2</td><td rowspan=1 colspan=1>-0.020 eV</td></tr><tr><td rowspan=1 colspan=1>γ3</td><td rowspan=1 colspan=1>0.315 eV</td></tr><tr><td rowspan=1 colspan=1>γ4</td><td rowspan=1 colspan=1>γ4-0.044 eV</td></tr><tr><td rowspan=1 colspan=1>γ5</td><td rowspan=1 colspan=1>0.038 eV</td></tr><tr><td rowspan=1 colspan=1>Δ</td><td rowspan=1 colspan=1>-0.008 eV</td></tr></table>

# 1. Electronic structure of bulk graphite

The tight-binding description of graphene described earlier can be extended to systems with an infinite number of layers. The coupling between layers leads to hopping terms between $\pi$ orbitals in different layers, leading to the so called Slonczewski-Weiss-McClure model (Slonczewski and Weiss, 1958). This model describes the band structure of bulk graphite with the Bernal stacking order in terms of seven parameters, $\gamma _ { 0 } , \gamma _ { 1 } , \gamma _ { 2 } , \gamma _ { 3 } , \gamma _ { 4 } , \gamma _ { 5 }$ and $\Delta$ . The parameter $\gamma _ { 0 }$ describes the hopping within each layer, and it has been considered previously. The coupling between orbitals in atoms that are nearest neighbors in successive layers is $\gamma _ { 1 }$ , which we called $t _ { \perp }$ earlier. The parameters $\gamma _ { 3 }$ and $\gamma _ { 4 }$ describe the hopping between orbitals at next nearest neighbors in successive layers and were discussed in the case of the bilayer. The coupling between orbitals at next nearest neighbor layers are $\gamma _ { 2 }$ and $\gamma _ { 5 }$ . Finally, $\Delta$ is an on site energy which reflects the inequivalence between the two sublattices in each graphene layer once the presence of neighboring layers is taken into account. The values of these parameters, and their dependence with pressure, or, equivalently, the interatomic distances, have been extensively studied (Brandt et al., 1988; Dillon et al., 1977; Dresselhaus and Mavroides, 1964; McClure, 1957, 1964; Nozières, 1958; Soule et al., 1964). A representative set of values is shown in Table[I]. It is unknown, however, how these parameters may vary in graphene stacks with a small number of layers.

The unit cell of graphite with Bernal stacking includes two layers, and two atoms within each layer. The tightbinding Hamiltonian described previously can be represented as a $4 \times 4$ matrix. In the continuum limit, the two inequivalent corners of the BZ can be treated separately, and the in plane terms can be described by the Dirac equation. The next terms in importance for the low energy electronic spectrum are the nearest neighbor couplings $\gamma _ { 1 }$ and $\gamma _ { 3 }$ . The influence of the parameter $\gamma _ { 4 }$ on the low energy bands is much smaller, as discussed below. Finally, the fine details of the spectrum of bulk graphite are determined by $\Delta$ , which breaks the electronhole symmetry of the bands preserved by $\gamma _ { 0 } , \gamma _ { 1 }$ and $\gamma _ { 3 }$ It is usually assumed to be much smaller than the other terms.

We label the two atoms from the unit cell in one layer as 1 and 2, and 3 and 4 correspond to the second layer. Atoms 2 and 3 are directly on top of each other. Then, the matrix elements of the Hamiltonian can be written

Table I Band structure parameters of graphite (Dresselhaus and Dresselhaus, 2002).

as:

$$
\begin{array} { r l } { H _ { 1 1 } ^ { K } } & { = 2 \gamma _ { 2 } \cos ( 2 \pi \tilde { x } _ { 2 } / c ) } \\ { H _ { 1 2 } ^ { I } } & { = \nu _ { F } ( k _ { x } + i k _ { y } ) } \\ { H _ { 1 3 } ^ { K } } & { = \frac { 3 \gamma _ { 4 } \alpha } { 2 } \left( 1 + e ^ { i k _ { x } c } \right) \left( k _ { x } + i k _ { y } \right) } \\ { H _ { 1 4 } ^ { K } } & { = \frac { 3 \gamma _ { 1 } \alpha } { 2 } \left( 1 + e ^ { i k _ { x } c } \right) \left( k _ { x } - i k _ { y } \right) } \\ { H _ { 2 2 } ^ { K } } & { = \Delta + 2 \gamma _ { 1 5 } \cos ( 2 \pi k _ { x } / c ) } \\ { H _ { 2 3 } ^ { K } } & { = \nu _ { 1 } \left( 1 + e ^ { i k _ { x } c } \right) } \\ { H _ { 2 4 } ^ { I } } & { = \frac { 3 \gamma _ { 2 } \alpha } { 2 } \left( 1 + e ^ { i k _ { x } c } \right) \left( k _ { x } + i k _ { y } \right) } \\ { H _ { 3 4 } ^ { I } } & { = \Delta + 2 \gamma _ { 1 6 } \cos ( 2 \pi k _ { x } / c ) } \\ { H _ { 3 4 } ^ { I } } & { = \nu _ { F } ( k _ { x } + i k _ { y } ) } \\ { H _ { 4 5 } ^ { I } } & { = 2 \gamma _ { 2 0 } \cos ( 2 \pi k _ { x } / c ) } \end{array}
$$

where $c$ is the lattice constant in the out of plane direction, equal to twice the interlayer spacing. The matrix elements of $H ^ { K ^ { \prime } }$ can be obtained by replacing $k _ { x }$ by $- k _ { x }$ (other conventions for the unit cell and the orientation of the lattice lead to different phases). Recent ARPES experiments (Bostwick et al., 2007b; Ohta et al., 2006; Zhou et al., 2006a,c) performed in epitaxially grown graphene stacks (Berger et al., 2004) confirm the main features of this model, formulated mainly on the basis of Fermi surface measurements (McClure, 1957; Soule et al., 1964). The electronic spectrum of the model can also be calculated in a magnetic field (de Gennes, 1964; Nakao, 1976), and the results are also consistent with STM on graphite surfaces (Kobayashi et al., 2005; Li and Andrei, 2007; Matsui et al., 2005; Niimi et al., 2006), epitaxially grown graphene stacks (Mallet et al., 2007), and with optical measurements in the infrared range (Li et al., 2006).

# F.Surface states in graphene

So far, we have discussed the basic bulk properties of graphene. Nevertheless, graphene has very interesting surface (edge) states that do not occur in other systems. A semi-infinite graphene sheet with a zigzag edge has a band of zero energy states localized at the surface (Fujita et al., 1996; Nakada et al., 1996;

![](images/d1473e4baeb17af2fde9fc2768a5253e7a4027758f04108257e3d71ce8313176.jpg)  
Figure 14 (Color online) Ribbon geometry with zigzag edges.

Wakabayashi et al., 1999). In section II.H we will discuss the existence of edge states using the Dirac equation. Here will discuss the same problem using the tightbinding Hamiltonian. To see why these edge states exist we consider the ribbon geometry with zigzag edges shown in Fig. 14. The ribbon width is such that it has $N$ unit cells in the transverse cross section ( $y$ direction). We will assume that the ribbon has infinite length in the longitudinal direction ( $x$ direction).

Let us rewrite (5), with $t ^ { \prime } = 0$ , in terms of the integer indices $m$ and $n$ , introduced in Fig. 14, and labeling the unit cells:

$$
\begin{array} { r c l } { { \displaystyle H ~ = ~ - t \sum _ { m , n , \sigma } [ a _ { \sigma } ^ { \dagger } ( m , n ) b _ { \sigma } ( m , n ) + a _ { \sigma } ^ { \dagger } ( m , n ) b _ { \sigma } ( m - 1 , n ) } } & { { } } & { { } } \\ { { } } & { { } } & { { ~ + a _ { \sigma } ^ { \dagger } ( m , n ) b _ { \sigma } ( m , n - 1 ) + \mathrm { h . c . } ] . } } \end{array}
$$

Given that the ribbon is infinite in the $\mathbf { a } _ { 1 }$ direction one can introduce a Fourier decomposition of the operators leading to

$$
\begin{array} { r c l } { { { \cal H } ~ = ~ - t \displaystyle \int \displaystyle \frac { \mathrm { d } k } { 2 \pi } \sum _ { n , \sigma } [ a _ { \sigma } ^ { \dagger } ( k , n ) b _ { \sigma } ( k , n ) + e ^ { i k a } a _ { \sigma } ^ { \dagger } ( k , n ) b _ { \sigma } ( k , n ) } } \\ { { } } & { { } } & { { } } \\ { { } } & { { + a _ { \sigma } ^ { \dagger } ( k , n ) b _ { \sigma } ( k , n - 1 ) + \mathrm { h . c . } ] , } } \end{array}
$$

where $c _ { \sigma } ^ { \dagger } ( k , n ) \left| 0 \right. = \left| c , \sigma , k , n \right.$ , and $c = a , b$ . The oneparticle Hamiltonian can be written as:

$$
\begin{array} { r c l } { { H ^ { 1 p } ~ = ~ - t \displaystyle \int \mathrm { d } k ~ \sum _ { n , \sigma } [ ( 1 + e ^ { i k a } ) \vert a , k , n , \sigma \rangle \langle b , k , n , \sigma \vert } } \\ { { ~ } } & { { ~ } } & { { ~ + \vert a , k , n , \sigma \rangle \langle b , k , n - 1 , \sigma \vert + \mathrm { h . c . } ] . } } \end{array}
$$

The solution of the Schrödinger equation, $H ^ { 1 p } \left| \mu , k , \sigma \right. =$ $E _ { \mu , k } \left| \mu , k , \sigma \right.$ , can be generally expressed as:

$$
| \mu , k , \sigma \rangle = \sum _ { n } [ \alpha ( k , n ) | a , k , n , \sigma \rangle + \beta ( k , n ) | b , k , n , \sigma \rangle ] ,
$$

where the coefficients $\alpha$ and $\beta$ satisfy the following equations:

$$
\begin{array} { r c l } { { } } & { { } } & { { E _ { \mu , k } \alpha ( k , n ) ~ = ~ - t [ ( 1 + e ^ { i k a } ) \beta ( k , n ) + \beta ( k , n - 1 ) ] , ( } } \\ { { } } & { { } } & { { E _ { \mu , k } \beta ( k , n ) ~ = ~ - t [ ( 1 + e ^ { - i k a } ) \alpha ( k , n ) + \alpha ( k , n + 1 ) [ } } \end{array}
$$

As the ribbon has a finite width we have to be careful with the boundary conditions. Since the ribbon only exists between $n = 0$ and $n \ = \ N \mathrm { ~ - ~ } 1$ at the boundary Eqs. (48) and (49) read:

$$
\begin{array} { r c l } { { } } & { { } } & { { E _ { \mu , k } \alpha ( k , 0 ) ~ = ~ - t ( 1 + e ^ { i k a } ) \beta ( k , 0 ) , } } \\ { { } } & { { } } & { { E _ { \mu , k } \beta ( k , N - 1 ) ~ = ~ - t ( 1 + e ^ { - i k a } ) \alpha ( k , N - 1 ) . } } \end{array}
$$

The surface (edge) states are solutions of Eqs. (48-51) with $E _ { \mu , k } = 0$ .

$$
\begin{array} { l } { { 0 = ( 1 + e ^ { i k a } ) \beta ( k , n ) + \beta ( k , n - 1 ) , } } \\ { { 0 = ( 1 + e ^ { - i k a } ) \alpha ( k , n ) + \alpha ( k , n + 1 ) , } } \\ { { 0 = \beta ( k , 0 ) , } } \\ { { 0 = \alpha ( k , N - 1 ) . } } \end{array}
$$

Equations (52) and (55) are easily solved giving:

$$
\begin{array} { l } { { \alpha ( k , n ) ~ = ~ [ - 2 \cos ( k a / 2 ) ] ^ { n } e ^ { i \frac { k a } { 2 } n } \alpha ( k , 0 ) , \nonumber } } \\ { { \beta ( k , n ) ~ = ~ [ - 2 \cos ( k a / 2 ) ] ^ { N - 1 - n } e ^ { - i \frac { k a } { 2 } ( N - 1 - n ) } \beta ( k , N - ( \sharp ) 7 ) } } \end{array}
$$

Let us consider, for simplicity, a semi-infinite system with a single edge. We must require the convergence condition $| - 2 \cos ( k a / 2 ) | < 1$ , in (57) because otherwise the wavefunction would diverge in the semi-infinite graphene sheet. Therefore, the semi-infinite system has edge states for $k a$ in the region $2 \pi / 3 < k a < 4 \pi / 3$ , which corresponds to $1 / 3$ of the possible momenta. Note that the amplitudes of the edge states are given by,

$$
\begin{array} { l } { \displaystyle | \alpha ( k , n ) | ~ = ~ \sqrt { \frac { 2 } { \lambda ( k ) } } e ^ { - n / \lambda ( k ) } , } \\ { \displaystyle | \beta ( k , n ) | ~ = ~ \sqrt { \frac { 2 } { \lambda ( k ) } } e ^ { - ( N - 1 - n ) / \lambda ( k ) } , } \end{array}
$$

where the penetration length is given by:

$$
\lambda ( k ) = - 1 / \ln | 2 \cos ( k a / 2 ) | .
$$

It is easily seen that the penetration length diverges when $k a$ approaches the limits of the region $] 2 \pi / 3 , 4 \pi / 3 [$ .

Although the boundary conditions defined by Eqs. (54) and (55) are satisfied for solutions (56) and (57) in the semi-infinite system, they are not in the ribbon geometry. In fact, Eqs. (58) and (59) are eigenstates only in the semi-infinite system. In the graphene ribbon the two edge states, which come from both sides of the edge, will overlap with each other. The bonding and anti-bonding states formed by the two edge states will then be the ribbon eigenstates (Wakabayashi et al., 1999) (note that at zero energy there are no other states with which the edge states could hybridize). As bonding and anti-bonding states result in a gap in energy the zero energy fat bands of edge states will become slightly dispersive, depending on the ribbon width $N$ . The overlap between the two edge states is larger as $k a$ approaches $2 \pi / 3$ and $4 \pi / 3$ .

![](images/645f1f20b6ea632363e98aa06a13b2387468a07ebad798ff47289c9df747c897.jpg)  
Figure 15 (Color online) Sketch of a zigzag termination of a graphene bilayer. As discussed in (Castro et al., 2007b), there is a band of surface states completely localized in the bottom layer, and another surface band which alternates between the two.

This means that deviations from zero energy flatness will be stronger near these points.

Edge states in graphene nanoribbons, just as the case of carbon nanotubes, are predicted to be Luttinger liquids, that is, interacting one-dimensional electron systems (Castro Neto et al., 2006b). Hence, clean nanoribbons must have 1D square root singularities in their density of states (Nakada et al., 1996) that can be probed by Raman spectroscopy. Disorder may smooth out these singularities, however. In the presence of a magnetic field, when the bulk states are gapped, the edge states are responsible for the transport of spin and charge (Abanin et al., 2006, 2007a; Abanin and Levitov, 2007; Abanin et al., 2007b).

# G. Surface states in graphene stacks

Single layer graphene can be considered a zero gap semiconductor, which leads to the extensively studied possibility of gap states, at $\epsilon = 0$ , as discussed in the previous section. The most studied such states are those localized near a graphene zigzag edge (Fujita et al., 1996; Wakayabashi and Sigrist, 2000). It can be shown analytically (Castro et al., 2007b) that a bilayer zigzag edge, like that shown in Fig. 15, analyzed within the nearest neighbor tight-binding approximation described before, has two bands of localized states, one completely localized in the top layer and indistinguishable from similar states in single layer graphene, and another band which alternates between the two layers. These states, as they lie at $\epsilon = 0$ , have finite amplitudes on one half of the sites only.

These bands, as in single layer graphene, occupy one third of the BZ of a stripe bounded by zigzag edges. They become dispersive in a biased bilayer. As graphite can be described in terms of effective bilayer systems, one for each value of the perpendicular momentum, $k _ { z }$ , bulk graphite with a zigzag termination should show one surface band per layer.

![](images/da00daf505b26ad269781701f37a8f706b5e51a057fea661d2ec70a7b02ed307.jpg)  
Figure 16 (Color online) A piece of a honeycomb lattice displaying both zigzag and armchair edges.

# H. The spectrum of graphene nanoribbons

The spectrum of graphene nanoribbons depend very much on the nature of their edges - zigzag or armchair (Brey and Fertig, 2006a,b; Nakada et al., 1996). In Fig. 16 we show a honeycomb lattice having zigzag edges along the $x$ direction and armchair edges along the $y$ direction. If we choose the ribbon to be infinite in the $x$ direction we produce a graphene nanoribbon with zigzag edges; conversely choosing the ribbon to be macroscopically large along the $y$ but finite in the $x$ direction we produce a graphene nanoribbon with armchair edges.

In Fig. 17 we show the fourteen energy levels, calculated in the tight-binding approximation, closest to zero energy for a nanoribbon with zigzag and armchair edges and of width $N = 2 0 0$ unit cells. We can see that they are both metallic, and that the zigzag ribbon presents a band of zero energy modes that is absent in the armchair case. This band at zero energy is the surface states living near the edge of the graphene ribbon. More detailed ab initio calculations of the spectra of graphene nanoribbons show that interaction effects can lead to electronic gaps (Son et al., 2006b) and magnetic states close to the graphene edges, independent of their nature (Son et al., 2006a; Yang et al., 2007a,b).

From the experimental point of view, however, graphene nanoribbons currently have a high degree of roughness at the edges. Such edge disorder can change significantly the properties of edge states (Areshkin and White, 2007; Gunlycke et al., 2007), leading to Anderson localization, and anomalies in the quantum Hall effect (Castro Neto et al., 2006b; Martin and Blanter, 2007) as well as Coulomb blockade effects (Sols et al., 2007). Such effects have already been observed in lithographically engineered graphene nanoribbons (Han et al., 2007; Özyilmaz et al., 2007). Furthermore, the problem of edge passivation by hydrogen or other elements is not clearly understood experimentally at this time. Passivation can be modeled in the tight-binding approach by modifications of the hopping energies (Novikov, 2007c) or via additional phases in the boundary conditions (Kane and Mele, 1997). Theoretical modeling of edge passivation indicate that those have a strong effect on the electronic properties at the edge of graphene nanoribbons (Barone et al., 2006; Hod et al., 2007).

![](images/52b5586414cf5f99de97f4b14a3811aec47eb0a36b771707a2f0874e39466a47.jpg)  
Figure 17 (Color online) Left: Energy spectrum, as calculated from the tight-binding equations, for a nanoribbon with armchair(top) and zigzag(bottom) edges. The width of the nanoribbon is $N = 2 0 0$ unit cells. Only fourteen eigenstates are depicted. Right: Zoom of the low energy states shown on the right.

In what follows we derive the spectrum for both zigzag and armchair edges directly from the Dirac equation. This was originally done both with and without a magnetic field (Brey and Fertig, 2006a,b; Nakada et al., 1996).

# 1. Zigzag nanoribbons

In the geometry of Fig. 16 the unit cell vectors are ${ \bf { a } } _ { 1 } =$ $a _ { 0 } ( 1 , 0 )$ and $a _ { 2 } = a _ { 0 } \left( 1 / 2 , \sqrt { 3 } / 2 \right)$ , which generate the unit vectors of the BZ given by $b _ { 1 } = { 4 \pi } / { ( a _ { 0 } \sqrt { 3 } ) \left( \sqrt { 3 } / 2 , - 1 / 2 \right) }$ and $\smash { b _ { 2 } ^ { } ~ = ~ 4 \pi / ( a _ { 0 } \sqrt { 3 } ) ( 0 , 1 ) }$ . From these two vectors we find two inequivalent Dirac points given by $\kappa =$ $( 4 \pi / 3 a _ { 0 } , 0 ) = ( K , 0 )$ and ${ \cal K } ^ { \prime } = ( - 4 \pi / 3 a _ { 0 } , 0 ) = ( - K , 0 )$ , with $a _ { 0 } = \sqrt { 3 } a$ . The Dirac Hamiltonian around the Dirac point $\kappa$ reads in momentum space:

$$
H _ { K } = v _ { F } \left( \begin{array} { c c } { 0 } & { p _ { x } - i p _ { y } } \\ { p _ { x } + i p _ { y } } & { 0 } \end{array} \right) ,
$$

and around the $\pmb { K } ^ { \prime }$ as:

$$
H _ { K ^ { \prime } } = v _ { F } \left( \begin{array} { c c } { 0 } & { p _ { x } + i p _ { y } } \\ { p _ { x } - i p _ { y } } & { 0 } \end{array} \right) .
$$

The wavefunction, in real space, for the sublattice $A$ is given by:

$$
\Psi _ { A } ( \pmb { r } ) = e ^ { i \pmb { K } \cdot \pmb { r } } \psi _ { A } ( \pmb { r } ) + e ^ { i \pmb { K } ^ { \prime } \cdot \pmb { r } } \psi _ { A } ^ { \prime } ( \pmb { r } ) ,
$$

and for sublattice $B$ is given by

$$
\Psi _ { B } ( \pmb { r } ) = e ^ { i \pmb { K } \cdot \pmb { r } } \psi _ { B } ( \pmb { r } ) + e ^ { i \pmb { K } ^ { \prime } \cdot \pmb { r } } \psi _ { B } ^ { \prime } ( \pmb { r } ) ,
$$

where $\psi _ { A }$ and $\psi _ { B }$ are the components of the spinor wavefunction of Hamiltonian (61) and $\psi _ { A } ^ { \prime }$ and $\psi _ { B } ^ { \prime }$ have identical meaning but relatively to (62). Let us assume that the edges of the nanoribbons are parallel to the $x -$ axis. In this case, the translational symmetry guarantees that the spinor wavefunction can be written as:

$$
\psi ( \pmb { r } ) = e ^ { i k _ { x } x } \left( \begin{array} { c } { \phi _ { A } ( y ) } \\ { \phi _ { B } ( y ) } \end{array} \right) ,
$$

and a similar equation for the spinor of Hamiltonian (62). For zigzag edges the boundary conditions at the edge of the ribbon (located at $y = 0$ and $y = L$ , where $L$ is the ribbon width) are:

$$
\Psi _ { A } ( y = L ) = 0 , \qquad \Psi _ { B } ( y = 0 ) = 0 ,
$$

leading to:

$$
\begin{array} { r c l } { { } } & { { 0 } } & { { = } } \\ { { } } & { { } } & { { e ^ { i K x } e ^ { i k _ { x } x } \phi _ { A } ( L ) + e ^ { - i K x } e ^ { i k _ { x } x } \phi _ { A } ^ { \prime } ( L ) , } } \\ { { } } & { { } } & { { 0 } } \\ { { } } & { { } } & { { = } } \end{array}
$$

The boundary conditions (67) and (68) are satisfied for any $x$ by the choice:

$$
\phi _ { A } ( L ) = \phi _ { A } ^ { \prime } ( L ) = \phi _ { B } ( 0 ) = \phi _ { B } ^ { \prime } ( 0 ) = 0 .
$$

We need now to find out the form of the envelope functions. The eigenfunction around the point $\mathbf { K }$ has the form:

$$
\begin{array} { r } { \left( \begin{array} { c c } { 0 } & { k _ { x } - \partial _ { y } } \\ { k _ { x } + \partial _ { y } } & { 0 } \end{array} \right) \left( \begin{array} { l } { \phi _ { A } ( y ) } \\ { \phi _ { B } ( y ) } \end{array} \right) = \tilde { \epsilon } \left( \begin{array} { l } { \phi _ { A } ( y ) } \\ { \phi _ { B } ( y ) } \end{array} \right) , } \end{array}
$$

with $\tilde { \epsilon } = \epsilon / v _ { F }$ and $\epsilon$ the energy eigenvalue. The eigenproblem can be written as two linear differential equations of the form:

$$
\left\{ \begin{array} { l } { ( k _ { x } - \partial _ { y } ) \phi _ { B } = \tilde { \epsilon } \phi _ { A } , } \\ { ( k _ { x } + \partial _ { y } ) \phi _ { A } = \tilde { \epsilon } \phi _ { B } . } \end{array} \right.
$$

Applying the operator $( k _ { x } + \partial _ { y } )$ to the first of Eqs. leads to:

$$
( - \partial _ { y } ^ { 2 } + k _ { x } ^ { 2 } ) \phi _ { B } = \tilde { \epsilon } ^ { 2 } \phi _ { B } ,
$$

with $\phi _ { A }$ given by:

$$
\phi _ { A } = \frac { 1 } { \widetilde { \epsilon } } ( k _ { x } - \partial _ { y } ) \phi _ { B } .
$$

The solution of (72) has the form:

$$
\phi _ { B } = A e ^ { z y } + B e ^ { - z y } ,
$$

leading to an eigenenergy $\tilde { \epsilon } ^ { 2 } = k _ { x } ^ { 2 } - z ^ { 2 }$ . The boundary conditions for a zigzag edge require that $\phi _ { A } ( y = L ) = 0$ and $\phi _ { B } ( y = 0 ) = 0$ , leading to:

$$
\left\{ \begin{array} { c } { { \phi _ { B } ( y = 0 ) = 0 \Leftrightarrow A + B = 0 , } } \\ { { \phi _ { A } ( y = L ) = 0 \Leftrightarrow ( k _ { x } - z ) A e ^ { z L } + ( k _ { x } + z ) B e ^ { - z L } = 0 } } \end{array} , \right.
$$

which leads to an eigenvalue equation of the form:

$$
e ^ { - 2 z L } = \frac { k _ { x } - z } { k _ { x } + z } .
$$

Equation (76) has real solutions for $z$ , whenever $k _ { x }$ is positive; these solutions correspond to surface waves (edge states) existing near the edge of the graphene ribbon. In section II.F we discussed these states from the point of view of the tight-binding model. In addition to real solutions for $z$ , (76) also supports complex ones, of the form $z = i k _ { n }$ , leading to:

$$
k _ { x } = { \frac { k _ { n } } { \tan ( k _ { n } L ) } } .
$$

The solutions of (77) correspond to confined modes in the graphene ribbon.

If we apply the same procedure to the Dirac equation around the Dirac point $\pmb { K } ^ { \prime }$ we obtain a different eigenvalue equation given by:

$$
e ^ { - 2 z L } = \frac { k _ { x } + z } { k _ { x } - z } .
$$

This equation supports real solutions for $z$ if $k _ { x }$ is negative. Therefore we have edge states for negative values $k _ { x }$ , with momentum around $\pmb { K } ^ { \prime }$ . As in the case of $\kappa$ , the system also supports confined modes, given by:

$$
k _ { x } = - { \frac { k _ { n } } { \tan ( k _ { n } L ) } } .
$$

One should note that the eigenvalue equations for $\pmb { K } ^ { \prime }$ are obtained from those for $\kappa$ by inversion, $k _ { x }  - k _ { x }$ .

We finally notice that the edge states for zigzag nanoribbons are dispersionless (localized in real space) when $t ^ { \prime } = 0$ . When electron-hole symmetry is broken ( $t ^ { \prime } \neq 0$ ) these states become dispersive with a Fermi velocity $v _ { e } \approx t ^ { \prime } a$ (Castro Neto et al., 2006b).

# Armchair nanoribbons

Let us now consider an armchair nanoribbon with armchair edges along the $y$ direction. The boundary conditions at the edges of the ribbon (located at $x = 0$ and

$x = L$ , where $L$ is the width of the ribbon):

$$
\Psi _ { A } ( x = 0 ) = \Psi _ { B } ( x = 0 ) = \Psi _ { A } ( x = L ) = \Psi _ { B } ( x = L ) = 0 .
$$

Translational symmetry guarantees that the spinor wavefunction of Hamiltonian (61) can be written as:

$$
\psi ( \pmb { r } ) = e ^ { i k _ { y } y } \left( \begin{array} { c } { { \phi _ { A } ( x ) } } \\ { { \phi _ { B } ( x ) } } \end{array} \right) ,
$$

and a similar equation for the spinor of the Hamiltonian (62). The boundary conditions have the form:

$$
\begin{array} { l } { { 0 = e ^ { i k _ { y } y } \phi _ { A } ( 0 ) + e ^ { i k _ { y } y } \phi _ { A } ^ { \prime } ( 0 ) , } } \\ { { 0 = e ^ { i k _ { y } y } \phi _ { B } ( 0 ) + e ^ { i k _ { y } y } \phi _ { B } ^ { \prime } ( 0 ) , } } \\ { { 0 = e ^ { i K L } e ^ { i k _ { y } y } \phi _ { A } ( L ) + e ^ { - i K L } e ^ { i k _ { y } y } \phi _ { A } ^ { \prime } ( L ) , } } \\ { { 0 = e ^ { i K L } e ^ { i k _ { y } y } \phi _ { B } ( L ) + e ^ { - i K L } e ^ { i k _ { y } y } \phi _ { B } ^ { \prime } ( L ) , } } \end{array}
$$

and are satisfied for any $y$ if:

$$
\phi _ { \mu } ( 0 ) + \phi _ { \mu } ^ { \prime } ( 0 ) = 0 ,
$$

and

$$
e ^ { i K L } \phi _ { \mu } ( L ) + e ^ { - i K L } \phi _ { \mu } ^ { \prime } ( L ) = 0 ,
$$

with $\mu = A , B$ . It is clear that these boundary conditions mix states from the two Dirac points. Now we must find the form of the envelope functions obeying the boundary conditions (86) and (87). As before, the functions $\phi _ { B }$ and $\phi _ { B } ^ { \prime }$ obey the second order differential equation (72) (with $y$ replaced by $x$ ) and the function $\phi _ { A }$ and $\phi _ { A } ^ { \prime }$ are determined from (73). The solutions of (72) have the form:

$$
\begin{array} { l } { { \phi _ { B } ~ = ~ A e ^ { i k _ { n } x } + B e ^ { - i k _ { n } x } , } } \\ { { \phi _ { B } ^ { \prime } ~ = ~ C e ^ { i k _ { n } x } + D e ^ { - i k _ { n } x } . } } \end{array}
$$

Applying the boundary conditions: (86) and (87), one obtains:

$$
\begin{array} { r l } { 0 } & { = \ A + B + C + D , } \\ { 0 } & { = \ A e ^ { i ( k _ { n } + K ) L } + D e ^ { - i ( k _ { n } + K ) L } } \\ & { + \ B e ^ { - i ( k _ { n } - K ) L } + C e ^ { i ( k _ { n } - K ) L } . } \end{array}
$$

The boundary conditions are satisfied with the choice:

$$
A = - D , B = C = 0 ,
$$

which leads to $\sin [ ( k _ { n } + K ) L ] = 0$ . Therefore the allowed values of $k _ { n }$ are given by

$$
k _ { n } = \frac { n \pi } { L } - \frac { 4 \pi } { 3 a _ { 0 } } ,
$$

and the eigenenergies are given by:

$$
\tilde { \epsilon } ^ { 2 } = k _ { y } ^ { 2 } + k _ { n } ^ { 2 } .
$$

No surface states exist in this case.

# I. Dirac fermions in a magnetic field

Let us now consider the problem of a uniform magnetic field $B$ applied perpendicular to the graphene plane 2. We use the Landau gauge: $\mathbf { A } = B ( - y , 0 )$ . Notice that the magnetic field introduces a new length scale in the problem:

$$
\ell _ { B } = \sqrt { \frac { c } { e B } } ,
$$

which is the magnetic length. The only other scale in the problem is the Fermi-Dirac velocity. Dimensional analysis shows that the only quantity with dimensions of energy we can make is $v _ { F } / \ell _ { B }$ . In fact, this determines the cyclotron frequency of the Dirac fermions:

$$
\omega _ { c } = \sqrt { 2 } \frac { v _ { F } } { \ell _ { B } }
$$

(the $\sqrt { 2 }$ factor comes from the quantization of the problem, see below). Eqs. (96) and (95) show that the cyclotron energy scales like $\sqrt { B }$ , in clear contrast with the non-relativistic problem where the cyclotron energy is linear in $B$ . This implies that the energy scale associated with the Dirac fermions is rather different from the one find in the ordinary 2D electron gas. For instance, for fields of the order $B \approx 1 0$ T the cyclotron energy in the 2D electron gas is of the order of 10 K. In contrast, for the Dirac fermion problem, for the same fields, the cyclotron energy is of the order of $1 , 0 0 0$ K, that is, two orders of magnitude bigger. This has strong implications for the observation of the quantum Hall effect at room temperature (Novoselov et al., 2007). Furthermore, for $B = 1 0$ T the Zeeman energy is relatively small, $g \mu _ { B } B \approx 5$ K, and can be disregarded.

Let us now consider the Dirac equation in more detail. Using the minimal coupling in (19) (i.e., replacing $- i \nabla$ by $- i \nabla + e \mathbf { A } / c$ )we find:

$$
v _ { F } \left[ \vec { \sigma } \cdot \left( - i \nabla + e { \bf A } / c \right) \right] \psi ( { \bf r } ) = E \psi ( { \bf r } ) ,
$$

in the Landau gauge the generic solution for the wavefunction has the form $\psi ( x , y ) = e ^ { i k x } \phi ( y )$ , and the Dirac equation reads:

$$
v _ { F } \left[ \begin{array} { c c } { { 0 } } & { { \partial _ { y } - k + B e y / c } } \\ { { - \partial _ { y } - k + B e y / c } } & { { 0 } } \end{array} \right] \phi ( y ) = E \phi ( y ) ( 9 \beta \updownarrow
$$

that can be rewritten as:

$$
\omega _ { c } \left[ \begin{array} { l l } { 0 } & { \mathcal { O } } \\ { \mathcal { O } ^ { \dagger } } & { 0 } \end{array} \right] \phi ( \xi ) = E \phi ( \xi ) ,
$$

![](images/b8330856e434578515dc7fe0283ecb64535d21b6cdc35673ac14dc4ca8143abe.jpg)  
Figure 18 (Color online) SdH oscillations observed in longitudinal resistivity $\rho _ { x x }$ of graphene as a function of the charge carrier concentration $n$ . Each peak corresponds to the population of one Landau level. Note that the sequence is not interrupted when passing through the Dirac point, between electrons and holes. The period of oscillations $\Delta n = 4 B / \Phi _ { 0 }$ , where $B$ is the applied field and $\Phi _ { 0 }$ is the flux quantum (Novoselov et al., 2005a).

or equivalently:

$$
( \mathcal { O } \sigma ^ { + } + \mathcal { O } ^ { \dagger } \sigma ^ { - } ) \phi = ( 2 E / \omega _ { c } ) \phi ,
$$

where $\sigma ^ { \pm } = \sigma _ { x } \pm i \sigma _ { y }$ , and we have defined the dimensionless length scale:

$$
\xi \ = \ \frac { y } { \ell _ { B } } - \ell _ { B } k ,
$$

and 1D harmonic oscillator operators:

$$
\begin{array} { r c l } { { } } & { { \mathcal { O } = \displaystyle \frac { 1 } { \sqrt { 2 } } \left( \partial _ { \xi } + \xi \right) , } } \\ { { } } & { { } } & { { \mathcal { O } ^ { \dagger } = \displaystyle \frac { 1 } { \sqrt { 2 } } \left( - \partial _ { \xi } + \xi \right) , } } \end{array}
$$

that obey canonical commutation relations: $[ \mathcal { O } , \mathcal { O } ^ { \dagger } ] = 1$ The number operator is simply: $N = \mathcal { O } ^ { \dagger } \mathcal { O }$ .

Firstly, we notice that (100) allows for a solution with zero energy:

$$
( \mathcal { O } \sigma ^ { + } + \mathcal { O } ^ { \dagger } \sigma ^ { - } ) \phi _ { 0 } = 0 ,
$$

and since the Hilbert space generated by $\vec { \sigma }$ is of dimension 2, and the spectrum generated by $\mathcal { O } ^ { \dagger }$ is bounded from below, we just need to ensure that:

$$
\begin{array} { r } { { \mathcal { O } } \phi _ { 0 } = 0 , } \\ { \sigma ^ { - } \phi _ { 0 } = 0 , } \end{array}
$$

in order for (103) to be fulfilled. The obvious zero mode solution is:

$$
\phi _ { 0 } ( \xi ) = \psi _ { 0 } ( \xi ) \otimes | \downarrow \rangle ,
$$

where $\mid \Downarrow \rangle$ indicates the state localized on sublattice $A$ and $| \Uparrow \rangle$ indicates the state localized on sublattice $B$ . Furthermore,

$$
\mathcal { O } \psi _ { 0 } ( \xi ) = 0 ,
$$

is the ground states of the 1D harmonic oscillator. All the solutions can now be constructed from the zero mode:

$$
\begin{array} { l } { { \phi _ { N , \pm } ( \xi ) ~ = ~ \psi _ { N - 1 } ( \xi ) \otimes \mid \Uparrow \rangle \pm \psi _ { N } ( \xi ) \otimes \mid \Downarrow \rangle } } \\ { { ~ = ~ \left( \begin{array} { l } { { \psi _ { N - 1 } ( \xi ) } } \\ { { \pm \psi _ { N } ( \xi ) } } \end{array} \right) ~ , } } \end{array}
$$

and their energy is given by (McClure, 1956):

$$
E _ { \pm } ( N ) \ = \ \pm \omega _ { c } \sqrt { N } ,
$$

where $N ~ = ~ 0 , 1 , 2 , \ldots$ is a positive integer, $\psi _ { N } ( \xi )$ is the solution of the 1D Harmonic oscillator (explicitly: $\psi _ { N } ( \xi ) = 2 ^ { - N / 2 } ( N ! ) ^ { - 1 / 2 } \exp \{ - \xi ^ { 2 } / 2 \} H _ { N } ( \xi )$ where $H _ { N } ( \xi )$ is a Hermite polynomial). The Landau levels at the opposite Dirac point, $\mathrm { K } '$ , have exactly the same spectrum and hence each Landau level is doubly degenerate. Of particular importance for the Dirac problem discussed here is the existence of a zero energy state $N = 0$ which is responsible, as we are going to show, to the anomalies observed in the quantum Hall effect. This particular Landau level structure has been observed by many different experimental probes, from Shubnikov-de Haas oscillations in single layer graphene (see Fig. 18) (Novoselov et al., 2005a; Zhang et al., 2005), to infrared spectroscopy (Jiang et al., 2007a), and to scanning tunneling spectroscopy (Li and Andrei, 2007) (STS) on a graphite surface.

# J. The anomalous integer quantum Hall effect

In the presence of disorder Landau levels get broadened and mobility edges appear (Laughlin, 1981). Notice that there will be a Landau level at zero energy that separates states with hole character ( $\mu < 0$ ) from states with electron character ( $\mu > 0$ ). The components of the resistivity and conductivity tensors are related by:

$$
\begin{array} { r } { \rho _ { x x } \ = \ \frac { \sigma _ { x x } } { \sigma _ { x x } ^ { 2 } + \sigma _ { x y } ^ { 2 } } , } \\ { \rho _ { x y } \ = \ \frac { \sigma _ { x y } } { \sigma _ { x x } ^ { 2 } + \sigma _ { x y } ^ { 2 } } , } \end{array}
$$

where $\sigma _ { x x } \ \left( \rho _ { x x } \right)$ is the longitudinal component and $\sigma _ { x y }$ $( \rho _ { x y } )$ is the Hall component of the conductivity (resistivity). When the chemical potential is inside of a region of localized states the longitudinal conductivity vanishes, $\sigma _ { x x } = 0$ , and hence: $\rho _ { x x } = 0$ , $\rho _ { x y } = 1 / \sigma _ { x y }$ . On the other hand, when the chemical potential is a region of delocalized states, when the chemical potential is crossing a Landau level, we have $\sigma _ { x x } \neq 0$ and $\sigma _ { x y }$ varies continuously (Sheng et al., 2006, 2007).

![](images/8aae65fab100f2657dd4aca38fe0b6949deb1aca0aeb929f21ea68771201e08f.jpg)  
Figure 19 (Color online) Geometry of Laughlin's thought experiment with a graphene ribbon: a magnetic field $B$ is applied normal to the surface of the ribbon, a current $I$ circles the loop, generating a Hall voltage $\mathrm { V _ { H } }$ , and a magnetic flux $\Phi$ .

The value of $\sigma _ { x y }$ in the region of localized states can be obtained from Laughlin's gauge invariance argument (Laughlin, 1981): one imagines making a graphene ribbon such as the one in Fig. 19 with a magnetic field $B$ normal through its surface and a current $I$ circling its loop. Due to the Lorentz force the magnetic field produces a Hall voltage $V _ { H }$ perpendicular to the field and current. The circulating current generates a magnetic flux $\Phi$ that threads the loop. The current is given by:

$$
I = c { \frac { \delta E } { \delta \Phi } } ,
$$

where $E$ is the total energy of the system. The localized states do not respond to changes in $\Phi$ , only the delocalized ones. When the flux is changed by a flux quantum $\delta \Phi = \Phi _ { 0 } = h c / e$ the extended states remain the same by gauge invariance. If the chemical potential is in the region of localized states, all the extended states below the chemical potential will be filled both before and after the change of flux by $\Phi _ { 0 }$ . However, during the change of flux an integer number of states enter the cylinder at one edge and leave at the opposite edge.

The question is how many occupied states are transferred between edges. Let us consider a naive and as shown further incorrect calculation in order to show the importance of the zero mode in this problem. Each Landau level contributes with one state times its degeneracy $g$ . In the case of graphene we have $g \ : = \ : 4$ since there are 2 spin states and 2 Dirac cones. Hence, we would expect that when the flux changes by one flux quantum the change in energy would be $\delta E _ { \mathrm { i n c . } } = \pm 4 N e V _ { H }$ , where $N$ is an integer. The plus sign applies to electron states (charge $+ e$ ) and the minus sign to hole states (charge $- e )$ . Hence, we would conclude that $I _ { \mathrm { i n c . } } = \pm 4 ( e ^ { 2 } / h ) V _ { H }$ and hence $\sigma _ { x y , \mathrm { i n c . } } ~ = ~ I / V _ { H } ~ = ~ \pm 4 N e ^ { 2 } / h$ , which is the naive expectation. The problem with this result is that when the chemical potential is exactly at half-filling, that is, at the Dirac point, it would predict a Hall plateau at $N = 0$ with $\sigma _ { x y , \mathrm { i n c . } } = 0$ which is not possible since there is a $N = 0$ Landau level, with extended states at this energy. The solution for this paradox is rather simple: because of the presence of the zero mode which is shared by the two Dirac points, there are exactly $2 \times ( 2 N + 1 )$ occupied states that are transferred from one edge to another. Hence, the change in energy is $\delta E = \pm 2 ( 2 N + 1 ) e V _ { H }$ for a change of flux of $\delta \Phi = h c / e$ . Therefore, the Hall conductivity is (Gusynin and Sharapov, 2005; Herbut, 2007; Peres et al., 2006c,d; Schakel, 1991):

![](images/cb5614b70410a773aafb6c5514c3c5ed50b7447b1a171e55398b3df4f82ae8d0.jpg)  
Figure 20 (Color online) Quantum Hall effect in graphene as a function of charge carrier concentration. The peak at $n = 0$ shows that in high magnetic fields there appears a Landau level at zero energy where no states exist in zero field. The field draws electronic states for this level from both conduction and valence bands. The dashed line indicate plateaus in $\sigma _ { x y }$ described by Eq. (111). Adapted from (Novoselov et al., 2005a).

$$
\sigma _ { x y } = { \frac { I } { V _ { H } } } = { \frac { c } { V _ { H } } } { \frac { \delta E } { \delta \Phi } } = \pm 2 ( 2 N + 1 ) { \frac { e ^ { 2 } } { h } } ,
$$

without any Hall plateau at $N = 0$ . This amazing result has been observed experimentally (Novoselov et al., 2005a; Zhang et al., 2005) as shown in Fig.20.

# K. Tight-binding model in a magnetic field

In the tight-binding approximation the hopping integrals are replaced by a Peierls substitution:

$e ^ { i e \int _ { \mathbf { R } } ^ { \mathbf { R ^ { \prime } } } \mathbf { A \cdot d r } } t _ { \mathbf { R } , \mathbf { R ^ { \prime } } } = e ^ { i { \frac { 2 \pi } { \Phi _ { 0 } } } \int _ { \mathbf { R } } ^ { \mathbf { R ^ { \prime } } } \mathbf { A \cdot d r } } t _ { \mathbf { R } , \mathbf { R ^ { \prime } } } ,$ (112) where $t _ { \mathbf { R } , \mathbf { R ^ { \prime } } }$ represents the hopping integral between the sites $\mathbf { R }$ and $\mathbf { R } ^ { \prime }$ , with no field present. The tight-binding Hamiltonian for a single graphene layer, in a constant magnetic field perpendicular to the plane, is conveniently written as,

$$
T = - t \sum _ { m , n , \sigma } [ e ^ { i \pi \frac { \Phi } { \Phi _ { 0 } } n \frac { 1 + \ z } { 2 } } a _ { \sigma } ^ { \dag } ( m , n ) b _ { \sigma } ( m , n ) + e ^ { - i \pi \frac { \Phi } { \Phi _ { 0 } } n } a _ { \sigma } ^ { \dag } ( m , n ) b _ { \sigma } ( m - 1 , n - ( 1 - z ) / 2 ) + e ^ { i \pi \frac { \Phi } { \Phi _ { 0 } } n \frac { z - 1 } { 2 } } a _ { \sigma } ^ { \dag } ( m , n ) ]
$$

holding for a graphene stripe with a zigzag ( $z = 1$ ) and armchair ( $z = - 1$ ) edges oriented along the $x -$ direction.

Fourier transforming along the $x$ direction gives,

$$
H = - t \sum _ { k , n , \sigma } [ e ^ { i \pi { \frac { \phi } { \Phi _ { 0 } } } n { \frac { 1 + z } { 2 } } } a _ { \sigma } ^ { \dagger } ( k , n ) b _ { \sigma } ( k , n ) + e ^ { - i \pi { \frac { \phi } { \Phi _ { 0 } } } n } e ^ { i k a } a _ { \sigma } ^ { \dagger } ( k , n ) b _ { \sigma } ( k , n - ( 1 - z ) / 2 ) + e ^ { i \pi { \frac { \phi } { \Phi _ { 0 } } } n { \frac { z - 1 } { 2 } } }
$$

Let us now consider the case of zigzag edges. The eigenproblem can be rewritten in terms of Harper's equations (Harper, 1955), and for zigzag edges we obtain

(Rammal, 1985):

$$
\begin{array} { r l r } { \displaystyle E _ { \mu , k } \alpha ( k , n ) ~ = ~ - t [ e ^ { i k a / 2 } 2 \cos ( \pi \frac \Phi { \Phi _ { 0 } } n - \frac { k a } 2 ) \beta ( k , n ) } & { } & \\ { ~ + ~ \beta ( k , n - 1 ) ] , } & { ( 1 \mathrm { ! } } \\ { \displaystyle E _ { \mu , k } \beta ( k , n ) ~ = ~ - t [ e ^ { - i k a / 2 } 2 \cos ( \pi \frac \Phi { \Phi _ { 0 } } n - \frac { k a } 2 ) \alpha ( k , n ) } & { } & \\ { ~ + ~ \alpha ( k , n + 1 ) ] , } & { ( 1 \mathrm { ! } } \end{array}
$$

![](images/f7b19c65d5ad0a4e9d54e952ba09392188cc27710a0ed88127693d2ac3073736.jpg)  
Figure 21 (Color online) Fourteen Energy levels of tightbinding electrons in graphene in the presence of a magnetic flux $\Phi = \Phi _ { 0 } / 7 0 1$ , for a finite stripe with $N = 2 0 0$ unit cells. The bottom panels are zoom in images of the top ones. The dashed line represents the chemical potential $\mu$ .

where the coefficients $\alpha ( k , n )$ and $\beta ( k , n )$ show up in Hamiltonian's eigenfunction, $| \psi ( k ) \rangle$ , written in terms of lattice-position-state states as:

$$
| \psi ( k ) \rangle = \sum _ { n , \sigma } ( \alpha ( k , n ) | a ; k , n , \sigma \rangle + \beta ( k , n ) | b ; k , n , \sigma \rangle ) ~ .
$$

Eqs. (114) and (115) hold in the bulk. Considering that the zigzag ribbon has $N$ unit cells along its width, from $n = 0$ to $n = N - 1$ , the boundary conditions at the edges are obtained from Eqs. (114) and (115), and read

$$
E _ { \mu , k } \alpha ( k , 0 ) = - t e ^ { i k a / 2 } 2 \cos \left( \frac { k a } { 2 } \right) \beta ( k , 0 ) ,
$$

$$
E _ { \mu , k } \beta ( k , N - 1 ) { = } { - } 2 t e ^ { - i { k a / 2 } } \mathrm { c o s } \left[ \pi { \frac { \Phi } { \Phi _ { 0 } } } ( N { - } 1 ) { - } { \frac { k a } { 2 } } \right] \alpha ( k , N { - } 1 ) { = } \nonumber
$$

Similar equations hold for a graphene ribbon with armchair edges.

In Fig. 21 we show fourteen energy levels, around zero energy, for both the zigzag and the armchair cases. The formation of the Landau levels is signaled by the presence of flat energy bands, following the bulk energy spectrum. From Fig. 21 it is straightforward to obtain the value of the Hall conductivity in the quantum Hall effect regime. Let us assume that the chemical potential is in between two Landau levels at positive energies, as represented by the dashed line in Fig. 21. The Landau level structure shows two zero energy modes, one of them is electron-like (hole-like), since close to the edge of the sample its energy is shifted upwards (downwards). The other Landau levels are doubly degenerate. The determination of the values for the Hall conductivity is done by counting how many energy levels (of electron-like nature) are below chemical potential. This counting produces the value $( 2 N + 1 )$ , with $N = 0 , 1 , 2 , \ldots$ (for the case of Fig. 21 one has $N \ = \ 2$ ). From this counting the Hall conductivity is given, including an extra factor of two accounting for the spin degree of freedom, by

![](images/dd9c4c9d954c80e1efa33de348259b737e33348e3ef23a3b1a6b425d25c44029.jpg)  
Figure 22 (Color online) Landau levels of the graphene stacks shown in Fig.13. The applied magnetic field is $1 \textrm { T }$ .

$$
\sigma _ { x y } = \pm 2 \frac { e ^ { 2 } } { h } ( 2 N + 1 ) = \pm 4 \frac { e ^ { 2 } } { h } \left( N + \frac { 1 } { 2 } \right) .
$$

Eq. (119) leads to the anomalous integer quantum Hall effect discussed previously, which is the hallmark of single layer graphene.

# L. Landau levels in graphene stacks

The dependence of the Landau level energies on the Landau index $N$ roughly follows the dispersion relation of the bands, as shown in Fig. 22. Note that, in a trilayer with Bernal stacking, two sets of levels have a $\sqrt { N }$ dependence, while the energies of other two depend linearly on $N$ . In an infinite rhombohedral stack, the Landau lev'els remain discrete and quasi-2D (Guinea et al., 2006). Note that, even in an infinite stack with the Bernal structure, the Landau level closest to $E = 0$ forms a band which does not overlap with the other Landau levels, leading to the possibility of a 3D integer quantum Hall effect (Bernevig et al., 2007; Kopelevich et al., 2006, 2003; Luk'yanchuk and Kopelevich, 2004).

The optical transitions between Landau levels can also be calculated. The selection rules are the same as for a graphene single layer, and only transitions between subbands with Landau level indices $M$ and $N$ such that $\vert N \vert = \vert M \pm 1 \vert$ are allowed. The resulting transitions, with their respective spectral strengths, are shown in Fig. 23. The transitions are grouped into subbands, which give rise to a continuum when the number of layers tends to infinity. In Bernal stacks with an odd number of layers, the transitions associated to Dirac subbands with linear dispersion have the largest spectral strength, and they give a significant contribution to the total absorption even if the number of layers is large, $N _ { \mathrm { L } } \lesssim 3 0 - 4 0$ (Sadowski et al., 2006).

![](images/0ff7abc352f376720ed7e0e1186d519efe8929693ba9ef37dee6058491d78e81.jpg)  
Figure 23 (Color online) Relative spectral strength of the low energy optical transitions between Landau levels in graphene stacks with Bernal ordering and an odd number of layers. The applied magnetic field is 1 T. Top left: 3 layers. Top right: 11 layers. Bottom: 51 layers. The large red circles are the transitions in a single layer.

# M. Diamagnetism

Back in 1952 Mrozowski (Mrozowski, 1952) studied diamagnetism of polycrystalline graphite and other condensed-matter molecular-ring systems. It was concluded that in such ring systems diamagnetism has two contributions: (1) diamagnetism from the filled bands of electrons; (2) Landau diamagnetism of free electrons and holes. For graphite the second source of diamagnetism is by far the largest of the two.

McClure (McClure, 1956) computed diamagnetism of a 2D honeycomb lattice using the theory introduced by Wallace (Wallace, 1947), a calculation he later generalized to three dimensional graphite (McClure, 1960). For the honeycomb plane the magnetic susceptibility in units of emu/g is

$$
\chi = - \frac { 0 . 0 0 1 4 } { T } \gamma _ { 0 } ^ { 2 } \mathrm { s e c h } ^ { 2 } \left( \frac { \mu } { 2 k _ { B } T } \right) ,
$$

where $\mu$ is the Fermi energy, $T$ the temperature, and $k _ { B }$ the Boltzmann constant. For graphite the magnetic susceptibility is anisotropic and the difference between the susceptibility parallel to the principal axis and that perpendicular to the principal axis is $- 2 1 . 5 { \times } 1 0 ^ { - 6 }$ emu/g. The susceptibility perpendicular to the principal axis is about the free-atom susceptibility of $\cdot 0 . 5 \times 1 0 ^ { - 6 }$ emu/g. In the 2D model the susceptibility turns out to be large due to the presence of fast moving electrons, with a velocity of the order of $v _ { F } \simeq 1 0 ^ { 6 } ~ \mathrm { m / s }$ , which in turn is a consequence of the large value of the hopping parameter $\gamma _ { 0 }$ . In fact the susceptibility turns out to be proportional to the square of $\gamma _ { 0 }$ . Later Sharma et al. extended the calculation of McClure for graphite by including the effect of trigonal warping and showed that the low temperature diamagnetism increases (Sharma et al., 1974).

Safran and DiSalvo (Safran and DiSalvo, 1979), interested in the susceptibility of graphite intercalation compounds, recalculated, in a tight-binding model, the susceptibility perpendicular to a graphite plane using Fukuyama's theory (Fukuyama, 1971), which includes interband transitions. The results agree with those published by McClure (McClure, 1956). Later, Safran computed the susceptibility of a graphene bilayer showing that this system is diamagnetic at small values of the Fermi energy, but there appears a paramagnetic peak when the Fermi energy is of the order of the interlayer coupling (Safran, 1984).

The magnetic susceptibility of other carbon based materials, as carbon nanotubes and $C _ { 6 0 }$ molecular solids was measured (Heremans et al., 1994) showing a diamagnetic response at finite magnetic fields different from that of graphite. The study of the magnetic response of nanographite ribbons with both zig-zag and arm-chair edges was done by Wakabayashi et al. using a numerical differentiation of the free energy (Wakabayashi et al., 1999). From these two systems, the zig-zag edge state is of particular interest since the system supports edge states even in the presence of a magnetic field, leading to very high density of states near the edge of the ribbon. The high temperature response of these nanoribbons was found to be diamagnetic whereas the low temperature susceptibility is paramagnetic.

The Dirac-like nature of the electronic quasiparticles in graphene led (Ghosal et al., 2007) to consider in general the problem of the diamagnetism of nodal fermions and Nakamura to study the orbital magnetism of Dirac fermions in weak magnetic fields(Nakamura, 2007). Koshino and Ando considered the diamagnetism of disordered graphene in the self consistent Born approximation, with a disorder potential of the form $U ( { \pmb r } ) =$ ${ \bf 1 } u _ { i } \delta ( { \pmb r } - { \pmb R } )$ (Koshino and Ando, 2007). At the neutrality point and zero temperature the susceptibility of disordered graphene is given by

$$
\chi ( 0 ) = - \frac { g _ { v } g _ { s } } { 3 \pi ^ { 2 } } e ^ { 2 } \gamma _ { 0 } ^ { 2 } \frac { 2 W } { \Gamma _ { 0 } } ,
$$

where $g _ { s } = g _ { v } = 2$ is the spin and valley degeneracies, $W$ is a dimensionless parameter for the disorder strength, defined as $W = n _ { i } u _ { i } ^ { 2 } / 4 \pi \gamma _ { 0 } ^ { 2 }$ , $n _ { i }$ the impurity concentration, and $\Gamma _ { 0 }$ is given by $\Gamma _ { 0 } = \epsilon _ { c } \exp [ - 1 / ( 2 W ) ]$ with $\epsilon _ { c }$ a parameter defining a cut-off function used in the theory. At finite Fermi energy $\epsilon _ { F }$ and zero temperature the magnetic susceptibility is given by

$$
\chi ( \epsilon _ { F } ) = - { \frac { g _ { v } g _ { s } } { 3 \pi } } e ^ { 2 } \gamma _ { 0 } ^ { 2 } { \frac { W } { \vert \epsilon _ { F } \vert } } .
$$

In the clean limit the susceptibility is given by (Koshino and Ando, 2007; McClure, 1956; Safran and DiSalvo, 1979):

$$
\chi ( \epsilon _ { F } ) = - { \frac { g _ { v } g _ { s } } { 6 \pi } } e ^ { 2 } \gamma _ { 0 } ^ { 2 } \delta ( \epsilon _ { F } ) .
$$

Spin-orbit coupling describes a process in which an electron changes simultaneously its spin and angular momentum or, in general, moves from one orbital wavefunction to another. The mixing of the spin and the orbital motion is a relativistic effect, which can be derived from Dirac's model of the electron. It is large in heavy ions, where the average velocity of the electrons is higher. Carbon is a light atom, and the spin orbit interaction is expected to be weak.

The symmetries of the spin orbit interaction, however, allow the formation of a gap at the Dirac points in clean graphene. The spin orbit interaction leads to a spin dependent shift of the orbitals, which is of a different sign for the two sublattices, acting as an effective mass within each Dirac point (Dresselhaus and Dresselhaus, 1965; Kane and Mele, 2005; Wang and Chakraborty, 2007a). The appearance of this gap leads to a non trivial spin Hall conductance, in similar way to other models which study the parity anomaly in relativistic field theory in (2+1) dimensions (Haldane, 1988). When the inversion symmetry of the honeycomb lattice is broken, either because the graphene layer is curved or because an external electric field is applied (Rashba interaction) additional terms, which modulate the nearest neighbor hopping, are induced (Ando, 2000). The intrinsic and extrinsic spin orbit interactions can be written as (Dresselhaus and Dresselhaus, 1965; Kane and Mele, 2005):

$$
\begin{array} { l } { { \displaystyle \mathcal { H } _ { S O ; i n t } \equiv \Delta _ { s o } \int d ^ { 2 } { \bf r } \hat { \Psi } ^ { \dagger } ( { \bf r } ) \hat { s } _ { z } \hat { \sigma } _ { z } \hat { \tau } _ { z } \hat { \Psi } ( { \bf r } ) } , } \\ { { \displaystyle \mathcal { H } _ { S O ; e x t } \equiv \lambda _ { R } \int d ^ { 2 } { \bf r } \hat { \Psi } ^ { \dagger } ( { \bf r } ) ( - \hat { s } _ { x } \hat { \sigma } _ { y } + \hat { s } _ { y } \hat { \sigma } _ { x } \hat { \tau } _ { z } ) \hat { \Psi } ( { \bf r } ) , } } \end{array}
$$

where $\hat { \sigma }$ and $\hat { \tau }$ are Pauli matrices which describe the sublattice and valley degrees of freedom, and $\hat { s }$ are Pauli matrices acting on actual spin space. $\Delta _ { s o }$ is the spin-orbit coupling and $\lambda _ { R }$ is the Rashba coupling.

The magnitude of the spin orbit coupling in graphene can be inferred from the known spin orbit coupling in the carbon atom. This coupling allows for transitions between the $p _ { z }$ and $p _ { x }$ and $p _ { y }$ orbitals. An electric field induces also transitions between the $p _ { z }$ and $s$ orbitals. These intra-atomic processes mix the $\pi$ and $\sigma$ bands in graphene. Using second order perturbation theory, one obtains a coupling between the low energy states in the $\pi$ band. Tight-binding (Huertas-Hernando et al., 2006; Zarea and Sandler, 2007) and band structure calculations (Min et al., 2006; Yao et al., 2007) give estimates for the intrinsic and extrinsic spin-orbit interactions in the range $0 . 0 1 - 0 . 2 $ K, and hence much smaller than the other energy scales in the problem (kinetic, interaction, and disorder).

# III. FLEXURAL PHONONS, ELASTICITY, AND CRUMPLING

Graphite, in the Bernal stacking configuration, is a layered crystalline solid with 4 atoms per unit cell. Its basic structure is essentially a repetition of the bilayer structure discussed earlier. The coupling between the layers, as we discussed, is weak and, therefore, graphene has been always the starting point for the discussion of phonons in graphite (Wirtz and Rubio, 2004). Graphene has two atoms per unit cell and if we consider graphene as strictly 2D it should have 2 acoustic modes (with dispersion $\omega _ { \mathrm { a c } } ( k ) \propto k$ as $k  0$ ) and 2 optical modes (with dispersion $\omega _ { \mathrm { o p } } ( k ) \propto$ constant, as $k  0$ ) solely due to the inplane translation and stretching of the graphene lattice. Nevertheless, graphene exists in the 3D space and hence the atoms can oscillate out-of-plane leading to 2 out-ofplane phonons (one acoustic and another optical) called flexural modes. The acoustic flexural mode has dispersion $\omega _ { \mathrm { { f l e x } } } ( k ) \propto k ^ { 2 }$ as $k  0$ which represents the translation of the whole graphene plane (essentially a one atom thick membrane) in the perpendicular direction (free particle motion). The optical flexural mode represents the out-of-phase out-of-plane oscillation of the neighboring atoms. In first approximation, if we neglect the coupling between graphene planes, graphite has essentially the same phonon modes, albeit they are degenerate. The coupling between planes has two main effects: (1) it lifts the degeneracy of the phonon modes, and (2) leads to a strong suppression of the energy of the flexural modes. Theoretically, flexural modes should become ordinary acoustic and optical modes in a fully covalent 3D solid, but in practice, the flexural modes survive due to the fact the planes are coupled by weak van der Waals-like forces. These modes have been measured experimentally in graphite (Wirtz and Rubio, 2004). Graphene can also be obtained as a suspended membrane, that is, without a substrate, being supported only by a scaffold or bridging micron-scale gaps (Bunch et al., 2007; Meyer et al., 2007a,b). Figure 24 shows a suspended graphene sheet and an atomic resolution image of its crystal lattice.

Because the flexural modes disperse like $k ^ { 2 }$ they dominate the behavior of structural fluctuations in graphene at low energies (low temperatures) and long wavelengths. It is instructive to understand how these modes appear from the point of view of elasticity theory (Chaikin and Lubensky, 1995; Nelson et al., 2004). Consider, for instance, a graphene sheet in 3D and let us parameterize the position of the sheet relative of a fixed coordinate frame in terms of the in-plane vector r and the height variable $h ( \mathbf { r } )$ so that a position in the graphene is given by the vector $\mathbf { R } = ( \mathbf { r } , h ( \mathbf { r } ) )$ . The unit vector normal to the surface is given by:

$$
{ \bf N } = \frac { { \bf z } - \nabla h } { \sqrt { 1 + ( \nabla h ) ^ { 2 } } } ,
$$

where $\nabla \ : = \ : ( \partial _ { x } , \partial _ { y } )$ is the 2D gradient operator, and $\mathbf { z }$ is the unit vector in the third direction. In a flat graphene configuration all the normal vectors are aligned and therefore $\nabla \cdot \mathbf { N } = 0$ . Deviations from the flat configuration requires misalignment of the normal vectors and costs elastic energy. This elastic energy can be written as:

![](images/4e27ba0d83b435a148a5b942e479464b879a41931ac0608d05485f183069eb69.jpg)  
Figure 24 (Color online) Suspended graphene sheet. (a) Bright-field transmission-electron-microscope image of a graphene membrane. Its central part (homogeneous and featureless region) is monolayer graphene. Adapted from (Meyer et al., 2007a). (b) Despite being only one atom thick, graphene remains a perfect crystal as this atomic resolution image shows. The image is obtained in a scanning transmission electron microscope. The visible periodicity is given by the lattice of benzene rings. Adapted from (Nair et al., 2008).

$$
E _ { 0 } = { \frac { \kappa } { 2 } } \int d ^ { 2 } \mathbf { r } \left( \nabla \cdot \mathbf { N } \right) ^ { 2 } \approx { \frac { \kappa } { 2 } } \int d ^ { 2 } \mathbf { r } \left( \nabla ^ { 2 } h \right) ^ { 2 }
$$

where $\kappa$ is the bending stiffness of graphene, and the expression in terms of $h ( \mathbf { r } )$ is valid for smooth distortions of the graphene sheet $( ( \nabla h ) ^ { 2 } \ll 1$ ). The energy (126) is valid in the absence of a surface tension or a substrate which break the rotational and translational symmetry of the problem, respectively. In the presence of tension there is an energy cost for solid rotations of the graphene sheet ( $\nabla h \neq 0$ ) and hence a new term has to be added

to the elastic energy:

$$
E _ { \mathrm { T } } = \frac { \gamma } { 2 } \int d ^ { 2 } { \bf r } \left( \nabla h \right) ^ { 2 } ,
$$

where $\gamma$ is the interfacial stiffness. A substrate, described by a height variable $s ( \mathbf { r } )$ , pins the graphene sheet through van der Waals and other electrostatic potentials so that the graphene configuration tries to follow the substrate $h ( \mathbf { r } ) \sim s ( \mathbf { r } )$ . Deviations from this configuration cost extra elastic energy that can be approximated by a harmonic potential (Swain and Andelman, 1999):

$$
E _ { \mathrm { S } } = \frac { v } { 2 } \int d ^ { 2 } { \bf r } \left( s ( { \bf r } ) - h ( { \bf r } ) \right) ^ { 2 } ,
$$

where $\boldsymbol { v }$ characterizes the strength of the interaction potential between substrate and graphene.

Firstly, let us consider the free floating graphene problem (126) that we can rewrite in momentum space as:

$$
E _ { 0 } = \frac { \kappa } { 2 } \sum _ { { \bf k } } k ^ { 4 } h _ { - { \bf k } } h _ { { \bf k } } .
$$

We now canonically quantize the problem by introducing a momentum operator $P _ { \mathbf { k } }$ that has the following commutator with $h _ { \mathbf { k } }$ :

$$
\left[ h _ { \mathbf { k } } , P _ { \mathbf { k ^ { \prime } } } \right] = i \delta _ { \mathbf { k } , \mathbf { k ^ { \prime } } } ,
$$

and write the Hamiltonian as:

$$
H = \sum _ { { \bf k } } \left\{ \frac { P _ { - { \bf k } } P _ { \bf k } } { 2 \sigma } + \frac { \kappa k ^ { 4 } } { 2 } h _ { - { \bf k } } h _ { \bf k } \right\} ,
$$

where $\sigma$ is graphene's 2D mass density. From the Heisenberg equations of motion for the operators it is trivial to find that $h _ { \mathbf { k } }$ oscillates harmonically with a frequency given by:

$$
\omega _ { \mathrm { f l e x } } ( { \bf k } ) = \left( \frac { \kappa } { \sigma } \right) ^ { 1 / 2 } k ^ { 2 } ,
$$

which is the long wavelength dispersion of the flexural modes. In the presence of tension it is easy to see that the dispersion is modified to:

$$
\omega ( \mathbf { k } ) = k \sqrt { \frac { \kappa } { \sigma } k ^ { 2 } + \frac { \gamma } { \sigma } } ,
$$

indicating that the dispersion of the flexural modes becomes linear in $k$ , as $k  0$ , under tension. That is what happens in graphite where the interaction between layers breaks the rotational symmetry of the graphene layers.

Eq. (132) also allows us to relate the bending energy of graphene with the Young modulus, $Y$ , of graphite. The fundamental resonance frequency of a macroscopic graphite sample of thickness $t$ is given by (Bunch et al., 2007):

$$
\nu ( { \bf k } ) = \left( \frac { Y } { \rho } \right) ^ { 1 / 2 } t k ^ { 2 } ,
$$

where $\rho = \sigma / t$ is the 3D mass density. Assuming that (134) works down to the single plane level, that is, when $t$ is the distance between planes, we find:

$$
\kappa = Y t ^ { 3 } ,
$$

which provides a simple relationship between the bending stiffness and the Young modulus. Given that $Y \approx 1 0 ^ { 1 2 }$ N/m and $t \approx 3 . 4$ $\mathrm { \AA }$ we find, $\kappa \approx 1$ eV. This result is in good agreement with ab initio calculations of the bending rigidity (Lenosky et al., 1992; Tu and Ou-Yang, 2002) and experiments in graphene resonators (Bunch et al., 2007).

The problem of structural order of a "free floating" graphene sheet can be fully understood from the existence of the flexural modes. Consider, for instance, the number of flexural modes per unit of area at a certain temperature $T$ :

$$
N _ { \mathrm { p h } } = \int \frac { d ^ { 2 } { \bf k } } { ( 2 \pi ) ^ { 2 } } n _ { \bf k } = \frac { 1 } { 2 \pi } \int _ { 0 } ^ { \infty } d k \frac { k } { e ^ { \beta \sqrt { \kappa / \sigma } k ^ { 2 } } - 1 }
$$

where $n _ { \mathbf { k } }$ is the Bose-Einstein occupation number ( $\beta =$ $1 / ( k _ { B } T ) )$ . For $T \neq 0$ the above integral is logarithmically divergent in the infrared ( $k  0$ ) indicating a divergent number of phonons in the thermodynamic limit. For a system with finite size $L$ the smallest possible wave vector is of the order of $k _ { \mathrm { m i n } } \sim 2 \pi / L$ . Using $k _ { \mathrm { m i n } }$ as a lower cutoff in the integral (136) we find:

$$
N _ { \mathrm { p h } } = \frac { \pi } { L _ { \mathrm { T } } ^ { 2 } } \ln \left( \frac { 1 } { 1 - e ^ { - L _ { \mathrm { T } } ^ { 2 } / L ^ { 2 } } } \right) ,
$$

where

$$
L _ { \mathrm { T } } = \frac { 2 \pi } { \sqrt { k _ { B } T } } \left( \frac { \kappa } { \sigma } \right) ^ { 1 / 4 } ,
$$

is the thermal wavelength of the flexural modes. Notice that that when $L \gg L _ { \mathrm { T } }$ the number of flexural phonons in (137) diverges logarithmically with the size of the system:

$$
N _ { \mathrm { p h } } \approx \frac { 2 \pi } { L _ { \mathrm { T } } ^ { 2 } } \ln \left( \frac { L } { L _ { \mathrm { T } } } \right) ,
$$

indicating that the system cannot be structurally ordered at any finite temperature. This is nothing but the crumpling instability of soft membranes (Chaikin and Lubensky, 1995; Nelson et al., 2004). For $L \ll L _ { \mathrm { T } }$ one finds that $N _ { \mathrm { p h } }$ goes to zero exponentially with the size of the system indicating that systems with finite size can be flat at sufficiently low temperatures. Notice that for $\kappa \approx 1 \ \mathrm { e V }$ , $\rho \approx 2 2 0 0 ~ \mathrm { k g / m ^ { 3 } }$ b $t \ : = \ : 3 . 4$ $\mathring { \mathrm { A } } \left( \sigma \ \approx \ 7 . 5 \times 1 0 ^ { - 7 } \ \mathrm { k g / m ^ { 2 } } \right)$ , and $T \approx 3 0 0 \mathrm { ~ K ~ }$ , we find $L _ { \mathrm { T } } \approx 1 \mathrm { ~ \AA ~ }$ indicating that free floating graphene should always crumple at room temperature due to thermal fluctuations associated with flexural phonons. Nevertheless, the previous discussion only involves the harmonic (quadratic part) of the problem. Non-linear effects such as large bending deformations (Peliti and Leibler, 1985), the coupling between flexural and in-plane modes (or phonon-phonon interactions (Bonini et al., 2007; Radzihovsky and Le Doussal, 1992)) and the presence of topological defects (Nelson and Peliti, 1987) can lead to strong renormalizations of the bending rigidity, driving the system toward a flat phase at low temperatures (Chaikin and Lubensky, 1995). This situation has been confirmed in numerical simulations of free graphene sheets (Adebpour et al., 2007; Fasolino et al., 2007).

The situation is rather different if the system is under tension or in the presence of a substrate or scaffold that can hold the graphene sheet. In fact, static rippling of graphene flakes suspended on scaffolds have been observed for single layer as well as bilayers (Meyer et al., 2007a,b). In this case the dispersion, in accordance with (133), is at least linear in $k$ , and the integral in (136) converges in the infrared ( $k \  \ 0$ ) indicating that the number of flexural phonons is finite and graphene does not crumple. We should notice that these thermal fluctuations are dynamic and hence average to zero over time, therefore, the graphene sheet is expected to be flat under these circumstances. Obviously, in the presence of a substrate or scaffold described by (128) static deformations of the graphene sheet are allowed. Also, hydrocarbon molecules that are often present on top of free hanging graphene membranes could quench flexural fluctuations making them static.

Finally, one should notice that in the presence of a metallic gate the electron-electron interactions lead to the coupling of the phonon modes to the electronic excitations in the gate. This coupling could be partially responsible to the damping of the phonon modes due to dissipative effects (Seoanez et al., 2007) as observed in graphene resonators (Bunch et al., 2007).

# IV. DISORDER IN GRAPHENE

Graphene is a remarkable material from the electronic point of view. Because of the robustness and specificity of the sigma bonding, it is very hard for alien atoms to replace the carbon atoms in the honeycomb lattice. This is one of the reasons why the electron mean free path in graphene can be so long, reaching up to one micrometer in the existing samples. Nevertheless, graphene is not immune to disorder and its electronic properties are controlled by extrinsic as well as intrinsic effects that are unique to this system. Among the intrinsic sources of disorder we can highlight: surface ripples and topological defects. Extrinsic disorder can come about in many different forms: adatoms, vacancies, charges on top of graphene or in the substrate, and extended defects such as cracks and edges.

It is easy to see that from the point of view of single electron physics (that is, terms that can be added to (5)), there are two main terms that disorder couples to. The first one is a local change in the single site energy,

$$
H _ { \mathrm { d d } } = \sum _ { i } V _ { i } \left( a _ { i } ^ { \dagger } a _ { i } + b _ { i } ^ { \dagger } b _ { i } \right) ,
$$

where $V _ { i }$ is the strength of the disorder potential on site $\mathbf { R } _ { i }$ , which is diagonal in the sublattice indices and hence, from the point of view of the Dirac Hamiltonian (18), can be written as:

$$
H _ { \mathrm { d d } } = \int d ^ { 2 } r \sum _ { a = 1 , 2 } V _ { a } ( { \bf r } ) \hat { \Psi } _ { a } ^ { \dagger } ( { \bf r } ) \hat { \Psi } _ { a } ( { \bf r } ) ,
$$

which acts as a chemical potential shift for the Dirac fermions, that is, shifts locally the Dirac point.

Because of the vanishing of the density of states in single layer graphene, and by consequence the lack of electrostatic screening, charge potentials may be rather important in determining the spectroscopic and transport properties (Adam et al., 2007; Ando, 2006b; Nomura and MacDonald, 2007). Of particular importance is the Coulomb impurity problem where,

$$
V _ { a } ( r ) = \frac { e ^ { 2 } } { \epsilon _ { 0 } } \frac { 1 } { r } ,
$$

where $\epsilon _ { 0 }$ is the dielectric constant of the medium. The solution of the Dirac equation for the Coulomb potential in 2D can be studied analytically (Biswas et al., 2007; DiVincenzo and Mele, 1984; Novikov, 2007a; Pereira et al., 2007b; Shytov et al., 2007). Its solution has many of the features of the 3D relativistic hydrogen atom problem (Baym, 1969). Just as in the case of the 3D problem the nature of the eigenfunctions depends strongly on graphene's dimensionless coupling constant:

$$
g = \frac { Z e ^ { 2 } } { \epsilon _ { 0 } v _ { F } } .
$$

Notice, therefore, that the coupling constant can be varied by either changing the charge of the impurity, $Z$ , or modifying the dielectric environment and changing $\epsilon _ { \mathrm { 0 } }$ . For $g < g _ { c } = 1 / 2$ the solutions of this problem are given in terms of Coulomb wavefunctions with logarithmic phase shifts. The local density of states (LDOS) is affected close to the impurity due the electron-hole asymmetry generated by the Coulomb potential. The local charge density decays like $1 / r ^ { 3 }$ plus fast oscillations of the order of the lattice spacing (in the continuum limit this would give rise to a Dirac delta function for the density (Kolezhuk et al., 2006)). Just like in 3D QED, the 2D problem becomes unstable for $g > g _ { c } = 1 / 2$ leading to super-critical behavior and the so-called fall of electron to the center (Landau and Lifshitz, 1981). In this case the LDOS is strongly affected by the presence of the Coulomb impurity with the appearance of bound states outside the band and scattering resonances within the band (Pereira et al., 2007b) and the local electronic density decays monotonically like $1 / r ^ { 2 }$ at large distances.

It has been argued (Schedin et al., 2007) that without high vacuum environment these Coulomb effects can be strongly suppressed by large effective dielectric constants due to the presence of a nanometer thin layer of absorbed water (Sabio et al., 2007). In fact, experiments in ultra-high vacuum conditions (Chen et al., 2007b) display strong scattering features in the transport that can be associated to charge impurities. Screening effects that affect the strength and range of the Coulomb interaction, are rather non-trivial in graphene (Fogler et al., 2007b; Shklovskii, 2007) and, therefore, important for the interpretation of transport data (Bardarson et al., 2007; Lewenkopf et al., 2007; Nomura et al., 2007; San-Jose et al., 2007).

Another type of disorder is the one that changes the distance or angles between the $p _ { z }$ orbitals. In this case, the hopping energies between different sites are modified leading to a new term to the original Hamiltonian (5):

$$
\begin{array} { r c l } { { H _ { \mathrm { o d } } } } & { { = } } & { { \displaystyle \sum _ { i , j } \left\{ \delta t _ { i j } ^ { ( a b ) } \left( a _ { i } ^ { \dagger } b _ { j } + \mathrm { h . c . } \right) \right. } } \\ { { } } & { { + } } & { { \displaystyle \left. \delta t _ { i j } ^ { ( a a ) } \left( a _ { i } ^ { \dagger } a _ { j } + b _ { i } ^ { \dagger } b _ { j } \right) \right\} , } } \end{array}
$$

or in Fourier space:

$$
\begin{array} { r l } { { \displaystyle H _ { \mathrm { o d } } ~ = ~ \sum _ { { \bf k } , { \bf k } ^ { \prime } } a _ { { \bf k } } ^ { \dagger } b _ { { \bf k } ^ { \prime } } \sum _ { i , { \vec { \delta } } _ { a b } } \delta t _ { i } ^ { ( a b ) } e ^ { i ( { \bf k } - { \bf k } ^ { \prime } ) \cdot { \bf R } _ { i } - i \vec { \delta } _ { a a } \cdot { \bf k } ^ { \prime } } + \mathrm { h . c . } } } \\ { { \displaystyle ~ + \left( a _ { { \bf k } } ^ { \dagger } a _ { { \bf k } ^ { \prime } } + b _ { { \bf k } } ^ { \dagger } b _ { { \bf k } ^ { \prime } } \right) \sum _ { i , \vec { \delta } _ { a a } } \delta t _ { i } ^ { ( a a ) } e ^ { i ( { \bf k } - { \bf k } ^ { \prime } ) \cdot { \bf R } _ { i } - i \vec { \delta } _ { a b } \cdot { \bf k } ^ { \prime } } ( 1 4 \mathrm { ~ } } } \end{array}
$$

where δt(ab) where $\delta t _ { i j } ^ { ( a b ) } ~ ( \delta t _ { i j . } ^ { ( a a ) } ) _ { . }$ is the chhange o he ing n $\mathbf { R } _ { i }$ ${ \bf R } _ { j }$ same (different) sublattices (we have written ${ \bf R } _ { j } = { \bf R } _ { i } + { \vec { \delta } }$ where $\vec { \delta } _ { a b }$ is the nearest neighbor vector, and $\vec { \delta } _ { a a }$ is the next nearest neighbor vector). Following the procedure of Sec. II.B we project out the Fourier components of the operators close to the $\mathrm { K }$ and $\mathrm { K } '$ points of the BZ using (17). If we assume that $\delta t _ { i j }$ is smooth over the lattice spacing scale, that is, it does not have an Fourier component with momentum $\mathbf { K } - \mathbf { K } ^ { \prime }$ (so the two Dirac cones are not coupled by disorder), we can rewrite (145) in real space as:

$$
\begin{array} { l } { { \displaystyle H _ { \mathrm { o d } } ~ = ~ \int d ^ { 2 } r \left[ { \cal A } ( { \bf r } ) a _ { 1 } ^ { \dagger } ( { \bf r } ) b _ { 1 } ( { \bf r } ) + \mathrm { h . c . } \right. } } \\ { { \displaystyle ~ + ~ \left. \phi ( { \bf r } ) \left( a _ { 1 } ^ { \dagger } ( { \bf r } ) a _ { 1 } ( { \bf r } ) + b _ { 1 } ^ { \dagger } ( { \bf r } ) b _ { 1 } ( { \bf r } ) \right) \right] , } } \end{array}
$$

with a similar expression for the cone 2 but with $\mathcal { A }$ replaced by $A ^ { * }$ , where,

$$
\begin{array} { l } { { \displaystyle { \mathcal A } ( { \bf r } ) ~ = ~ \sum _ { \vec { \delta } _ { a b } } \delta t ^ { ( a b ) } ( { \bf r } ) e ^ { - i \vec { \delta } _ { a b } \cdot { \bf K } } ~ } , } \\ { { \displaystyle ~ \phi ( { \bf r } ) ~ = ~ \sum _ { \vec { \delta } _ { a a } } \delta t ^ { ( a a ) } ( { \bf r } ) e ^ { - i \vec { \delta } _ { a a } \cdot { \bf K } } ~ . } } \end{array}
$$

Notice that whereas $\phi ( \mathbf { r } ) \mathbf { \Psi } = \mathbf { \Psi } \phi ^ { * } ( \mathbf { r } )$ , because of the inversion symmetry of the two triangular sublattices that

make up the honeycomb lattice, $\mathcal { A }$ is complex because of lack of inversion symmetry for nearest neighbor hopping. Hence,

$$
\begin{array} { r } { \mathcal { A } ( \mathbf { r } ) = \mathcal { A } _ { x } ( \mathbf { r } ) + i \mathcal { A } _ { y } ( \mathbf { r } ) . } \end{array}
$$

In terms of the Dirac Hamiltonian (18) we can rewrite (146) as:

$$
\begin{array} { l } { { \displaystyle H _ { \mathrm { o d } } ~ = ~ \int d ^ { 2 } r \left[ \hat { \Psi } _ { 1 } ^ { \dagger } ( { \bf r } ) \pmb { \sigma } \cdot \vec { \mathcal { A } } ( { \bf r } ) \hat { \Psi } _ { 1 } ( { \bf r } ) \right. } } \\ { { \displaystyle ~ + ~ \left. \phi ( { \bf r } ) \hat { \Psi } _ { 1 } ^ { \dagger } ( { \bf r } ) \hat { \Psi } _ { 1 } ( { \bf r } ) \right] ~ , } } \end{array}
$$

where $\vec { \mathcal { A } } = ( \mathcal { A } _ { x } , \mathcal { A } _ { y } )$ . This result shows that changes in the hopping amplitude lead to the appearance of vector, $\vec { A }$ , and scalar, $\Phi$ , potentials in the Dirac Hamiltonian. The presence of a vector potential in the problem indicates that an effective magnetic field $\vec { B } = c / ( e v _ { F } ) \nabla \times \vec { A }$ should also be present, naively implying a broken time reversal symmetry, although the original problem was time reversal invariant. This broken time reversal symmetry is not real since (150) is the Hamiltonian around only one of the Dirac cones. The second Dirac cone is related to the first by time reversal symmetry indicating that the effective magnetic field is reversed in the second cone. Therefore, there is no global broken symmetry but a compensation between the two cones.

# A. Ripples

Graphene is a one atom thick system, the extreme case of a soft membrane. Hence, just like soft membranes, it is subject to distortions of its structure either due to thermal fluctuations (as we discussed in Sec. III) or interaction with a substrate, scaffold, and absorbands (Swain and Andelman, 1999). In the first case the fluctuations are time dependent (although with time scales much longer than the electronic ones), while in the second case the distortions act as quenched disorder. In both cases, the disorder comes about because of the modification of the distance and relative angle between the carbon atoms due to the bending of the graphene sheet. This type of off-diagonal disorder does not exist in ordinary 3D solids, or even in quasi-1D or quasi-2D systems, where atomic chains and atomic planes, respectively, are embedded in a 3D crystalline structure. In fact, graphene is also very different from other soft membranes because it is (semi) metallic, while previously studied membranes were insulators.

The problem of the bending of graphitic systems and its effect on the hybridization of the $\pi$ orbitals has been studied a great deal in the context of classical minimal surfaces (Lenosky et al., 1992) and applied to fullerenes and carbon nanotubes (Kane and Mele, 1997; Tersoff, 1992; Tu and Ou-Yang, 2002; Xin et al., 2000; Zhong-can et al., 1997w). In order to understand the effect of bending on graphene, consider the situation shown in Fig.25. The bending of the graphene sheet has three main effects: the decrease of the distance between carbon atoms, a rotation of the $p _ { Z }$ orbitals (compression or dilation of the lattice are energetically costly due to the large spring constant of graphene $\approx 5 7 \ : \mathrm { \ e V / \AA ^ { 2 } }$ (Xin et al., 2000)), and a re-hybridization between $\pi$ and $\sigma$ orbitals (Eun-Ah Kim and Castro Neto, 2007). Bending by a radius $R$ decreases the distance between the orbitals from $\ell$ to $d = 2 R \sin [ \ell / ( 2 R ) ] \approx \ell - \ell ^ { 3 } / ( 2 4 R ^ { 2 } )$ for $R \ \gg \ \ell$ . The decrease in the distance between the orbitals increases the overlap between the two lobes o f the $p _ { Z }$ orbital Harrison, 1980): $V _ { p p a } ~ \approx ~ V _ { p p a } ^ { 0 } [ 1 ~ +$ $\ell ^ { 2 } / ( 1 2 R ^ { 2 } ) ]$ , where $a = \pi , \sigma$ , and $V _ { p p a } ^ { 0 }$ is the overlap for a flat graphene sheet. The rotation of the $p _ { Z }$ orbitals can be understood within the Slater-Koster formalism, namely, the rotation can be decomposed into a $p _ { z } - p _ { z }$ ( $\pi$ bond) plus a $p _ { x } \mathrm { ~ - ~ } p _ { x }$ $\sigma$ bond) hybridization with energies $V _ { p p \pi }$ and $V _ { p p \sigma }$ , respectively (Harrison, 1980): $V ( \theta ) \ = \ V _ { p p \pi } \cos ^ { 2 } ( \theta ) - V _ { p p \sigma } \sin ^ { 2 } ( \theta ) \ \approx \ V _ { p p \pi }$ (Vppπ + $V _ { p p \sigma } ) ( \ell / ( 2 R ) ) ^ { 2 }$ , leading to a decrease in the overlap. Furthermore, the rotation leads to re-hybridization between $\pi$ and $\sigma$ orbitals leading to a further shift in energy of the order of (Eun-Ah Kim and Castro Neto, 2007): $\delta \epsilon _ { \pi } \approx ( V _ { s p \sigma } ^ { 2 } + V _ { p p \sigma } ^ { 2 } ) / ( \epsilon _ { \pi } - \epsilon _ { a } )$ b

![](images/16475eb72639412d9bf349360350786e9541cf580273d8581196aff2aea748a2.jpg)  
Figure 25 Bending of the surface of graphene by a radius $R$ and its effect on the $p _ { z }$ orbitals.

In the presence of a substrate, as we discussed in Sec.III, elasticity theory predicts that graphene can be expected to follow the substrate in a smooth way. Indeed, by minimizing the elastic energy (126), (127), and (128) with respect to the height $h$ we get (Swain and Andelman, 1999):

$$
\kappa \nabla ^ { 4 } h ( { \bf r } ) - \gamma \nabla ^ { 2 } h ( { \bf r } ) + v h ( { \bf r } ) = v s ( { \bf r } ) ,
$$

that can be solved by Fourier transform:

$$
h ( { \bf k } ) = \frac { s ( { \bf k } ) } { 1 + ( \ell _ { t } k ) ^ { 2 } + ( \ell _ { c } k ) ^ { 4 } } ,
$$

where

$$
\begin{array} { l } { { \ell _ { t } ~ = ~ \left( \frac { \gamma } { v } \right) ^ { 1 / 2 } , } } \\ { { \ell _ { c } ~ = ~ \left( \frac { \kappa } { v } \right) ^ { 1 / 4 } . } } \end{array}
$$

Eq. (152) gives the height configuration in terms of the substrate profile, and $\ell _ { t }$ and $\ell _ { c }$ provide the length scales for elastic distortion of graphene on a substrate. Hence, disorder in the substrate translates into disorder in the graphene sheet (albeit restricted by elastic constraints). This picture has been confirmed by STM measurements on graphene (Ishigami et al., 2007; Stolyarova et al., 2007) in which strong correlations were found between the roughness of the substrate and the graphene topography. Ab initio band structure calculations also give support to this scenario (Dharma-Wardana, 2007).

The connection between the ripples and the electronic problem comes from the relation between the height field $h ( \mathbf { r } )$ and the local curvature of the graphene sheet, $R$ :

$$
\frac { 2 } { R ( { \bf r } ) } \approx \nabla ^ { 2 } h ( { \bf r } ) ,
$$

and, hence we see that due to bending the electrons are subject to a potential which depends on the structure of a graphene sheet (Eun-Ah Kim and Castro Neto, 2007):

$$
V ( { \bf r } ) \approx V ^ { 0 } - \alpha \left( \nabla ^ { 2 } h ( { \bf r } ) \right) ^ { 2 } ,
$$

where $\alpha$ $\alpha \approx 1 0 \ \mathrm { e V } \ \mathrm { \AA } ^ { 2 }$ ) is the constant that depends on microscopic details. The conclusion from (155) is that the Dirac fermions are scattered by ripples of the graphene sheet through a potential which is proportional to the square of the local curvature. The coupling between geometry and electron propagation is unique to graphene, and results in additional scattering and resistivity due to ripples (Katsnelson and Geim, 2008).

# B. Topological lattice defects

Structural defects of the honeycomb lattice like pentagons, heptagons and their combinations such as StoneWales defect (a combination of two pentagon-heptagon pairs) are also possible in graphene and can lead to scattering (Cortijo and Vozmediano, 2007a,b). These defects induce long range deformations, which modify the electron trajectories.

Let us consider first a disclination. This defect is equivalent to the deletion or inclusion of a wedge in the lattice. The simplest one in the honeycomb lattice is the absence of a $6 0 ^ { \circ }$ wedge. The resulting edges can be glued in such a way that all sites remain three-fold coordinated. The honeycomb lattice is recovered everywhere, except at the apex of the wedge, where a fivefold ring, a pentagon, is formed. One can imagine a situation where the nearest neighbor hoppings are unchanged. Nevertheless, the existence of a pentagon implies that the two sublattices of the honeycomb structure can no longer be defined. A trajectory around the pentagon after a closed circuit has to change the sublattice index.

The boundary conditions imposed at the edges of a disclination are sketched in Fig. 26, showing the identification of sites from different sublattices. In addition, the wavefunctions at the $K$ and $K ^ { \prime }$ points are exchanged when moving around the pentagon.

![](images/10e18c57e2281f677f7759432e3cd4a9200c394dd93ec61230e36cab21d36193.jpg)  
Figure 26 (Color online) Sketch of the boundary conditions associated to a disclination (pentagon) in the honeycomb lattice.

Far away from the defect, a slow rotation of the components of the spinorial wavefunction can be described by a gauge field which acts on the valley and sublattice indices (González et al., 1992, 1993b). This gauge field is technically non-abelian, although a transformation can be defined which makes the resulting Dirac equation equivalent to one with an effective abelian gauge field (González et al., 1993b). The final continuum equation gives a reasonable description of the electronic spectrum of fullerenes of different sizes (González et al., 1992, 1993b), and other structures which contain pentagons (Kolesnikov and Osipov, 2004, 2006; Lammert and Crespi, 2004; LeClair, 2000; Osipov et al., 2003). It is easy to see that an heptagon leads to the opposite effective field.

An in-plane dislocation, that is, the inclusion of a semiinfinite row of sites, can be considered as induced by a pentagon and a heptagon together. The non-abelian field described above is canceled away from the core. A closed path implies a shift by one (or more) lattice spacings. The wavefunctions at the $K$ and $K ^ { \prime }$ points acquire pphases, $e ^ { \pm 2 \pi i / 3 }$ , under a translation by one lattice unit. Hence, the description of a dislocation in the continuum limit requires an (abelian) vortex of charge $\pm 2 \pi / 3$ at its core. Dislocations separated over distances of the order of $d$ lead to an effective flux through an area of perimeter $\it l$ of the order of (Morpurgo and Guinea, 2006):

$$
\Phi \sim \frac { d } { k _ { \mathrm { { F } } } l ^ { 2 } }
$$

where $k _ { \mathrm { { F } } }$ is the Fermi vector of the electrons.

In general, a local rotation of the axes of the honeycomb lattice induces changes in the hopping which lead to mixing of the $K$ and $K ^ { \prime }$ wavefunctions, leading to a gauge field like the one induced by a pentagon (González et al., 2001). Graphene samples with disclinations and dislocations are feasible in particular substrates (Couraux et al., 2008), and gauge fields related to the local curvature are then expected to play a crucial role in such structures. The resulting electronic structure can be analyzed using the theory of quantum mechanics in curved space (Birrell and Davies, 1982; Cortijo and Vozmediano, 2007a,b; de Juan et al., 2007).

# C. Impurity states

Point defects, similar to impurities and vacancies, can nucleate electronic states in their vicinity. Hence, a concentration of $n _ { i }$ impurities per carbon atom leads to a change in the electronic density of the order of $n _ { i }$ . The corresponding shift in the Fermi energy is $\epsilon _ { \mathrm { F } } \simeq v _ { \mathrm { F } } \sqrt { n _ { i } }$ . In addition, impurities lead to a finite elastic mean free path, $l _ { \mathrm { e l a s } } \simeq a n _ { i } ^ { - 1 / 2 }$ , and to an elastic scattering time $\tau _ { \mathrm { e l a s } } \simeq ( v _ { \mathrm { F } } n _ { i } ) ^ { - 1 }$ . Hence, the regions with few impurities can be considered low-density metals in the dirty limit, as Telas $\tau _ { \mathrm { e l a s } } ^ { - 1 } \simeq \epsilon _ { \mathrm { F } }$ .

The Dirac equation allows for localized solutions that satisfy many possible boundary conditions. It is known that small circular defects result in localized and semilocalized states (Dong et al., 1998), that is, states whose wavefunction decays as $1 / r$ as a function of the distance from the center of the defect. A discrete version of these states can be realized in a nearest neighbor tightbinding model with unitary scatterers such as vacancies (Pereira et al., 2006). In the continuum, the Dirac equation (19) for the wavefunction, $\psi ( \mathbf { r } ) = ( \phi _ { 1 } ( \mathbf { r } ) , \phi _ { 2 } ( \mathbf { r } ) )$ , can be written as:

$$
\begin{array} { r l r } { \partial _ { w } \phi _ { 1 } ( w , w ^ { * } ) } & { { } = } & { 0 , } \\ { \partial _ { w ^ { * } } \phi _ { 2 } ( w , w ^ { * } ) } & { { } = } & { 0 , } \end{array}
$$

where $w = x + i y$ is a complex number. These equations indicate that at zero energy the components of the wavefunction can only be either holomorphic or antiholomorphic with respect to the variable $w$ (that is, $\phi _ { 1 } ( w , w ^ { * } ) = \phi _ { 1 } ( w ^ { * } )$ and $\phi _ { 2 } ( w , w ^ { * } ) = \phi _ { 2 } ( w ) )$ . Since the boundary conditions require that the wavefunction vanishes at infinity the only possible solutions have the form: $\Psi _ { K } ( \tilde { \mathbf { r } } ) \propto ( 1 / ( x + i y ) ^ { n } , 0 )$ or $\Psi _ { K ^ { \prime } } ( \tilde { \mathbf { r } } ) \propto ( 0 , 1 / ( x - i y ) ^ { n } )$ . The wavefunctions in the discrete lattice must be real, and at large distances the actual solution found near a vacancy tends to a superposition of two solutions formed from wavefunctions from the two valleys with equal weight, in a way similar to the mixing at armchair edges (Brey and Fertig, 2006b).

The construction of the semi-localized state around a vacancy in the honeycomb lattice can be extended to other discrete models which leads to the Dirac equation in the continuum limit. A particular case is the nearest neighbor square lattice with half flux per plaquette, or the nearest neighbor square lattice with two flavors per site. The latter has been extensively studied in relation to the effects of impurities on the electronic structure of d-wave superconductors (Balatsky et al., 2006), and numerical results are in agreement with the existence of this solution. As the state is localized on one sublattice

$$
3 0 0 - 3 0 0 = 3 0 0
$$

Figure 27 (Color online) Sketch of a rough graphene surface. The full line gives the boundary beyond which the lattice can be considered undistorted. The number of mid-gap states changes depending on a difference in the number of removed sites for two sublattices.

only, the solution can be generalized for the case of two vacancies.

# D. Localized states near edges, cracks, and voids

Localized states can be defined at edges where the number of atoms in the two sublattices is not compensated. The number of them depend on details of the edge. The graphene edges can be strongly deformed, due to the bonding of other atoms to carbon atoms at the edges. These atoms should not induce states in the graphene $\pi$ band. In general, a boundary inside the graphene material will exist, as sketched in Fig. 27, beyond which the $s p ^ { 2 }$ hybridization is well defined. If this is the case, the number of mid-gap states near the edge is roughly proportional to the difference in sites between the two sublattices near this boundary.

Along a zigzag edge there is one localized state per three lattice units. This implies that a precursor structure for localized states at the Dirac energy can be found in ribbons or constrictions of small lengths (Muñoz-Rojas et al., 2006), which modifies the electronic structure and transport properties.

Localized solutions can also be found near other defects that contain broken bonds or vacancies. These states do not allow an analytical solution, although, as discussed above, the continuum Dirac equation is compatible with many boundary conditions, and it should describe well localized states that vary slowly over distances comparable to the lattice spacing. The existence of these states can be investigated by analyzing the scaling of the spectrum near a defect as a function of the size of the system, $L$ (Vozmediano et al., 2005). A number of small voids and elongated cracks show states whose energy scales as $L ^ { - 2 }$ , while the energy of extended states scales as $L ^ { - 1 }$ . A state with energy scaling $L ^ { - 2 }$ is compatible with continuum states for which the modulus of the wavefunction decays as $r ^ { - 2 }$ as a function of the distance from the defect.

# E. Self-doping

The band structure calculations discussed in previous sections show that the electronic structure of a single graphene plane is not strictly symmetrical in energy (Reich et al., 2002). The absence of electron-hole symmetry shifts the energy of the states localized near impurities above or below the Fermi level, leading to a transfer of charge from/to the clean regions. Hence, the combination of localized defects and the lack of perfect electronhole symmetry around the Dirac points leads to the possibility of self-doping, in addition to the usual scattering processes.

Extended lattice defects, like edges, grain boundaries, or micro-cracks, are likely to induce a number of electronic states proportional to their length, $L / a$ , where $a$ is of the order of the lattice constant. Hence, a distribution of extended defects of length $L$ at a distance equal to $L$ itself gives rise to a concentration of $L / a$ carriers per carbon in regions of size of the order of $( L / a ) ^ { 2 }$ . The resulting system can be considered a metal with a low density of carriers, $n _ { \mathrm { c a r r i e r } } \propto a / L$ per unit cell, and an elastic mean free path $l _ { \mathrm { e l a s } } \simeq L$ . Then, we obtain:

$$
\begin{array} { r c l } { { \epsilon _ { \mathrm { { F } } } } } & { { \simeq } } & { { \frac { v _ { \mathrm { { F } } } } { \sqrt { a L } } } } \\ { { } } & { { } } & { { { } } } \\ { { \frac { 1 } { \tau _ { \mathrm { { e l a s } } } } } } & { { \simeq } } & { { \frac { v _ { \mathrm { { F } } } } { L } } } \end{array}
$$

and, therefore, $( \tau _ { \mathrm { e l a s } } ) ^ { - 1 } \ll \epsilon _ { \mathrm { F } }$ when $a / L \ll 1$ . Hence, the existence of extended defects leads to the possibility of self-doping but maintaining most of the sample in the clean limit. In this regime, coherent oscillations of transport properties are expected, although the observed electronic properties may correspond to a shifted Fermi energy with respect to the nominally neutral defect-free system.

One can describe the effects that break electron-hole symmetry near the Dirac points in terms of a finite nextnearest neighbor hopping between $\pi$ orbitals, $t ^ { \prime }$ , in (148). Consider, for instance, electronic structure of a ribbon of width $L$ terminated by zigzag edges, which, as discussed, lead to surface states for $t ^ { \prime } = 0$ . The translational symmetry along the axis of the ribbon allows us to define bands in terms of a wavevector parallel to this axis. On the other hand, the localized surface bands, extending from $k _ { \parallel } = ( 2 \pi ) / 3$ to $k _ { \parallel } = - ( 2 \pi ) / 3$ acquire a dispersion of order $t ^ { \prime }$ . Hence, if the Fermi energy remains unchanged at the position of the Dirac points ( $\epsilon _ { \mathrm { { D i r a c } } } = - 3 t ^ { \prime }$ ), this band will be filled, and the ribbon will no longer be charge neutral. In order to restore charge neutrality, the Fermi level needs to be shifted by an amount of the order of $t ^ { \prime }$ . As a consequence, some of the extended states near the Dirac points are filled, leading to the phenomenon of self-doping. The local charge is a function of distance to the edges, setting the Fermi energy so that the ribbon is globally neutral. Note that the charge transferred to the surface states is mostly localized near the edges of the system.

![](images/959e1322492f1abcdb190d3a0aa22650ce50d8247ae90d9fae14d10e359ca374.jpg)  
Figure 28 (Color online) Top: Self-consistent analysis of the displaced charge density (in units of number of electrons per carbon) is shown as a continuous line, and the corresponding electrostatic potential (in units of $t$ ) is shown as a dashed line, for a graphene ribbon with periodic boundary conditions along the zig-zag edge (with a length of $L = 9 6 0 a$ ) and with a circumference of size $W = 8 0 \sqrt { 3 } a$ . The inset shows the charge density near the edge. Due to the presence of the edge, there is a displaced charge in the bulk (bottom panel) that is shown as a function of width $W$ . Notice that the displaced charge vanishes in the bulk limit ( $W \to \infty$ ), in agreement with (161). Adapted from Peres et al., 2006c.

The charge transfer is suppressed by electrostatic effects, as large deviations from charge neutrality have an associated energy cost (Peres et al., 2006c). In order to study these charging effects we add to the free-electron Hamiltonian (5) the Coulomb energy of interaction between electrons:

$$
H _ { I } = \sum _ { i , j } U _ { i , j } n _ { i } n _ { j } ,
$$

where $\begin{array} { r } { n _ { i } = \sum _ { \sigma } ( a _ { i , \sigma } ^ { \dagger } a _ { i , \sigma } + b _ { i , \sigma } ^ { \dagger } b _ { i , \sigma } ) } \end{array}$ is the number operator at site $\mathbf { R } _ { i }$ , and

$$
U _ { i , j } = \frac { e ^ { 2 } } { \epsilon _ { 0 } | { \bf R } _ { i } - { \bf R } _ { j } | } ,
$$

is the Coulomb interaction between electrons. We expect, on physics grounds, that an electrostatic potential builds up at the edges, shifting the position of the surface states, and reducing the charge transferred to/from them. The potential at the edge induced by a constant doping $\delta$ per carbon atom is roughly, $\sim ( \delta e ^ { 2 } / a ) ( W / a )$ $\delta e ^ { 2 } / a$ is the Coulomb energy per carbon), and $W$ the width of the ribbon ( $W / a$ is the number of atoms involved). The charge transfer is stopped when the potential shifts the localized states to the Fermi energy, that is, when $t ^ { \prime } \approx$ $( e ^ { 2 } / a ) ( W / a ) \delta$ . The resulting self-doping is therefore

$$
\delta \sim \frac { t ^ { \prime } a ^ { 2 } } { e ^ { 2 } W } ,
$$

![](images/783e3ad49419d9a76068eb0b07883c5c9fe7df912db845d508622257fb0fac86.jpg)  
Figure 29 (Color online) Gauge field induced by a simple elastic strain. Top: The hopping along the horizontal bonds is assumed to be changed on the right hand side of the graphene lattice, defining a straight boundary between the unperturbed and perturbed regions (green dashed line). Bottom: The modified hopping acts like a constant gauge field, which displaces the Dirac cones in opposite directions at the $K$ and $K ^ { \prime }$ points of the Brillouin zone. The conservation of energy and momentum parallel to the boundary leads to a deflection of electrons by the boundary.

that vanishes when $W \to \infty$

We treat Hamiltonian (159) within the Hartree approximation (that is, we replace $H _ { I }$ by $\begin{array} { r } { H _ { \mathrm { M . F . } } \ = \ \sum _ { i } V _ { i } n _ { i } } \end{array}$ where $\begin{array} { r } { V _ { i } ~ = ~ \sum _ { j } U _ { i , j } \langle n _ { j } \rangle } \end{array}$ , and solve the problem selfconsistently for $\left. { n _ { i } } \right.$ ). Numerical results for graphene ribbons of length $L = 8 0 \sqrt { 3 } a$ and different widths are shown in Fig. 28 ( $t ^ { \prime } / t \ = \ 0 . 2$ and $e ^ { 2 } / a \ = \ 0 . 5 t$ . The largest width studied is $\sim 0 . 1 \mu \mathrm { m }$ , and the total number of carbon atoms in the ribbon is $\approx 1 0 ^ { 5 }$ . Notice that as $W$ increases, the self-doping decreases indicating that, for a perfect graphene plane ( $W \to \infty$ ), the self-doping effect disappears. For realistic parameters, we find that the amount of self-doping is $1 0 ^ { - 4 } - 1 0 ^ { - 5 }$ electrons per unit cell for sizes $0 . 1 - 1 \mu \mathrm { m }$ .

# F. Vector potential and gauge field disorder

As discussed in Sec. IV, lattice distortions modify the Dirac equation that describes the low energy band structure of graphene. We consider here deformations that change slowly on the lattice scale, so that they do not mix the two inequivalent valleys. As shown earlier, perturbations that hybridize the two sublattices lead to terms that change the Dirac Hamiltonian from $v _ { F } \sigma \cdot { } \mathbf { k }$ into $v _ { \mathrm { F } } \pmb { \sigma } \cdot \mathbf { k } + \pmb { \sigma } \cdot \mathbf { A }$ . Hence, the vector A can be thought of as if induced by an effective gauge field, A. In order to preserve time reversal symmetry, this gauge field must have opposite signs at the two Dirac cones, $\mathbf { A } _ { K } = - \mathbf { A } _ { K ^ { \prime } }$ .

A simple example is a distortion that changes the hopping between all bonds along a given axis of the lattice. Let us assume that the sites at the ends of those bonds define the unit cell, as sketched in Fig. 29. If the distortion is constant, its only effect is to displace the Dirac points away from the BZ corners. The two inequivalent points are displaced in opposite directions. This uniform distortion is the equivalent of a constant gauge field, which does not change the electronic spectrum. The situation changes if one considers a boundary that separates two domains where the magnitude of the distortion is different. The shift of the Dirac points leads to a deflection of the electronic trajectories that cross the boundary, as also sketched in Fig. 29. The modulation of the gauge field leads to an effective magnetic field, which is of opposite sign for the two valleys.

Table II Estimates of the effective magnetic length, and effective magnetic fields generated by the deformations considered in this section. The intrinsic curvature entry also refers to the contribution from topological defects.   

<table><tr><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>lB</td><td rowspan=1 colspan=1>Bh=1 nm, l=10nm, a=0.1nm</td></tr><tr><td rowspan=1 colspan=1>Intrinsic curvature</td><td rowspan=1 colspan=1>l(</td><td rowspan=1 colspan=1>0.06T</td></tr><tr><td rowspan=1 colspan=1>Extrinsic curvature</td><td rowspan=1 colspan=1>\5 2</td><td rowspan=1 colspan=1>0.006T</td></tr><tr><td rowspan=1 colspan=1>Elastic strains</td><td rowspan=1 colspan=1>1 al</td><td rowspan=1 colspan=1>6T</td></tr></table>

We have shown in Section IV.B how topological lattice defects, such as disclinations and dislocations, can be described by an effective gauge field. Those defects can only exist in graphene sheets that are intrinsically curved, and the gauge field only depends on topology of the lattice. Changes in the nearest neighbor hopping also lead to effective gauge fields. We consider next two physical processes that induce effective gauge fields: (i) changes in the hopping induced by hybridization between $\pi$ and $\sigma$ bands, which arise in curved sheets, and $( i i )$ changes in the hopping due to modulation in the bond length, which is associated with elastic strain. The strength of these fields depends on parameters that describe the value of the $\boldsymbol { \mathscr { n } }$ $\sigma$ hybridization, and the dependence of hopping on the bond length.

A comparison of the relative strengths of the gauge fields induced by intrinsic curvature, $\pi - \sigma$ hybridization (extrinsic curvature), and elastic strains, arising from a ripple of typical height and size is given in Table II.

# 1Gauge field induced by curvature

As we discussed in Sec. IV.A, when the $\pi$ orbitals are not parallel, the hybridization between them depends on their relative orientation. The angle $\theta _ { i }$ determines the relative orientation of neighboring orbitals at some position $\mathbf { r } _ { i }$ in the graphene sheet. The value of $\theta _ { i }$ depends on the local curvature of the layer. The relative angle of rotation of two $p _ { z }$ orbitals at position $\mathbf { r } _ { i }$ and $\mathbf { r } _ { j }$ can be written as: $\cos ( \theta _ { i } - \theta _ { j } ) = { \bf N } _ { i } \cdot { \bf N } _ { j }$ , where $\mathbf { N } _ { i }$ is the unit vector perpendicular to the surface, defined in (125). If $\mathbf { r } _ { j } = \mathbf { r } _ { i } + \mathbf { u } _ { i j }$ we can write:

$$
\mathbf { N } _ { i } \cdot \mathbf { N } _ { j } \approx 1 + \mathbf { N } _ { i } \cdot \left[ ( \mathbf { u } _ { i j } \cdot \nabla ) \mathbf { N } _ { i } \right] + \frac { 1 } { 2 } \mathbf { N } _ { i } \cdot \left[ ( \mathbf { u } _ { i j } \cdot \nabla ) ^ { 2 } \mathbf { N } _ { i } \right] ,
$$

where we assume smoothly varying $\mathbf { N } ( \mathbf { r } )$ . We use (125) in terms of the height field $h ( \mathbf { r } )$ ${ \bf N } ( { \bf r } ) \approx { \bf z } - \nabla h ( { \bf r } ) -$ $( \nabla h ) ^ { 2 } \mathbf { z } / 2 )$ to rewrite (162) as:

$$
{ \bf N } _ { i } \cdot { \bf N } _ { j } \approx 1 - \frac { 1 } { 2 } [ ( { \bf u } _ { i j } \cdot \nabla ) \nabla h ( { \bf r } _ { i } ) ] ^ { 2 } .
$$

Hence, bending of the graphene sheet leads to a modification of the hopping amplitude between different sites in the form:

$$
\delta t _ { i j } \approx - \frac { t _ { i j } ^ { 0 } } { 2 } [ ( \mathbf { u } _ { i j } \cdot \nabla ) \nabla h ( \mathbf { r } _ { i } ) ] ^ { 2 } ,
$$

where $t _ { i j } ^ { 0 }$ is the bare hopping energy. A similar effect leads to changes the electronic states in a carbon nanotubes (Kane and Mele, 1997). Using the results of Sec. IV, namely (147), we can now see that a vector potential is generated for nearest neighbor hopping ( $\mathbf { u } = \vec { \delta } _ { a b }$ ) (Eun-Ah Kim and Castro Neto, 2007):

$$
\begin{array} { l } { { A _ { x } ^ { ( h ) } ~ = ~ - \frac { 3 E _ { a b } a ^ { 2 } } { 8 } \left[ ( \partial _ { x } ^ { 2 } h ) ^ { 2 } - ( \partial _ { y } ^ { 2 } h ) ^ { 2 } \right] } } \\ { { A _ { y } ^ { ( h ) } ~ = ~ \frac { 3 E _ { a b } a ^ { 2 } } { 4 } \left( \partial _ { x } ^ { 2 } h + \partial _ { y } ^ { 2 } h \right) \partial _ { x } h \partial _ { y } h } } \end{array}
$$

where the coupling constant $E _ { a b }$ depends on microscopic details (Eun-Ah Kim and Castro Neto, 2007). The flux of effective magnetic field through a ripple of lateral dimension $\it l$ and height $h$ is given approximately by:

$$
\Phi \approx { \frac { E _ { a b } a ^ { 2 } h ^ { 2 } } { v _ { \mathrm { { F } } } l ^ { 3 } } }
$$

where the radius of curvature is $R ^ { - 1 } \approx h l ^ { - 2 }$ . For a ripple with $l \approx 2 0 \mathrm { n m } , h \approx 1 \mathrm { n m }$ , taking $E _ { a b } / v _ { \mathrm { F } } \approx 1 0 ~ \mathrm { \AA } ^ { - 1 }$ , we find $\Phi \approx 1 0 ^ { - 3 } \Phi _ { 0 }$ .

# 2. Elastic strain

The elastic free energy for graphene can be written in terms of the in-plane displacement $\mathbf { u } ( \mathbf { r } ) = ( u _ { x } , u _ { y } )$ as:

$$
F [ { \bf u } ] = \frac { 1 } { 2 } \int d ^ { 2 } r \left[ ( B - \mathcal { G } ) ( \sum _ { i = 1 , 2 } u _ { i i } ) ^ { 2 } + 2 \mathcal { G } \sum _ { i , j = 1 , 2 } u _ { i j } ^ { 2 } \right] (
$$

where $\boldsymbol { \beta }$ is the bulk modulus, $\mathcal { G }$ is the shear modulus, and

$$
u _ { i j } = \frac { 1 } { 2 } \left( \frac { \partial u _ { i } } { \partial x _ { j } } + \frac { \partial u _ { j } } { \partial x _ { j } } \right) ,
$$

is the strain tensor ( $x _ { 1 } = x$ and $x _ { 2 } = y$ ).

There are many types of static deformation of the honeycomb lattice which can affect the propagation of Dirac fermions. The simplest one is due to changes in the area of the unit cell either due to dilation or contraction. Changes in the unit cell area lead to local changes in the density of electrons and, therefore, local changes in the chemical potential in the system. In this case, their effect is similar to the one found in (148), and we must have:

$$
\phi _ { \mathrm { d p } } ( \mathbf { r } ) = g ( u _ { x x } + u _ { y y } ) ,
$$

and their effect is diagonal in the sublattice index.

The nearest neighbor hopping depends on the length of the carbon bond. Hence, elastic strains that modify the relative orientation of the atoms also lead to an effective gauge field, which acts on each $K$ point separately, as first discussed in relation to carbon nanotubes (Mañes, 2007; Suzuura and Ando, 2002b). Consider two carbon atoms located in two different sublattices in the same unit cell at $\mathbf { R } _ { i }$ . The change in the local bond length can be written as:

$$
\delta u _ { i } = \frac { \vec { \delta } _ { a b } } { a } \cdot \left[ \mathbf { u } _ { \mathrm { A } } ( \mathbf { R } _ { i } ) - \mathbf { u } _ { \mathrm { B } } ( \mathbf { R } _ { i } + \vec { \delta } _ { a b } ) \right] .
$$

The local displacements of the atoms in the unit cell can be related to $\mathbf { u } ( \mathbf { r } )$ by (Ando, 2006a):

$$
( \vec { \delta } _ { a b } \cdot \nabla ) { \bf u } = \kappa ^ { - 1 } ( { \bf u } _ { \mathrm { A } } - { \bf u } _ { \mathrm { B } } ) ,
$$

where $\kappa$ is a dimensionless quantity that depends on microscopic details. Changes in the bond length lead to changes in the hopping amplitude:

$$
t _ { i j } \approx t _ { i j } ^ { 0 } + \frac { \partial t _ { i j } } { \partial a } \delta u _ { i } ,
$$

and we can write:

$$
\delta t ^ { ( a b ) } ( \mathbf { r } ) \approx \beta \frac { \delta u ( \mathbf { r } ) } { a } ,
$$

where

$$
\beta = \frac { \partial t ^ { ( a b ) } } { \partial \ln ( a ) } .
$$

Substituting (170) into (173) and the final result into (147), one finds (Ando, 2006a):

$$
\begin{array} { r c l } { { } } & { { } } & { { \displaystyle { \mathcal A } _ { x } ^ { ( s ) } = \frac { 3 } { 4 } \beta \kappa \left( u _ { x x } - u _ { y y } \right) , } } \\ { { } } & { { } } & { { \displaystyle { \mathcal A } _ { y } ^ { ( s ) } = \frac { 3 } { 2 } \beta \kappa u _ { x y } . } } \end{array}
$$

We assume that the strains induced by a ripple of dimension $\it { l }$ and height $h$ scale as $u _ { i j } \sim ( h / l ) ^ { 2 }$ . Then, using $\beta / v _ { \mathrm { F } } \approx a ^ { - 1 } \sim 1 \mathrm { \AA } ^ { - 1 }$ , we find that the total flux through a ripple is:

$$
\Phi \approx { \frac { h ^ { 2 } } { a l } } .
$$

For ripples such that $h \sim 1 \mathrm { n m }$ and $l \sim 2 0 \mathrm { n m }$ , this estimate gives $\Phi \sim 1 0 ^ { - 1 } \Phi _ { 0 }$ in reasonable agreement with the estimates in ref. (Morozov et al., 2006).

The strain tensor must satisfy some additional constraints, as it is derived from a displacement vector field. These constraints are called Saint Venant compatibility conditions (Landau and Lifshitz, 1959):

$$
W _ { i j k l } = \frac { \partial u _ { i j } } { \partial x _ { k } \partial x _ { l } } + \frac { \partial u _ { k l } } { \partial x _ { i } \partial x _ { j } } - \frac { \partial u _ { i l } } { \partial x _ { j } \partial x _ { k } } - \frac { \partial u _ { j k } } { \partial x _ { i } \partial x _ { l } } = 0 .
$$

An elastic deformation changes the distances in the crystal lattice and can be considered as a change in the metric:

$$
g _ { i j } = \delta _ { i j } + u _ { i j }
$$

The compatibility equations (177) are equivalent to the condition that the curvature tensor derived from (178) is zero. Hence, a purely elastic deformation cannot induce intrinsic curvature in the sheet, which only arises from topological defects. The effective fields associated with elastic strains can be large (Morozov et al., 2006), leading to significant changes in the electronic wavefunctions. An analysis of the resulting state, and the possible instabilities that may occur can be found in (Guinea et al., 2007).

# 3. Random gauge fields

The preceding discussion suggests that the effective fields associated with lattice defects can modify significantly the electronic properties. This is the case when the fields do not change appreciably on scales comparable to the (effective) magnetic length. The general problem of random gauge fields for Dirac fermions has been extensively analyzed before the current interest in graphene, as the topic is also relevant for the IQHE (Ludwig et al., 1994) and d-wave superconductivity (Nersesyan et al., 1994). The one electron nature of this two dimensional problem makes it possible, at the Dirac energy, to map it onto models of interacting electrons in one dimension, where many exact results can be obtained (Castillo et al., 1997). The low energy density of states, $\rho ( \omega )$ , acquires an anomalous exponent, $\rho ( \omega ) \propto | \omega | ^ { 1 - \Delta }$ , where $\Delta > 0$ . The density of states is enhanced near the Dirac energy, reflecting the tendency of disorder to close gaps. For sufficiently large values of the random gauge field, a phase transition is also possible (Chamon et al., 1996; Horovitz and Doussal, 2002).

Perturbation theory shows that random gauge fields are a marginal perturbation at the Dirac point, leading to logarithmic divergences. These divergences tend to have the opposite sign with respect to those induced by the Coulomb interaction (see Sec. V.B). As a result, a renormalization group (RG) analysis of interacting electrons in a random gauge field suggests the possibility of non-trivial phases (Aleiner and Efetov, 2006; Altland,

2006; Dell'Anna, 2006; Foster and Ludwig, 2006a,b; Khveshchenko, 2007; Nomura et al., 2007; Stauber et al., 2005), where interactions and disorder cancel each other.

# G. Coupling to magnetic impurities

Magnetic impurities in graphene can be introduced chemically by deposition and intercalation (Calandra and Mauri, 2007; Uchoa et al., 2007), or self-generated by the introduction of defects (Kumazaki and Hirashima, 2006, 2007). The energy dependence of the density of states in graphene leads to changes in the formation of a Kondo resonance between a magnetic impurity and the graphene electrons. The vanishing of the density of states at the Dirac energy implies that a Kondo singlet in the ground state is not formed unless the exchange coupling exceeds a critical value, of the order of the electron bandwidth, a problem already studied in connection with magnetic impurities in d-wave superconductors (Cassanello and Fradkin, 1996, 1997; Fritz et al., 2006; Polkovnikov, 2002; Polkovnikov et al., 2001). For weak exchange couplings, the magnetic impurity remains unscreened. An external gate changes the chemical potential, allowing for a tuning of the Kondo resonance (Sengupta and Baskaran, 2007). The situation changes significantly if the scalar potential induced by the magnetic impurity is taken into account. This potential that can be comparable to the bandwidth allows the formation of mid-gap states and changes the phase-shift associated to spin scattering (Hentschel and Guinea, 2007). These phase-shifts have a weak logarithmic dependence on the chemical potential, and a Kondo resonance can exist, even close to the Dirac energy.

The RKKY interaction between magnetic impurities is also modified in graphene. At finite fillings, the absence of intra-valley backscattering leads to a reduction of the Friedel oscillations, which decay as $\sin ( 2 k _ { \mathrm { F } } r ) / | r | ^ { 3 }$ (Ando, 2006b; Cheianov and Fal'ko, 2006; Wunsch et al., 2006). This effect leads to an RKKY interaction, at finite fillings, which oscillate and decay as $| r | ^ { - 3 }$ . When intervalley scattering is included, the interaction reverts to the usual dependence on distance in two dimensions, $| r | ^ { - 2 }$ (Cheianov and Fal'ko, 2006). At half-filling extended defects lead to an RKKY interaction with an $| r | ^ { - 3 }$ dependence (Dugaev et al., 2006; Vozmediano et al., 2005). This behavior is changed when the impurity potential is localized on atomic scales (Brey et al., 2007; Saremi, 2007), or for highly symmetrical couplings (Saremi, 2007).

# H. Weak and strong localization

In sufficiently clean systems, where the Fermi wavelength is much shorter than the mean free path, $k _ { \mathrm { F } } l \gg 1$ electronic transport can be described in classical terms, assuming that electrons follow well defined trajectories. At low temperatures, when electrons remain coherent over long distances, quantum effects lead to interference corrections to the classical expressions for the conductivity, the weak localization correction (Bergman, 1984; Chakravarty and Schmid, 1986). These corrections are usually due to the positive interference between two paths along closed loops, traversed in opposite directions. As a result, the probability that the electron goes back to the origin is enhanced, so that quantum corrections decrease the conductivity. These interferences are suppressed for paths longer than the dephasing length, $l _ { \phi }$ , determined by interactions between the electron and environment. Interference effects can also be suppressed by magnetic fields that break down time reversal symmetry and adds a random relative phase to the process discussed above. Hence, in most metals, the conductivity increases when a small magnetic field is applied (negative magnetoresistance).

Graphene is special in this respect, due to the chirality of its electrons. The motion along a closed path induces a change in the relative weight of the two components of the wavefunction, leading to a new phase, which contributes to the interference processes. If the electron traverses a path without being scattered from one valley to the other, this (Berry) phase changes the sign of the amplitude of one path with respect to the timereversed path. As a consequence, the two paths interfere destructively, leading to a suppression of backscattering (Suzuura and Ando, 2002a). Similar processes take place in materials with strong spin orbit coupling, as the spin direction changes along the path of the electron (Bergman, 1984; Chakravarty and Schmid, 1986). Hence, if scattering between valleys in graphene can be neglected, one expects a positive magnetoresistance, i. e., weak anti-localization. In general, intra- and intervalley elastic scattering can be described in terms of two different scattering times, $\tau _ { i n t r a }$ and $\tau _ { i n t e r }$ , so that if $\tau _ { i n t r a } \ll \tau _ { i n t e r }$ one expects weak anti-localization processes, while if $\tau _ { i n t e r } \ll \tau _ { i n t r a }$ ordinary weak localization will take place. Experimentally, localization effects are always strongly suppressed close to the Dirac point but can be partially or, in rare cases, completely recovered at high carrier concentrations, depending on a particular single-layer sample (Morozov et al., 2006; Tikhonenko et al., 2007). Multilayer samples exhibit an additional positive magnetoresistance in higher magnetic fields, which can be attribued to classical changes in the current distribution due to a vertical gradient of concentration (Morozov et al., 2006) and anti-localization effects (Wu et al., 2007).

The propagation of an electron in the absence of intervalley scattering can be affected by the effective gauge fields induced by lattice defects and curvature. These fields can suppress the interference corrections to the conductivity (Morozov et al., 2006; Morpurgo and Guinea, 2006). In addition, the description in terms of free Dirac electrons is only valid near the neutrality point. The

Fermi energy acquires a trigonal distortion away from the Dirac point, and backward scattering within each valley is no longer completely suppressed (McCann et al., 2006), leading to a further suppression of anti-localization effects at high dopings. Finally, the gradient of external potentials induce a small asymmetry between the two sublattices (Morpurgo and Guinea, 2006). This effect will also contribute to reduce anti-localization, without giving rise to localization effects.

The above analysis has to be modified for a graphene bilayer. Although the description of the electronic states requires a two component spinor, the total phase around a closed loop is $2 \pi$ , and backscattering is not suppressed (Kechedzhi et al., 2007). This result is consistent with experimental observations, which show the existence of weak localization effects in a bilayer (Gorbachev et al., 2007).

When the Fermi energy is at the Dirac point, a replica analysis shows that the conductivity approaches a universal value of the order of $e ^ { 2 } / h$ (Fradkin, 1986a,b). This result is valid when intervalley scattering is neglected (Ostrovsky et al., 2006, 2007; Ryu et al., 2007). Localization is induced when these terms are included (Aleiner and Efetov, 2006; Altland, 2006), as also confirmed by numerical calculations (Louis et al., 2007). Interaction effects tend to suppress the effects of disorder. The same result, namely a conductance of the order of $e ^ { 2 } / h$ , is obtained for disordered graphene bilayers where a self-consistent calculation leads to universal conductivity at the neutrality point (Katsnelson, 2007c; Nilsson et al., 2006a, 2007a). In a biased graphene bilayer, the presence of impurities leads to the appearance of impurity tails in the density of states due to the creation of midgap states which are sensitive to the applied electric field that opens the gap between the conduction and valence bands (Nilsson and Castro Neto, 2007).

One should point out that most of the calculations of transport properties assume self-averaging, that is, that one can exchange a problem with lack of translational invariance by an effective medium system with damping. Obviously this procedure only works when the disorder is weak and the system is in the metallic phase. Close to the localized phase this procedure breaks down, the system divides itself into regions of different chemical potential and one has to think about transport in real space, usually described in terms of percolation (Cheianov et al., 2007b; Shklovskii, 2007). Single electron transistor (SET) measurements of graphene show that this seems to be the situation in graphene at halffilling (Martin et al., 2007).

Finally, we should point out that graphene stacks suffer from another source of disorder, namely, c-axis disorder that is either due to impurities between layers or rotation of graphene planes relative to each other. In either case the in-plane and out-of-plane transport is directly affected. This kind of disorder has been observed experimentally by different techniques (Bar et al., 2007; Hass et al., 2007b). In the case of the bilayer, the rotation of planes changes substantially the spectrum restoring the Dirac fermion description (Lopes dos Santos et al., 2007). The transport properties in the out of plane direction are determined by the in: op re $\begin{array} { r } { \hat { \bf j } _ { n , n + 1 } = i t \sum ( c _ { A , n , s } ^ { \dagger } c _ { A , n + 1 , s } - } \end{array}$ $c _ { A , n + 1 , s } ^ { \dagger } c _ { A , n , s } )$ $n$ $A$ terlayer hopping $t$ . If we only consider hopping between nearest neighbor sites in consecutive layers, these sites belong to one of the two sublattices in each layer.

In a multilayer with Bernal stacking, these connected sites are the ones where the density of states vanishes at zero energy, as discussed above. Hence, even in a clean system, the number of conducting channels in the direction perpendicular to the layers vanishes at zero energy (Nilsson et al., 2006a, 2007a). This situation is reminiscent of the in plane transport properties of a single layer graphene. Similar to the latter case, a self-consistent Born approximation for a small concentration of impurities leads to a finite conductivity, which becomes independent of the number of impurities.

# I. Transport near the Dirac point

In clean graphene, the number of channels available for electron transport decreases as the chemical potential approaches the Dirac energy. As a result, the conductance through a clean graphene ribbon is, at most, $4 e ^ { 2 } / h$ , where the factor of 4 stands for the spin and valley degeneracy. In addition, only one out of every three possible clean graphene ribbons have a conduction channel at the Dirac energy. The other two thirds are semiconducting, with a gap of the order of $v _ { \mathrm { { F } } } / W$ , where $W$ is the width. This result is a consequence of the additional periodicity introduced by the wavefunctions at the $K$ and $K ^ { \prime }$ points of the Brillouin Zone, irrespective of the boundary conditions.

A wide graphene ribbon allows for many channels, which can be approximately classified by the momentum perpendicular to the axis of the ribbon, $k _ { y }$ . At the Dirac energy, transport through these channels is inhibited by the existence of a gap, $\Delta _ { k _ { y } } = v _ { \mathrm { { F } } } k _ { y }$ . Transport through these channels is suppressed by a factor of the order of $e ^ { - k _ { y } L }$ , where $L$ is the length of the ribbon. The number of transverse channels increases as $W / a$ , where $W$ is the width of the ribbon and $a$ is a length of the order of the lattice spacing. The allowed values of $k _ { y }$ are $\propto n _ { y } / W$ , where $n _ { y }$ is an integer. Hence, for a ribbon such that $W \gg L$ , there are many channels which satisfy $k _ { y } L \ll 1$ . Transport through these channels is not strongly inhibited, and their contribution dominates when the Fermi energy lies near the Dirac point. The conductance arising from these channels is given approximately by (Katsnelson, 2006b; Tworzydlo et al., 2006):

The transmission at normal incidence, $k _ { y } = 0$ , is one, in agreement with the absence of backscattering in graphene, for any barrier that does not induce intervalley scattering (Katsnelson et al., 2006). The transmission of a given channel scales as $T ( k _ { y } ) = 1 / \cosh ^ { 2 } ( k _ { y } L / 2 )$ .

Eq.(179) shows that the contribution from all transverse channels lead to a conductance which scales, similar to a function of the length and width of the system, as the conductivity of a diffusive metal. Moreover, the value of the effective conductivity is of the order of $e ^ { 2 } / h$ . It can also be shown that the shot noise depends on current in the same way as in a diffusive metal. A detailed analysis of possible boundary conditions at the contacts and their influence on evanescent waves can be found in (Robinson and Schomerus, 2007; Schomerus, 2007). The calculations leading to eq.(179) can be extended to a graphene bilayer. The conductance is, again, a summation of terms arising from evanescent waves between the two contacts, and it has the dependence on sample dimensions of a 2D conductivity of the order of $e ^ { 2 } / h$ (Snyman and Beenakker, 2007), although there is a prefactor twice bigger than the one in single layer graphene.

The calculation of the conductance of clean graphene in terms of transmission coefficients, using the Landauer method leads to an effective conductivity which is equal to the value obtained for bulk graphene using diagrammatic methods, the Kubo formula (Peres et al., 2006d), in the limit of zero impurity concentration and zero doping. Moreover, this correspondence remains valid for the case of a bilayer without and with trigonal warping effects (Cserti et al., 2007a; Koshino and Ando, 2006).

Disorder at the Dirac energy changes the conductance of graphene ribbons in two opposite directions (Louis et al., 2007): i) a sufficiently strong disorder, with short range (intervalley) contributions, lead to a localized regime, where the conductance depends exponentially on the ribbon length, and ii) at the Dirac energy, disorder allows mid-gap states that can enhance the conductance mediated by evanescent waves discussed above. A fluctuating electrostatic potential also reduces the effective gap for the transverse channels, enhancing further the conductance. The resonant tunneling regime mediated by mid-gap state was suggested by analytical calculations (Titov, 2007). The enhancement of the conductance by potential fluctuations can also be studied semi-analytically. In the absence of intervalley scattering, it leads to an effective conductivity which grows with ribbon length (San-Jose et al., 2007). In fact, analytical and numerical studies (Bardarson et al., 2007; Lewenkopf et al., 2007; Nomura et al., 2007; San-Jose et al., 2007) show that the conductivity obeys a universal scaling with the lattice size $L$ :

$$
\sigma ( L ) = \frac { 2 e ^ { 2 } } { h } \left( A \ln ( L / \xi ) + B \right) ,
$$

$$
G \sim \frac { e ^ { 2 } } { h } \frac { W } { 2 \pi } \int d k _ { y } e ^ { - k _ { y } L } \sim \frac { e ^ { 2 } } { h } \frac { W } { L } .
$$

where $\xi$ is a length scale associated with range of interactions and $A$ and $B$ are numbers of the order of unit ( $A \approx 0 . 1 7$ and $B \approx 0 . 2 3$ for a graphene lattice in the shape of a square of size $L$ (Lewenkopf et al., 2007)). Notice, therefore, that the conductivity is always of the order of $e ^ { 2 } / h$ and has a weak dependence on size.

# J. Boltzmann Equation description of DC transport in doped graphene

It was shown experimentally that the DC conductivity of graphene depends linearly on the gate potential (Novoselov et al., 2005a, 2004, 2005b), except very close to the neutrality point (see Fig.30). Since the gate potential depends linearly on the electronic density, $n$ , one has a conductivity $\sigma \propto n$ . As shown by Shon and Ando (Shon and Ando, 1998) if the scatterers are short range the DC conductivity should be independent of the electronic density, at odds with the experimental result. It has been shown (Ando, 2006b; Nomura and MacDonald, 2006, 2007) that by considering a scattering mechanism based on screened charged impurities it is possible to obtain from a Boltzmann equation approach a conductivity varying linearly with the density, in agreement with the experimental result (Ando, 2006b; Katsnelson and Geim, 2008; Novikov, 2007b; Peres et al., 2007b; Trushin and Schliemann, 2007).

![](images/c9c3d02632fd19679c706158acf46634854b80cc8dc401f49820256e5c405af6.jpg)  
Figure 30 (Color online) An example of changes in conductivity $\sigma$ of graphene with varying gate voltage, $V _ { g }$ , and, therefore, carrier concentration $n$ Here $\sigma$ is proportional to $n$ as discussed in the text. Note that samples with higher mobility $\mathrm { ~ ( > ~ 1 ~ } \mathrm { m } ^ { 2 } / \mathrm { V s } \mathrm { ) }$ normally show a sublinear dependence, presumably indicating the presence of different types of scatterers. Inset: scanning-electron micrograph of one of experimental devices (in false colors matching those seen in visible optics. The scale of the micrograph is given by the width of the Hall bar, which is one micrometer. Adapted from (Novoselov et al., 2005a).

The Boltzmann equation has the form (Ziman, 1972)

$$
- v _ { k } \cdot \nabla _ { r } f ( \epsilon _ { k } ) - e ( E + v _ { k } \times H ) \cdot \nabla _ { k } f ( \epsilon _ { k } ) = - \left. \frac { \partial f _ { k } } { \partial t } \right| _ { s c a t t . } .
$$

The solution of the Boltzmann equation in its general form is difficult and one needs therefore to rely upon some approximation. The first step in the usual approximation scheme is to write the distribution as $f ( \epsilon _ { k } ) =$ $f ^ { 0 } ( \epsilon _ { k } ) + g ( \epsilon _ { k } )$ where $f ^ { 0 } ( \epsilon _ { k } )$ is the steady state distribution function and $g ( \epsilon _ { k } )$ is assumed to be small. Inserting this ansatz in (181) and keeping only terms that are linear in the external fields one obtains the linearized Boltzmann equation (Ziman, 1972) which reads

$$
\begin{array} { l } { - } & { \displaystyle \frac { \partial f ^ { 0 } \left( \epsilon _ { k } \right) } { \partial \epsilon _ { k } } v _ { k } \cdot \left[ \left( - \frac { \epsilon _ { k } - \zeta } { T } \right) \nabla _ { r } T + e \left( E - \frac { 1 } { e } \nabla _ { r } \zeta \right) \right] = } \\ { - } & { \displaystyle \frac { \partial f _ { k } } { \partial t } \bigg | _ { s c a t t . } + v _ { k } \cdot \nabla _ { r } g _ { k } + e ( v _ { k } \times H ) \cdot \nabla _ { k } g _ { k } . \qquad ( 1 8 2 } \end{array}
$$

The second approximation has to do with the form of the scattering term. The simplest approach is to introduce a relaxation time approximation:

$$
- \left. \frac { \partial f _ { k } } { \partial t } \right| _ { s c a t t . } \longrightarrow \frac { g _ { k } } { \tau _ { k } } ,
$$

where $^ { \prime } k$ is the relaxation time, assumed to be momentum dependent. This momentum dependence is determined phenomenologically in such way that the dependence of the conductivity upon the electronic density agrees with experimental data. The Boltzmann equation is certainly not valid at the Dirac point, but since many experiments are performed at finite carrier density, controlled by an external gate voltage, we expect the Boltzmann equation to give reliable results if an appropriate form for $\tau _ { k }$ is used (Adam et al., 2007).

Let us compute the Boltzmann relaxation time, $\tau _ { k }$ for two different scattering potentials:(i) a Dirac delta function potential; (ii) a unscreened Coulomb potential. The relaxation time $\tau _ { k }$ is defined as:

$$
{ \frac { 1 } { \tau _ { k } } } = n _ { i } \int d \theta \int { \frac { k ^ { \prime } d k ^ { \prime } } { ( 2 \pi ) ^ { 2 } } } S ( k , k ^ { \prime } ) ( 1 - \cos \theta ) ,
$$

where $n _ { i }$ is impurity concentration per unit of area, and the transition rate $S ( k , k ^ { \prime } )$ is given, in the Born approximation, by:

$$
S ( { \pmb k } , { \pmb k } ^ { \prime } ) = 2 \pi | H _ { { \pmb k } ^ { \prime } , { \pmb k } } | ^ { 2 } \frac { 1 } { v _ { F } } \delta ( { \pmb k } ^ { \prime } - { \pmb k } ) ,
$$

where the $v _ { F } k$ is the dispersion of Dirac fermions in graphene and $H _ { k ^ { \prime } , k }$ is defined as

$$
H _ { k ^ { \prime } , k } = \int d \pmb { r } \psi _ { k ^ { \prime } } ^ { * } ( \pmb { r } ) U _ { S } ( \pmb { r } ) \psi _ { k } ( \pmb { r } ) ,
$$

with $U _ { S } ( \mathbf { r } )$ the scattering potential and $\psi _ { k } ( r )$ is the electronic spinor wavefunction of a clean graphene sheet. If

the potential is short range,(Shon and Ando, 1998) of the form $U _ { S } = v _ { 0 } \delta ( r )$ , the Boltzmann relaxation time turns out to be

$$
\tau _ { k } = \frac { 4 v _ { F } } { n _ { i } v _ { 0 } ^ { 2 } } \frac { 1 } { k } .
$$

On the other hand, if the potential is the Coulomb potential, given by $U _ { S } ( \pmb { r } ) = e Q / ( 4 \pi \epsilon _ { 0 } \epsilon r )$ for charged impurities of charge $Q$ , the relaxation time is given by

$$
\tau _ { k } = \frac { v _ { F } } { u _ { 0 } ^ { 2 } } k .
$$

where $u _ { 0 } ^ { 2 } = n _ { i } Q ^ { 2 } e ^ { 2 } / ( 1 6 \epsilon _ { 0 } ^ { 2 } \epsilon ^ { 2 } )$ . As we argue below, the phenomenology of Dirac fermions implies that the scattering in graphene must be of the form (188).

Within the relaxation time approximation the solution of the linearized Boltzmann equation when an electric field is applied to the sample is

$$
g _ { \bf k } = - \frac { \partial f ^ { 0 } ( \epsilon _ { k } ) } { \partial \epsilon _ { k } } e \tau _ { k } v _ { k } \cdot E ,
$$

and the electric current reads (spin and valley indexes included)

$$
J = \frac { 4 } { A } \sum _ { k } e v _ { k } g _ { k } .
$$

Since at low temperatures the following relation $- \partial f ^ { \mathrm { 0 } } ( \epsilon _ { k } ) / \partial \epsilon _ { k } \to \delta ( \mu - v _ { F } k )$ holds, one can easily see that assuming (188) where $k$ is measured relatively to the Dirac point, the electronic conductivity turns out to be

$$
\sigma _ { x x } = 2 \frac { e ^ { 2 } } { h } \frac { \mu ^ { 2 } } { u _ { 0 } ^ { 2 } } = 2 \frac { e ^ { 2 } } { h } \frac { \pi v _ { F } ^ { 2 } } { u _ { 0 } ^ { 2 } } n ,
$$

where $u _ { 0 }$ is the strength of the scattering potential (with dimensions of energy). The electronic conductivity depends linearly on the electron density, in agreement with the experimental data. We stress that the Coulomb potential is one possible mechanism of producing a scattering rate of the form (188) but we do not exclude that other mechanisms may exist (see, for instance, (Katsnelson and Geim, 2008)).

# K. Magnetotransport and universal conductivity

The description of the magnetotransport properties of electrons in a disordered honeycomb lattice is complex because of the interference effects associated with the Hofstadter problem (Gumbs and Fekete, 1997). We shall simplify our problem by describing electrons in the honeycomb lattice as Dirac fermions in the continuum approximation, introduced in Sec. II.B. Furthermore, we will only focus on the problem of short range scattering in the unitary limit since in this regime many analytical results are obtained (Kumazaki and Hirashima,

2006; Mariani et al., 2007; Pereira et al., 2006, 2007a; Peres et al., 2006c; Skrypnyk and Loktev, 2006, 2007). The problem of magnetotransport in the presence of Coulomb impurities, as discussed in the previous section is still an open research problem. A similar approach was considered by Abrikosov in the quantum magnetoresistance study of non-stoichiometric chalcogenides (Abrikosov, 1998). In the case of graphene, the effective Hamiltonian describing Dirac fermions in a magnetic field (including disorder) can be written as: $H = H _ { 0 } + H _ { i }$ where $H _ { 0 }$ is given by (5) and $H _ { i }$ is the impurity potential reading (Peres et al., 2006c):

$$
H _ { i } = V \sum _ { j = 1 } ^ { N _ { i } } \delta ( \boldsymbol { r } - \boldsymbol { r } _ { j } ) I
$$

The formulation of the problem in second quantization requires the solution of $H _ { 0 }$ , which was done in Section II.I. The field operators, close to the K point, are defined as (the spin index is omitted for simplicity):

$$
\begin{array} { l } { \displaystyle \Psi ( \pmb { r } ) = \sum _ { k } \frac { e ^ { i k x } } { \sqrt { L } } \left( { 0 \atop \phi _ { 0 } ( y ) } \right) c _ { k , - 1 } } \\ { + \sum _ { n , k , \alpha } \frac { e ^ { i k x } } { \sqrt { 2 L } } \left( { \phi _ { n } ( y - k l _ { B } ^ { 2 } ) \atop \phi _ { n + 1 } ( y - k l _ { B } ^ { 2 } ) } \right) c _ { k , n , \alpha } , } \end{array}
$$

where $c _ { k , n , \alpha }$ destroys an electron in band $\alpha = \pm 1$ , with energy level $n$ and guiding center $k l _ { B } ^ { 2 }$ ; $c _ { k , - 1 }$ destroys an electron in the zero Landau level; the cyclotron frequency is given by (96). The sum over $n = 0 , 1 , 2 , \ldots$ , is cut off at $n _ { 0 }$ given by $E ( 1 , n _ { 0 } ) = W$ , where $W$ is of the order of the electronic bandwidth. In this representation $H _ { 0 }$ becomes diagonal, leading to Green's functions of the form (in Matsubara representation):

$$
G _ { 0 } ( k , n , \alpha ; i \omega ) = \frac { 1 } { i \omega - E ( \alpha , n ) } ,
$$

where

$$
E ( \alpha , n ) = \alpha \omega _ { c } \sqrt { n }
$$

are the Landau levels for this problem ( $\alpha = \pm 1$ labels the two bands). Notice that $G _ { 0 } ( k , n , \alpha ; i \omega )$ is effectively $k$ . independent, and $E ( \alpha , - 1 ) = 0$ is the zero energy Landau level. When expressed in the Landau basis, the scattering Hamiltonian (192) connects Landau levels of negative and positive energy.

# 1The full self-consistent Born approximation (FSBA)

In order to describe the effect of impurity scattering on the magnetoresistance of graphene, the Green's function for Landau levels in the presence of disorder needs to be computed. In the context of the 2D electron gas, an equivalent study was performed by Ohta and Ando, and $g _ { c } = A _ { c } / ( 2 \pi l _ { B } ^ { 2 } )$ is the degeneracy of a Landau level per unit cell. One should notice that the Green's functions do not depend upon $p$ explicitly. The self-consistent solution of Eqs. (197), (198), (199), (200) and (201) gives the density of states, the electron self-energy, and the change of Landau level energy position due to disorder.

![](images/a6484d74820752ddf1befe2d0364469de0a06df43adc421ad1d5629502ca2b3e.jpg)  
Figure 31 (Color online) Top: Electronic density of states (DOS), $\rho ( \omega )$ , as a function of $\omega / \omega _ { c }$ $\omega _ { c } = 0 . 1 4 ~ \mathrm { e V }$ ) in a magnetic field $B = 1 2$ T for different impurity concentrations $n _ { i }$ . Bottom: $\rho ( \omega )$ , as a function of $\omega / \omega _ { c }$ $\omega _ { c } = 0 . 1 \ \mathrm { e V }$ is the cyclotron frequency) in a magnetic field $B = 6$ T. The solid line shows the DOS in the absence of disorder. The position of the Landau levels in the absence of disorder are shown as vertical lines. The two arrows in the top panel show the position of the renormalized Landau levels (see Fig.32) given by the solution of Eq. (202). Adapted from Peres et al., 2006c.

The effect of disorder in the density of states of Dirac fermions in a magnetic field is shown in Fig. 31. For reference we note that $E ( 1 , 1 ) = 0 . 1 4 ~ \mathrm { e V }$ , for $B = 1 4 \mathrm { ~ T ~ }$ , and $E ( 1 , 1 ) = 0 . 1 ~ \mathrm { e V }$ , for $B = 6$ T. From Fig. 31 we see that, for a given $n _ { i }$ , the effect of broadening due to impurities is less effective as the magnetic field increases. It is also clear that the position of Landau levels is renormalized relatively to the non-disordered case. The renormalization of the Landau level position can be determined from poles of (197) and (198):

(Ando, 1974a,b,c, 1975; Ando and Uemura, 1974; Ohta, 1968, 1971) using the averaging procedure over impurity positions of Duke (Duke, 1968). Below the averaging procedure over impurity positions is performed in the standard way, namely, having determined the Green's function for a given impurity configuration $( r _ { 1 } , \dots r _ { N _ { i } } )$ , the position averaged Green's function is determined from:

$$
\begin{array} { r l r } & { } & { \langle G ( p , n , \alpha ; i \omega ; r _ { 1 } , \dots r _ { N _ { i } } ) \rangle \equiv G ( p , n , \alpha ; i \omega ) } \\ & { } & { = { \cal L } ^ { - 2 N _ { i } } \left[ \prod _ { j = 1 } ^ { N _ { i } } \int d { \pmb r } _ { j } \right] G ( p , n , \alpha ; i \omega ; { \pmb r } _ { 1 } , \dots { \pmb r } _ { N _ { i } } ) . } \end{array}
$$

In the presence of Landau levels the average over impurity positions involves the wavefunctions of the onedimensional harmonic oscillator. After a lengthy algebra, the Green's function in the presence of vacancies, in the FSBA, can be written as:

$$
\begin{array} { r l } & { G ( p , n , \alpha ; \omega + 0 ^ { + } ) \ = \ [ \omega - E ( n , \alpha ) - \Sigma _ { 1 } ( \omega ) ] ^ { - 1 } } \\ & { \ G ( p , - 1 ; \omega + 0 ^ { + } ) \ = \ [ \omega - \Sigma _ { 2 } ( \omega ) ] ^ { - 1 } , } \end{array}
$$

where

$$
\begin{array} { r c l } { { } } & { { } } & { { \Sigma _ { 1 } ( \omega ) ~ = ~ - n _ { i } [ Z ( \omega ) ] ^ { - 1 } , } } \\ { { } } & { { } } & { { \Sigma _ { 2 } ( \omega ) ~ = ~ - n _ { i } [ g _ { c } G ( p , - 1 ; \omega + 0 ^ { + } ) / 2 + Z ( \omega ) ] ^ { - 1 } , } } \\ { { } } & { { } } & { { Z ( \omega ) ~ = ~ g _ { c } G ( p , - 1 ; \omega + 0 ^ { + } ) / 2 } } \\ { { } } & { { } } & { { ~ + ~ g _ { c } \displaystyle \sum _ { n , \alpha } G ( p , n , \alpha ; \omega + 0 ^ { + } ) / 2 , } } \end{array}
$$

$$
\omega - E ( \alpha , n ) - \operatorname { R e } \Sigma ( \omega ) = 0 .
$$

Of course, due to the importance of scattering at low energies, the solution to Eq. (202) does not represent exact eigenstates of system since the imaginary part of the self-energy is non-vanishing. However, these energies do determine the form of the density of states, as we discuss below.

In Fig. 32, the graphical solution to Eq. (202) is given for two different energies $\boldsymbol { \mathscr { E } } ( - 1 , n )$ , with $n = 1 , 2$ ), its is clear that the renormalization is important for the first Landau level. This result is due to the increase in scattering at low energies, which is present already in the case of zero magnetic field. The values of $\omega$ satisfying Eq. (202) show up in the density of states as the energy values where the oscillations due to the Landau level quantization have a maximum. In Fig. 31, the position of the renormalized Landau levels is shown in the upper panel (marked by two arrows), corresponding to the bare energies $E ( - 1 , n )$ , with $n = 1 , 2$ . The importance of this renormalization decreases with the reduction of the number of impurities. This is clear in Fig. 31 where a visible shift toward low energies is evident when $n _ { i }$ has a small $1 0 \%$ change, from $n _ { i } = 1 0 ^ { - 3 }$ to $n _ { i } = 9 \times 1 0 ^ { - 4 }$ .

The study of the magnetoresistance properties of the system requires the calculation of the conductivity tensor. We compute the current-current correlation function and from it the conductivity tensor is derived. The details of the calculations are presented in (Peres et al., 2006c). If we however neglect the real part of the selfenergy, assume for $\mathrm { I m } \Sigma _ { i } ( \omega ) = \Gamma$ $( i = 1 , 2$ ) a constant value, and consider that $E ( 1 , 1 ) \gg \Gamma$ , these results reduce to those of (Gorbar et al., 2002).

It is instructive to consider first the case $\omega , T  0$ , which leads to $( \sigma _ { x x } ( 0 , 0 ) = \sigma _ { 0 }$ ):

$$
\begin{array} { r l r } { \sigma _ { 0 } } & { = } & { \frac { e ^ { 2 } } { h } \frac { 4 } { \pi } \left[ \frac { \mathrm { I m } \Sigma _ { 1 } ( 0 ) / \mathrm { I m } \Sigma _ { 2 } ( 0 ) - 1 } { 1 + ( \mathrm { I m } \Sigma _ { 1 } ( 0 ) / \omega _ { c } ) ^ { 2 } } \right. } \\ & { + } & { \left. \frac { n _ { 0 } + 1 } { n _ { 0 } + 1 + ( \mathrm { I m } \Sigma _ { 1 } ( 0 ) / \omega _ { c } ) ^ { 2 } } \right] , } \end{array}
$$

where we include a factor 2 due to the valley degeneracy. In the absence of a magnetic field ( $\omega _ { c }  0$ ) the above

![](images/d7859c7730bc6fd110a0689711b709c955e9b705223dc40ad9ed3dc034da86f0.jpg)  
Figure 32 (Color online) Imaginary (right) and real (left) parts of $\Sigma _ { 1 } ( \omega )$ (top) and $\Sigma _ { 2 } ( \omega )$ (bottom), in units of $\omega _ { c }$ ,as a function of $\omega / \omega _ { c }$ . The right panels also show the intercept of $\omega - E ( \alpha , n )$ with $\operatorname { R e } \Sigma ( \omega )$ as required by Eq. (202). Adapted from Peres et al., 2006c.

expression reduces to:

$$
\sigma _ { 0 } = { \frac { e ^ { 2 } } { h } } { \frac { 4 } { \pi } } \left[ 1 - { \frac { [ \mathrm { I m } \Sigma _ { 1 } ( 0 ) ] ^ { 2 } } { ( v _ { F } \Lambda ) ^ { 2 } + [ \mathrm { I m } \Sigma _ { 1 } ( 0 ) ] ^ { 2 } } } \right] ,
$$

where we have introduced the energy cut-off, $v _ { F } \Lambda$ . Either when $\mathrm { I m } \Sigma _ { 1 } ( 0 ) \simeq \mathrm { I m } \Sigma _ { 2 } ( 0 )$ and $\omega _ { c } \gg \mathrm { I m } \Sigma _ { 1 } ( 0 )$ (or $n _ { 0 } \gg$ $\mathrm { I m } \Sigma _ { 1 } ( 0 ) / \omega _ { c }$ , $\omega _ { c } = E ( 0 , 1 ) = \sqrt { 2 } v _ { F } / l _ { B } ^ { 2 }$ , or when $\Lambda v _ { F } \gg$ $\mathrm { I m } \Sigma _ { 1 } ( 0 )$ , in the absence of an applied field, Eqs. (203) and (204) reduce to:

$$
\sigma _ { 0 } = \frac { 4 } { \pi } \frac { e ^ { 2 } } { h } ,
$$

which is the so-called universal conductivity of graphene (Fradkin, 1986a,b; Katsnelson, 2006b; Lee, 1993; Ludwig et al., 1994; Nersesyan et al., 1994; Peres et al., 2006c; Tworzydlo et al., 2006; Yang and Nayak, 2002; Ziegler, 1998). This result was obtained previously by Ando and collaborators using the second order self-consistent Born approximation (Ando et al., 2002; Shon and Ando, 1998).

Because the DC magnetotransport properties of graphene are normally measured with the possibility of tuning its electronic density by a gate potential (Novoselov et al., 2004), it is important to compute the conductivity kernel, since this has direct experimental relevance. In the the case $\omega \to 0$ we write the conductivity $\sigma _ { x x } ( 0 , T )$ as:

$$
\sigma _ { x x } ( 0 , T ) = \frac { e ^ { 2 } } { \pi h } \int _ { - \infty } ^ { \infty } d \epsilon \frac { \partial f ( \epsilon ) } { \partial \epsilon } K _ { B } ( \epsilon ) ,
$$

where the conductivity kernel $K _ { B } ( \epsilon )$ is given in the Appendix of Ref. (Peres et al., 2006c). The magnetic field dependence of kernel $K _ { B } ( \epsilon )$ is shown in Fig. 33. Observe that the effect of disorder is the creation of a region where $K _ { B } ( \epsilon )$ remains constant before it starts to increase in energy with superimposed oscillations coming from the Landau levels. The same effect, but with the absence of the oscillations, was identified at the level of the self-consistent density of states plotted in Fig. 31. Together with $\sigma _ { x x } ( 0 , T )$ , the Hall conductivity $\sigma _ { x y } ( 0 , T )$ allows the calculation of the resistivity tensor (109).

![](images/6c6d7722771ea8a123b72bc39fb3f22a6944108621fccc9f5a53d2898611d296.jpg)  
Figure 33 (Color online) Conductivity kernel, $K ( \omega )$ (in units of $e ^ { 2 } / ( \pi h ) )$ , as a function of energy $\omega$ for different magnetic fields and for $n _ { i } = 1 0 ^ { - 3 }$ . The horizontal lines mark the universal limit of the conductivity per cone, $\sigma _ { 0 } = 2 e ^ { 2 } / ( \pi h )$ . The vertical lines show the position of the Landau levels in the absence of disorder. Adapted from Peres et al., 2006c.

Let us now focus on the optical conductivity, $\sigma _ { x x } ( \omega )$ (Gusynin et al., 2007; Peres et al., 2006c). This quantity can be probed by reflectivity experiments in the subterahertz to mid-infrared frequency range (Bliokh, 2005). This quantity is depicted in Fig. 34 for different magnetic fields. It is clear that the first peak is controlled by the $E ( 1 , 1 ) - E ( 1 , - 1 )$ , and we have checked that it does not obey any particular scaling form as a function of $\omega / B$ . On the other hand, as the effect of scattering becomes less important the high energy conductivity oscillations start obeying the scaling $\omega / \sqrt { B }$ , as we show in the lower right panel of Fig. 34.

# V. MANY-BODY EFFECTS

# A. Electron-phonon interactions

In Sec. IV.F.1 and Sec. IV.F.2 we discussed how static deformations of the graphene sheet due to bending and strain couple to the Dirac fermions via vector potentials. Just as bending has to do with the flexural modes of the graphene sheet (as discussed in Sec. III), strain fields are related to optical and acoustic modes (Wirtz and Rubio, 2004). Given the local displacements of the atoms in each sublattice, $\mathbf { u } _ { A }$ and $\mathbf { u } _ { B }$ , the electron-phonon coupling has essentially the form discussed previously for static fields.

The coupling to acoustic modes is the most straightforward one, since it already appears in the elastic theory. If $\mathbf { u } _ { \mathrm { a c } }$ is the acoustic phonon displacement, then the relation between this displacement and the atom displacement is given by equation (171), and its coupling to electrons is given by the vector potential (175) in the Dirac equation (150).

![](images/48d74209b3fa96b52c32d1042f5a36e21a08f41430296801149c9e3c9e988203.jpg)  
Figure 34 (Color online) Frequency dependent conductivity per cone, $\sigma ( \omega )$ (in units of $e ^ { 2 } / ( \pi h )$ at $T = 1 0$ K and $n _ { i } ~ = ~ 1 0 ^ { - 3 }$ , as a function of the energy $\omega$ (in units of $\omega _ { c }$ for different values of the magnetic field $B$ . The vertical arrows in the upper left panel, labeled a, b, and $\mathbf { c }$ , show the positions of the transitions between different Landau levels: $E ( 1 , 1 ) - E ( - 1 , 0 )$ , $E ( 2 , 1 ) - E ( - 1 , 0 )$ , and $E ( 1 , 1 ) - E ( 1 , - 1 )$ , respectively. The horizontal continuous lines show the value of the universal conductivity. The lower right panel shows the conductivity for different values of magnetic field as a function of $\omega / \sqrt { B }$ . Adapted from Peres et al., 2006c.

For optical modes the situation is slightly different since the optical mode displacement is (Ando, 2006a, 2007b):

$$
{ \bf u } _ { \mathrm { o p } } = \frac { 1 } { \sqrt { 2 } } ( { \bf u } _ { A } - { \bf u } _ { B } ) ,
$$

that is, the bond length deformation vector. To calculate the coupling to the electrons we can proceed as previously and compute the change in the nearest neighbor hopping energy due to the lattice distortion through (172), (173), (170), and (207). Once again the electron-phonon interaction becomes a problem of the coupling of the electrons with a vector potential as in (150) where the components of the vector potential are:

$$
\begin{array} { r c l } { { \displaystyle { \mathcal A } _ { x } ^ { ( \mathrm { o p } ) } ~ = ~ - \sqrt { \frac { 3 } { 2 } } \frac { \beta } { a ^ { 2 } } u _ { y } ^ { \mathrm { o p } } , } } \\ { { \displaystyle { \mathcal A } _ { y } ^ { ( \mathrm { o p } ) } ~ = ~ - \sqrt { \frac { 3 } { 2 } } \frac { \beta } { a ^ { 2 } } u _ { x } ^ { \mathrm { o p } } , } } \end{array}
$$

where $\beta = \partial t / \partial \ln ( a )$ was defined in (174). Notice that we can write: $\bar { \mathcal { A } } ^ { \mathrm { o p } } = - \sqrt { 3 / 2 } ( \beta / a ^ { 2 } ) \vec { \sigma } \times \mathbf { u } _ { \mathrm { o p } }$ . A similar expression is valid close to the $\mathrm { K } '$ point with $\vec { A }$ replaced by $- \bar { A }$ .

Optical phonons are particularly important in graphene research because of Raman spectroscopy. The latter has played a particularly important role in the study of carbon nanotubes (Saito et al., 1998) because of the 1D character of these systems, namely, the presence of van Hove singularities in the 1D spectrum lead to colossal enhancements of the Raman signal that can be easily detected, even for a single isolated carbon nanotube. In graphene the situation is rather different since its 2D character leads to a much smoother density of states (except for the van Hove singularity at high energies of the order of the hopping energy t ≈ 2.8 eV). Nevertheless, graphene is an open surface and hence is readily accessible by Raman spectroscopy. In fact, it has played a very important role because it allows the identification of the number of planes (Ferrari et al., 2006; Graf et al., 2007; Gupta et al., 2006; Pisana et al., 2007; Yan et al., 2007), and the study of the optical phonon modes in graphene, particularly the ones in the center of the BZ with momentum $q \approx 0$ . Similar studies have been performed in graphite ribbons (Cancado et al., 2004).

![](images/c2099f2976f1a0723abf83399fdc1f0863b2832d604bfd4d3f11c1ffaf55ce92.jpg)  
Figure 35 (Color online) Top: The continuous line is the relative phonon frequency shift as a function of $\mu / \omega _ { 0 }$ , and the dashed line is the damping of the phonon due to electron-hole pair creation; Bottom: (a) Electron-hole process that leads to phonon softening $\omega _ { 0 } > 2 \mu$ ), and (b) electron-hole process that leads to phonon hardening $\omega _ { 0 } < 2 \mu$ .

Let us consider the effect of the Dirac fermions on the optical modes. If one treats the vector potential, electron-phonon coupling, (150) and (208) up to second order perturbation theory, its main effect is the polarization of the electron system by creating electron-hole pairs. In the QED language, the creation electron-hole pairs is called pair (electron/anti-electron) production (Castro Neto, 2007). Pair production is equivalent to a renormalization of the phonon propagator by a selfenergy that is proportional to the polarization function of the Dirac fermions.

The renormalized phonon frequency, $\Omega _ { 0 } ( \mathbf { q } )$ , is given by (Ando, 2006a, 2007b; Castro Neto and Guinea, 2007; Lazzeri and Mauri, 2006; Saha et al., 2007):

$$
\Omega _ { 0 } ( { \bf q } ) \approx \omega _ { 0 } - \frac { 2 \beta ^ { 2 } } { a ^ { 2 } \omega _ { 0 } } \chi ( { \bf q } , \omega _ { 0 } ) ,
$$

where $\omega _ { 0 }$ is the bare phonon frequency, and the electronphonon polarization function is given by:

$$
\chi ( \mathbf { q } , \omega ) = \sum _ { s , s ^ { \prime } = \pm 1 } \int \frac { d ^ { 2 } \mathbf { k } \quad f [ E _ { s } ( \mathbf { k } + \mathbf { q } ) ] - f [ E _ { s ^ { \prime } } ( \mathbf { k } ) ] } { ( 2 \pi ) ^ { 2 } \omega _ { 0 } - E _ { s } ( \mathbf { k } + \mathbf { q } ) + E _ { s ^ { \prime } } ( \mathbf { k } ) + i \eta } ( \frac { 4 } { 3 { \pi } } \mathbf { \epsilon }
$$

where $E _ { s } ( \mathbf { q } )$ is the Dirac fermion dispersion ( $s = + 1$ for the upper band, and $s = - 1$ for the lower band), and $f [ E ]$ is the Fermi-Dirac distribution function. For Raman spectroscopy, the response of interest is at $q = 0$ where clearly only the interband processes such that $s s ^ { \prime } = - 1$ (that is, processes between the lower and upper cones) contribute. The electron-phonon polarization function can be easily calculated using the linearized Dirac fermion dispersion (7) and the low energy density of states (15):

$$
\begin{array} { l } { { \displaystyle \chi ( 0 , \omega _ { 0 } ) = \frac { 6 \sqrt { 3 } } { \pi v _ { F } ^ { 2 } } \int _ { 0 } ^ { v _ { F } \Lambda } d E E \left( f [ - E ] - f [ E ] \right) \left( \frac { 1 } { \omega _ { 0 } + 2 E + i \eta } \right. } } \\ { { \displaystyle - \left. \frac { 1 } { \omega - 2 E + i \eta } \right) , } } \end{array}
$$

where we have introduced the cut-off momentum $\Lambda$ $\approx$ $1 / a$ ) so that the integral converges in the ultraviolet. At zero temperature, $T = 0$ , we have $f [ E ] = \theta ( \mu - E )$ and we assume electron doping, $\mu > 0$ , so that $f [ - E ] = 1$ (for the case of hole doping, $\mu < 0$ , is obtained by electronhole symmetry). The integration in (211) gives:

$$
\begin{array} { r l r } { \chi ( 0 , \omega _ { 0 } ) } & { = } & { \frac { 6 \sqrt { 3 } } { \pi v _ { F } ^ { 2 } } \left[ v _ { F } \Lambda - \mu + \frac { \omega _ { 0 } } { 4 } \left( \ln \left| \frac { \omega _ { 0 } / 2 + \mu } { \omega _ { 0 } / 2 - \mu } \right| \right. \right. } \\ { ~ } & { + } & { \left. \left. i \pi \theta \left( \omega _ { 0 } / 2 - \mu \right) \right) \right] , } \end{array}
$$

where the cut-off dependent term is a contribution coming from the occupied states in the lower $\pi$ band and hence is independent of the value of the chemical potential. This contribution can be fully incorporated into the bare value of $\omega _ { 0 }$ in (209). Hence the relative shift in the phonon frequency can be written as:

$$
\frac { \delta \omega _ { 0 } } { \omega _ { 0 } } \approx - \frac { \lambda } { 4 } \left( - \frac { \mu } { \omega _ { 0 } } + \ln \left| \frac { \omega _ { 0 } / 2 + \mu } { \omega _ { 0 } / 2 - \mu } \right| + i \pi \theta \left( \omega _ { 0 } / 2 - \mu \right) \right) \Omega
$$

where

$$
\lambda = \frac { 3 6 \sqrt { 3 } } { \pi } \frac { \beta ^ { 2 } } { 8 M a ^ { 2 } \omega _ { 0 } } ,
$$

is the dimensionless electron-phonon coupling. Notice that (213) has a real and imaginary part. The real part represents the actual shift in frequency, while the imaginary part gives the damping of the phonon mode due to pair production (see Fig. 35). There is a clear change in behavior depending whether $\mu$ is larger or smaller than $\omega _ { 0 } / 2$ . For $\mu ~ < ~ \omega _ { 0 } / 2$ there is a decrease in the phonon frequency implying that the lattice is softening, while for $\mu > \omega _ { 0 } / 2$ the lattice hardens. The interpretation for this effect is also given in Fig. 35. On the one hand, if the frequency of the phonon is larger than twice the chemical potential, real electron-hole pairs are produced, leading to stronger screening of the electronion interaction and hence, to a softer phonon mode. At the same time the phonons become damped and decay. On the other hand, if the frequency of the phonon is smaller than the twice the chemical potential, the production of electron-hole pairs is halted by the Pauli principle and only virtual excitations can be generated leading to polarization and lattice hardening. In this case, there is no damping and the phonon is long lived. This amazing result has been observed experimentally by Raman spectroscopy (Pisana et al., 2007; Yan et al., 2007). Electron-phonon coupling has also been investigated theoretically in the case of a finite magnetic field (Ando, 2007a; Goerbig et al., 2007). In this case, resonant coupling occurs due to the large degeneracy of the Landau levels and different Raman transitions are expected as compared with the zero-field case. The coupling of electrons to flexural modes on a free standing graphene sheet was discussed in ref. (Mariani and von Oppen, 2007).

# B. Electron-electron interactions

Of all disciplines of condensed matter physics, the study of electron-electron interactions is probably one of the most complex since it involves the understanding of the behavior of a macroscopic number of variables. Hence, the problem of interacting systems is a field in constant motion and we shall not try to give here a comprehensive survey of the problem for graphene. Instead, we will focus on a small number of topics that are of current discussion in the literature.

Since graphene is a truly 2D system, it is informative to compare it with the more standard 2DEG that has been studied extensively in the last 25 years since the development of heterostructures and the discovery of the quantum Hall effect (for a review, see (Stone, 1992)). At the simplest level, metallic systems have two main kind of excitations: electron-hole pairs and collective modes such as plasmons.

Electron-hole pairs are incoherent excitations of the Fermi sea and a direct result of Pauli's exclusion principle: an electron inside the Fermi sea with momentum $\mathbf { k }$ is excited outside the Fermi sea to a new state with momentum $\mathbf { k } + \mathbf { q }$ , leaving a hole behind. The energy associated with such an excitation is simply: $\omega = \epsilon _ { \mathbf { k } + \mathbf { q } } - \epsilon _ { \mathbf { k } }$ and for states close to the Fermi surface ${ \bf k } \approx { \bf k } _ { F }$ )their energy scales linearly with the excitation momentum, $\omega _ { q } \approx v _ { F } q$ .

![](images/5fe0ed0bfe363677834de57b1d3b68373acb98e09a1243401510d875bae8cf90.jpg)  
Figure 36 (Color online) Electron-hole continuum and collective modes of: (a) a 2DEG; (b) undoped graphene; (c) doped graphene.

graphene reflect its screening properties as well. In fact, the polarization and dielectric functions of undoped graphene are rather different from the ones of the 2DEG (Lindhard function). In the random phase approximation (RPA), the polarization function can be calculated analytically (González et al., 1993a, 1994; Shung, 1986a):

$$
\Pi ( q , \omega ) = \frac { q ^ { 2 } } { 4 \sqrt { v _ { F } ^ { 2 } q ^ { 2 } - \omega ^ { 2 } } } ,
$$

and hence, for $\omega > v _ { F } q$ the polarization function is imaginary indicating the damping of electron-hole pairs. Notice that the static polarization function ( $\omega = 0$ )vanishes linearly with $q$ , indicating the lack of screening in the system. This polarization function has been also calculated in the presence of a finite chemical potential (Ando, 2006b; Hwang and Das Sarma, 2007; Shung, 1986a,b; Wunsch et al., 2006).

Undoped, clean graphene is a semimetal, with a vanishing density of states at the Fermi level. As a result the linear Fermi Thomas screening length diverges, and the long range Coulomb interaction is not screened. At finite electron density $n$ , the Thomas-Fermi screening length reads :

$$
\lambda _ { T F } \approx \frac { 1 } { 4 \alpha } \frac { 1 } { k _ { F } } = \frac { 1 } { 4 \alpha } \frac { 1 } { \sqrt { \pi n } } ,
$$

In a system with non-relativistic dispersion such as normal metals and semiconductors, the electron-hole continuum is made out of intra-band transitions only and exists even at zero energy since it is always possible to produce electron-hole pairs with arbitrarily low energy close to the Fermi surface, as shown in Fig. 36(a). Besides that, the 2DEG can also sustain collective excitations such as plasmons that have dispersion: $\omega _ { \mathrm { p l a s m o n } } ( q ) \propto \sqrt { q }$ , and exist outside the electron-hole continuum at sufficiently long wavelengths (Shung, 1986a).

In systems with relativistic-like dispersion, such as graphene, these excitations change considerably, especially when the Fermi energy is at the Dirac point. In this case the Fermi surface shrinks to a point and hence intraband excitations disappear and only interband transitions between the lower and upper cones can exist (see Fig.36(b)). Therefore, neutral graphene has no electronhole excitations at low energy, instead each electron-hole pair costs energy and hence the electron-hole occupies the upper part of the energy versus momentum diagram. In this case, plasmons are suppressed and no coherent collective excitations can exist. If the chemical potential is moved away from the Dirac point then intra-band excitations are restored and the electron-hole continuum of graphene shares features of the 2DEG and undoped graphene. The full electron-hole continuum of doped graphene is shown in Fig. 36(c), and in this case plasmon modes are allowed. As the chemical potential is raised away from the Dirac point, graphene resembles more and more the 2DEG.

These features in the elementary excitations of

where

$$
\alpha = \frac { e ^ { 2 } } { \epsilon _ { 0 } v _ { F } } ,
$$

is the dimensionless coupling constant in the problem (the analogue of (143) in the Coulomb impurity problem). Going beyond the linear ThomasFermi regime, it has been shown that the Coulomb law is modified (Fogler et al., 2007b; Katsnelson, 2006a; Zhang and Fogler, 2007).

The Dirac Hamiltonian in the presence of interactions can be written as:

$$
\begin{array} { l } { { \displaystyle { \mathcal { H } } \ \equiv \ - i v _ { \mathrm { F } } \int d ^ { 2 } { \bf r } \ \hat { \Psi } ^ { \dagger } ( { \bf r } ) \sigma \cdot \nabla \hat { \Psi } ( { \bf r } ) } \ ~ } \\ { { \displaystyle ~ + \frac { e ^ { 2 } } { 2 \epsilon _ { 0 } } \int d ^ { 2 } { \bf r } d ^ { 2 } { \bf r ^ { \prime } } \frac { 1 } { | { \bf r - r ^ { \prime } } | } \hat { \rho } ( { \bf r } ) \hat { \rho } ( { \bf r ^ { \prime } } ) , } \ ~ } \end{array}
$$

where

$$
\hat { \rho } ( { \bf r } ) = \hat { \Psi } ^ { \dagger } ( { \bf r } ) \hat { \Psi } ( { \bf r } ) ,
$$

is the electronic density. Observe that Coulomb interaction, unlike in QED, is assumed to be instantaneous since $v _ { F } / c ~ \approx ~ 1 / 3 0 0$ and hence retardation effects are very small. Moreover, the photons propagate in 3D space whereas the electrons are confined to the 2D graphene sheet. Hence, the Coulomb interaction breaks the Lorentz invariance of the problem and makes the many-body situation rather different from the one in QED (Baym and Chin, 1976). Furthermore, the problem depends on two parameters, $v _ { \mathrm { F } }$ and $e ^ { 2 } / \epsilon _ { 0 }$ . Under a dimensional scaling, $\mathbf { r }  \lambda \mathbf { r } , t  \lambda t , \Psi  \lambda ^ { - 1 } \Psi$ , both parameters remain invariant. In RG language, the Coulomb interaction is a marginal variable, whose strength relative to the kinetic energy does not change upon a change in scale. If the units are chosen in such a way that $v _ { \mathrm { F } }$ is dimensionless, the value of $e ^ { 2 } / \epsilon _ { 0 }$ will also be rendered dimensionless. This is the case in theories considered renormalizable in quantum field theory.

![](images/b9d25053d261f99f8f205dae92745cb34494eb6e54b117e68244293a15330cc0.jpg)  
Figure 37 Hartree-Fock self-energy diagram which leads to a logarithmic renormalization of the Fermi velocity.

The Fermi velocity in graphene is comparable to that in half-filled metals. In solids with lattice constant $a$ , the total kinetic energy per site, $1 / ( m a ^ { 2 } )$ , where $m$ is the bare mass of the electron, is of the same order of magnitude as the electrostatic energy, $e ^ { 2 } / ( \epsilon _ { 0 } a )$ . The Fermi velocity for fillings of the order of unity is $v _ { \mathrm { F } } \sim 1 / ( m a )$ . Hence, $e ^ { 2 } / ( \epsilon _ { 0 } v _ { \mathrm { F } } ) ~ \sim ~ 1$ . This estimate is also valid in graphene. Hence, unlike in QED, where $\alpha _ { \mathrm { Q E D } } = 1 / 1 3 7$ , the coupling constant in graphene is $\alpha \sim 1$ .

Despite the fact that the coupling constant is of the order of unity, a perturbative RG analysis can be applied. RG techniques allow us to identify stable fixed points of the model, which may be attractive over a broader range than the one where a perturbative treatment can be rigorously justified. Alternatively, an RG scheme can be reformulated as the process of piecewise integration of high energy excitations (Shankar, 1994). This procedure leads to changes in the effective low energy couplings. The scheme is valid if the energy of the renormalized modes is much larger than the scales of interest.

The Hartree-Fock correction due to Coulomb interactions between electrons (given by the diagram in Fig. 37) gives a logarithmic correction to the electron self-energy (González et al., 1994):

$$
\Sigma _ { H F } ( \mathbf { k } ) = \frac { \alpha } { 4 } k \ln \left( \frac { \Lambda } { k } \right)
$$

where $\Lambda$ is a momentum cutoff which sets the range of validity of the Dirac equation. This result remains true even to higher order in perturbation theory (Mishchenko, 2007) and is also obtained in large $N$ expansions (Rosenstein et al., 1989, 1991; Son, 2007) ( $N$ is the number of flavors of Dirac fermions), with the only modification being the prefactor in (220). This result implies that the Fermi velocity is renormalized towards higher values. As a consequence, the density of states near the Dirac energy is reduced, in agreement with the general trend of repulsive interactions to induce or increase gaps.

This result can be understood from the RG point of view by studying the effect of reducing the cut-off from $\Lambda$ to $\Lambda - d \Lambda$ and its effect on the effective coupling. It can be shown that $\alpha$ obeys the equation (González et al., 1994):

$$
\Lambda \frac { \partial \alpha } { \partial \Lambda } = - \frac { \alpha } { 4 } .
$$

Therefore, the Coulomb interaction becomes marginally irrelevant. These features are confirmed by a full relativistic calculation, although the Fermi velocity cannot, obviously, surpass the velocity of light (González et al., 1994). This result indicates that strongly correlated electronic phases, such as ferromagnetism (Peres et al., 2005) and Wigner crystals (Dahal et al., 2006) are suppressed in clean graphene.

A calculation of higher order self-energy terms leads to a wavefunction renormalization, and to a finite quasiparticle lifetime, which grows linearly with quasiparticle energy (González et al., 1994, 1996). The wavefunction renormalization implies that the quasiparticle weight tends to zero as its energy is reduced. A strong coupling expansion is also possible, assuming that the number of electronic flavors justifies an RPA expansion, keeping only electron-hole bubble diagrams (González et al., 1999). This analysis confirms that the Coulomb interaction is renormalized towards lower values.

The enhancement in the Fermi velocities leads to a widening of the electronic spectrum. This is consistent with measurements of the gaps in narrow single wall nanotubes, which show deviations from the scaling with $R ^ { - 1 }$ , where $R$ is the radius, expected from the Dirac equation (Kane and Mele, 2004). The linear dependence of the inverse quasiparticle lifetime with energy is consistent with photo-emission experiments in graphite, for energies larger with respect to the interlayer interactions (Bostwick et al., 2007b; Sugawara et al., 2007; Xu et al., 1996; Zhou et al., 2006a,c). Note that in graphite, band structure effects modify the lifetimes at low energies (Spataru et al., 2001). The vanishing of the quasiparticle peak at low energies can lead to an energy dependent renormalization of the interlayer hopping (Vozmediano et al., 2002, 2003). Other thermodynamic properties of undoped and doped graphene can also be calculated (Barlas et al., 2007; Vafek, 2007).

Non-perturbative calculations of the effects of the long range interactions in undoped graphene show that a transition to a gapped phase is also possible, when the number of electronic flavors is large (Khveshchenko, 2001; Khveshchenko and Shively, 2006; Luk'yanchuk and Kopelevich, 2004). The broken symmetry phase is similar to the excitonic transition found in materials where it becomes favorable to create electronhole pairs that then form bound excitons (excitonic transition).

Undoped graphene cannot have well defined plasmons, as their energies fall within the electron-hole continuum, and therefore have a significant Landau damping. At finite temperatures, however, thermally excited quasiparticles screen the Coulomb interaction, and an acoustic collective charge excitation can exist (Vafek, 2006).

Doped graphene shows a finite density of states at the Fermi level, and the long range Coulomb interaction is screened. Accordingly, there are collective plasma interactions near ${ \bf q } \to 0$ , which disperse as $\omega _ { p } \sim { \sqrt { | \mathbf { q } | } }$ , since the system is 2D (Campagnoli and Tosatti, 1989; Shung, 1986a,b). The fact that the electronic states are described by the massless Dirac equation implies that $\omega _ { P } \propto$ $n ^ { 1 / 4 }$ , where $n$ is the carrir density. The static dielecic constant has a continuous derivative at $2 k _ { \mathrm { F } }$ , unlike in the case of the 2D electron gas (Ando, 2006b; Sarma et al., 2007; Wunsch et al., 2006). This fact is associated with the suppressed backward scattering in graphene. The simplicity of the band structure of graphene allows analytical calculation of the energy and momentum dependence of the dielectric function (Sarma et al., 2007; Wunsch et al., 2006). The screening of the long-range Coulomb interaction implies that the low energy quasiparticles show a quadratic dependence on energy with respect to the Fermi energy (Hwang et al., 2007).

One way to probe the strength of the electronelectron interactions is via the electronic compressibility. Measurements of the compressibility using a single electron transistor (SET) show very little sign of interactions in the system, being well fitted by the non-interacting result that, contrary to the two-dimensional electron gas (2DEG) (Eisenstein et al., 1994; Giuliani and Vignale, 2005), is positively divergent (Martin et al., 2007; Polini et al., 2007). Bilayer graphene, on the other hand, shares characteristics of the single layer and the 2DEG with a non-monotonic dependence of the compressibility on the carrier density (Kusminskiy et al., 2007). In fact, bilayer graphene very close to half-filling has been predicted to be unstable towards Wigner crystallization (Dahal et al., 2007), just like the 2DEG. Furthermore, according to Hartree-Fock calculations, clean bilayer graphene is unstable towards ferromagnetism (Nilsson et al., 2006b).

# 1. Screening in graphene stacks

The electron-electron interaction leads to the screening of external potentials. In a doped stack, the charge tends to accumulate near the surfaces, and its distribution is determined by the dielectric function of the stack in the out-of-plane direction. The same polarizability describes the screening of an external field perpendicular to the layers, like the one induced by a gate in electrically doped systems (Novoselov et al., 2004). The selfconsistent distribution of charge in a biased graphene bilayer has been studied in ref. (McCann, 2006). From the observed charge distribution and self-consistent calculations, an estimate of the band structure parameters and their relation with the induced gap can be obtained (Castro et al., 2007a).

In the absence of interlayer hopping, the polarizability of a set of stacks of 2D electron gases can be written as a sum of the screening by the individual layers. Using the accepted values for the effective masses and carrier densities of graphene, this scheme gives a first approximation to screening in graphite (Visscher and Falicov, 1971). The screening length in the out of plane direction is of about 2 graphene layers (Morozov et al., 2005). This model is easily generalizable to a stack of semimetals described by the 2D Dirac equation (González et al., 2001). At half filling, the screening length in all directions diverges, and the screening effects are weak.

Interlayer hopping modifies this picture significantly. The hopping leads to coherence (Guinea, 2007). The out of plane electronic dispersion is similar to that of a one dimensional conductor. The out of plane polarizability of a multilayer contains intra- and interband contributions. The subbands in a system with the Bernal stacking have a parabolic dispersion, when only the nearest neighbor hopping terms are included. This band structure leads to an interband susceptibility described by a sum of terms like those in (228), which diverges at halffilling. In an infinite system, this divergence is more pronounced at $k _ { \perp } = \pi / c$ , that is, for a wave vector equal to twice the distance between layers. This effect greatly enhances Friedel like oscillations in the charge distribution in the out of plane direction, which can lead to the changes in the sign of the charge in neighboring layers (Guinea, 2007). Away from half-filling a graphene bilayer behaves, from the point of view of screening, in a way very similar to the 2DEG (Wang and Chakraborty, 2007b).

# C. Short range interactions

In this section we discuss the effect of short range Coulomb interactions on the physics of graphene. The simplest carbon system with a hexagonal shape is the benzene molecule. The value of the Hubbard interaction among $\pi$ -electrons was, for this system, computed long ago by Parr et al. (Parr et al., 1950), yielding a value of $U = 1 6 . 9 3 \mathrm { e V }$ For comparison purposes, in polyacetylene the value for the Hubbard interaction is $U \simeq 1 0 \ \mathrm { e V }$ and the hopping energy is $t \approx 2 . 5$ eV (Baeriswyl et al., 1986). These two examples just show that the value of the onsite Coulomb interaction is fairly large for $\pi -$ electrons. As a first guess for graphene, one can take $U$ to be of the same order as for polyacethylene, with the hopping integral $t \simeq 2 . 8 ~ \mathrm { e V }$ . Of course in pure graphene the electronelectron interaction is not screened, since the density of states is zero at the Dirac point, and one should work out the effect of Coulomb interactions by considering the bare Coulomb potential. On the other hand, as we have seen before, defects induce a finite density of states at the Dirac point, which could lead to an effective screening of the long-range Coulomb interaction. Let us assume that the bare Coulomb interaction is screened in graphene and that Coulomb interactions are represented by the Hubbard interaction. This means that we must add to the Hamiltonian (5) a term of the form:

$$
\begin{array} { l } { { \displaystyle H _ { U } ~ = ~ U \sum _ { { \cal R } _ { i } } \left[ a _ { \uparrow } ^ { \dagger } ( { \bf R } _ { i } ) a _ { \uparrow } ( { \bf R } _ { i } ) a _ { \downarrow } ^ { \dagger } ( { \bf R } _ { i } ) a _ { \downarrow } ( { \bf R } _ { i } ) \right. } } \\ { { \displaystyle ~ + ~ \left. b _ { \uparrow } ^ { \dagger } ( { \bf R } _ { i } ) b _ { \uparrow } ( { \bf R } _ { i } ) b _ { \downarrow } ^ { \dagger } ( { \bf R } _ { i } ) b _ { \downarrow } ( { \bf R } _ { i } ) \right] } } \end{array}
$$

The simplest question one can ask is whether this system shows a tendency toward some kind of magnetic order driven by the interaction $U$ . Within the simplest HartreeFock approximation (Peres et al., 2004), the instability line toward ferromagnetism is given by:

$$
U _ { F } ( \mu ) = \frac { 2 } { \rho ( \mu ) } ,
$$

which is nothing but the Stoner criterion. Similar results are obtained in more sophisticated calculations (Herbut, 2006). Clearly, at half-filling the value for the density of states is $\rho ( 0 ) = 0$ and the critical value for $U _ { F }$ is arbitrarily large. Therefore we do not expect a ferromagnetic ground state at the neutrality point of one electron per carbon atom. For other electronic densities, $\rho ( \mu )$ becomes finite producing a finite value for $U _ { F }$ . We note that the inclusion of $t ^ { \prime }$ does not change these findings, since the density of states remains zero at the neutrality point.

The line toward an antiferromagnetic ground state is given by (Peres et al., 2004)

$$
U _ { A F } ( \mu ) = \frac { 2 } { \frac { 1 } { N } \sum _ { k , \mu > 0 } \frac { 1 } { | E _ { + } ( k ) | } } ,
$$

where $E _ { + } ( k )$ is given in (6). This result gives a finite $U _ { A F }$ at the neutrality point (Martelo et al., 1997; Sorella and Tosatti, 1992):

$$
U _ { A F } ( 0 ) = 2 . 2 3 t .
$$

Quantum Monte Carlo calculations (Paiva et al., 2005; Sorella and Tosatti, 1992), raise however its value to:

$$
U _ { A F } ( 0 ) \simeq 5 t .
$$

Taking for graphene the same value for $U$ as in polyacetylene and $t = 2 . 8 \mathrm { e V }$ , one obtains $U / t \simeq 3 . 6$ , which put the system far from the transition toward an antiferromagnet ground state. Yet another possibility is that the system may be in a sort of a quantum spin liquid (Lee and Lee, 2005) (as originally proposed by Pauling (Pauling, 1972) in 1956) since mean field calculations give a critical value for $U$ to be of the order of $U / t \simeq 1 . 7$ . Whether this type of ground state really exists and whether quantum fluctuations pushes this value of $U$ toward larger values is not known.

# Bilayer graphene: exchange

The exchange interaction can be large in an unbiased graphene bilayer with a small concentration of carriers. It was shown previously that the exchange contribution to the electronic energy of a single graphene layer does not lead to a ferromagnetic instability (Peres et al., 2005). The reason for this is a significant contribution from the interband exchange, which is a term usually neglected in doped semiconductors. This contribution depends on the overlap of the conduction and valence wavefunctions, and it is modified in a bilayer. The interband exchange energy is reduced in a bilayer (Nilsson et al., 2006c), and a positive contribution that depends logarithmically on the bandwidth in graphene is absent in its bilayer. As a result, the exchange energy becomes negative, and scales as $n ^ { 3 / 2 }$ , where $n$ is the carrier density, similar to the 2DEG. The quadratic dispersion at low energies implies that the kinetic energy scales as $n ^ { 2 }$ , again as in the 2DEG. This expansion leads to:

$$
E = E _ { k i n } + E _ { e x c } \approx \frac { \pi v _ { F } ^ { 2 } n ^ { 2 } } { 8 t _ { \perp } } - \frac { e ^ { 2 } n ^ { 3 / 2 } } { 2 7 \sqrt { \pi } \epsilon _ { 0 } }
$$

Writing $n _ { \uparrow } = ( n + s ) / 2 , n _ { \downarrow } = ( n - s ) / 2$ , where $s$ is the magnetization, (227) predicts a second order transition to a ferromagnetic state for $n = ( 4 e ^ { 4 } t ^ { 2 } ) / ( 8 1 \pi ^ { 3 } v _ { F } ^ { 4 } \epsilon _ { 0 } )$ . Higher order corrections to (227) lead to a first order transition at slightly higher densities (Nilsson et al., 2006c). For a ratio $\gamma _ { 1 } / \gamma _ { 0 } \approx 0 . 1$ , this analysis implies that a graphene bilayer should be ferromagnetic for carrier densities such that $| n | \lesssim 4 \times 1 0 ^ { 1 0 } \mathrm { c m ^ { - 2 } }$ .

A bilayer is also the unit cell of Bernal graphite, and the exchange instability can also be studied in an infinite system. Taking into account nearest neighbor interlayer hopping only, bulk graphite should also show an exchange instability at low doping. In fact, there is some experimental evidence for a ferromagnetic instability in strongly disordered graphite (Esquinazi et al., 2002, 2003; Kopelevich and Esquinazi, 2006).

The analysis described above can be extended to the biased bilayer, where a gap separates the conduction and valence bands (Stauber et al., 2007). The analysis of this case is somewhat different, as the Fermi surface at low doping is a ring, and the exchange interaction can change its bounds. The presence of a gap reduces further the mixing of the valence and conduction band, leading to an enhancement of the exchange instability. At all doping levels, where the Fermi surface is ring shaped, the biased bilayer is unstable towards ferromagnetism.

# 2.Bilayer graphene: short range interactions

The band structure of a graphene bilayer, at half filling, leads to logarithmic divergences in different response functions at $\mathbf q = 0$ . The two parabolic bands that are tangent at $\mathbf k = 0$ lead to a susceptibility which is propor

![](images/b72a42c14ca214c4f6e1871c9e3289ca9b0f155598495d848a9a6ea9a995d7d0.jpg)  
Figure 38 (Color online) Sketch of the expected magnetization of a graphene bilayer at half-filling.

tional to:

$$
\chi ( \vec { q } , \omega ) \propto \int _ { | \vec { q } | < \Lambda } d ^ { 2 } { \bf k } \frac { 1 } { \omega - ( v _ { F } ^ { 2 } / t ) | { \bf k } | ^ { 2 } } \propto \log \left( \frac { \Lambda } { \sqrt { ( \omega t ) / v _ { F } ^ { 2 } } } \right)
$$

where $\Lambda \sim \sqrt { t ^ { 2 } / v _ { F } ^ { 2 } }$ is a high momentum cutoff. These logarithmic divergences are similar to the ones which show up when the Fermi surface of a 2D metal is near a saddle point in the dispersion relation (González et al., 1996). A full treatment of these divergences requires a RG approach (Shankar, 1994). Within a simpler mean field treatment, however, it is easy to notice that the divergence of the bilayer susceptibility gives rise to an instability towards an antiferromagnetic phase, where the carbon atoms which are not connected to the neighboring layers acquire a finite magnetization, while the magnetization of the atoms with neighbors in the contiguous layers remain zero. A scheme of the expected ordered state is shown in Fig. 38.

# D. Interactions in high magnetic fields

The formation of Landau levels enhances the effect of interactions due to the quenching of the kinetic energy. This effect is most pronounced at low fillings, when only the lowest levels are occupied. New phases may appear at low temperatures. We consider here phases different from the fractional quantum Hall effect, which has not been observed in graphene so far. The existence of new phases can be inferred from the splitting of the valley or spin degeneracy of the Landau levels, which can be observed in spectroscopy measurements (Jiang et al., 2007a; Sadowski et al., 2006), or in the appearance of new quantum Hall plateaus (Abanin et al., 2007b; Giesbers et al., 2007; Goswami et al., 2007; Jiang et al., 2007b; Zhang et al., 2006).

Interactions can lead to new phases when their effect overcomes that of disorder. An analysis of the competition between disorder and interactions is found in ref. (Nomura and MacDonald, 2007). The energy splitting of the different broken symmetry phases, in a clean system, is determined by lattice effects, so that it is reduced by factors of order $a / l _ { B }$ , where $a$ is a length of the order of the lattice spacing, and $l _ { B }$ is the magnetic length (Alicea and Fisher, 2006, 2007; Goerbig et al., 2006; Wang et al., 2007). The combination of disorder and a magnetic field may also lift the degeneracy between the two valleys, favoring valley polarized phases (Abanin et al., 2007a).

In addition to phases with enhanced ferromagnetism or with broken valley symmetry, interactions at high magnetic fields can lead to excitonic instabilities (Gusynin et al., 2006) and Wigner crystal phases (Zhang and Joglekar, 2007). When only the $n = 0$ state is occupied, the Landau levels have all their weight in a given sublattice. Then, the breaking of valley degeneracy can be associated with a charge density wave, which opens a gap (Fuchs and Lederer, 2007). It is interesting to note that in these phases new collective excitations are possible (Doretto and Morais Smith, 2007).

Interactions modify the edge states in the quantum Hall regime. A novel phase can appear when the $n = 0$ is the last filled level. The Zeeman splitting shifts the electron and hole like chiral states, which disperse in opposite directions near the boundary of the sample. The resulting level crossing between an electron like level with spin anti-parallel to the field, and a hole like level with spin parallel to the field, may lead to Luttinger liquid features in the edge states (Abanin et al., 2007b; Fertig and Brey, 2006).

# VI. CONCLUSIONS

Graphene is a unique system in many ways. It is truly 2D, has unusual electronic excitations described in terms of Dirac fermions that move in a curved space, it is an interesting mix of a semiconductor (zero density of states) and a metal (gaplessness), and has properties of soft matter. The electrons in graphene seem to be almost insensitive to disorder and electron-electron interactions and have very long mean free paths. Hence, graphene's properties are rather different from what is found in usual metals and semiconductors. Graphene has also a very robust but flexible structure with unusual phonon modes that do not exist in ordinary 3D solids. In some sense, graphene brings together issues in quantum gravity and particle physics, and also from soft and hard condensed matter. Interestingly enough, these properties can be easily modified with the application of electric and magnetic fields, addition of layers, by control of its geometry, and chemical doping. Moreover, graphene can be directly and relatively easily probed by various scanning probe techniques from mesoscopic down to atomic scales, because it is not buried inside a 3D structure. This makes graphene one of the most versatile systems in condensed matter research.

Besides the unusual basic properties, graphene has the potential for a large number of applications (Geim and Novoselov, 2007), from chemical sensors (Chen et al., 2007c; Schedin et al., 2007) to transistors (Nilsson et al., 2007b; Oostinga et al., 2007). Graphene can be chemically and/or structurally modified in order to change its functionality and henceforth its potential applications. Moreover, graphene can be easily obtained from graphite, a material that is abundant on earth's surface. This particular characteristic makes graphene one of the most readily available materials for basic research since it frees economically challenged research institutions in developing countries from the dependence of expensive sample growing techniques.

Many of graphene's properties are currently subject of intense research and debate. The understanding of the nature of the disorder and how it affects the transport properties (a problem of fundamental importance for applications), the effect of phonons on electronic transport, the nature of electron-electron interactions and how they modify its physical properties are research areas that are still in their infancy. In this review, we have touched only the surface of a very deep sea that still has to be explored.

Whereas hundreds of papers have been written on monolayer graphene in the last few years, only a small fraction actually deals with multilayers. The majority of the theoretical and experimental efforts have been concentrated on the single layer, perhaps because of its simplicity and the natural attraction that a one atom thick material, which can be produced by simple methods in almost any laboratory in the world, creates for human imagination. Nevertheless, few layer graphene is equally interesting and unusual with a technological potential perhaps bigger than the single layer. Indeed, the theoretical understanding and experimental exploration of multilayers is far behind the single layer. This is a fertile and open field of research for the future.

Finally, we have focused entirely on pure carbon graphene where the band structure is dominated by the Dirac description. Nevertheless, chemical modification of graphene can lead to entirely new physics. Depending on the nature of chemical dopants and how they are introduced into the graphene lattice (adsorption, substitution, or intercalation) the results can be many. Small concentrations of adsorbed alkali metal can be used to change the chemical potential while adsorbed transition elements can lead to strong hybridization effects that affect the electronic structure. In fact, the introducion of dand f-electron atoms in the graphene lattice may produce a significant enhancement of the electron-electron interactions. Hence, it is easy to envision a plethora of manybody effects that can be induced by doping and have to be studied in the context of Dirac electrons: Kondo effect, ferromagnetism, antiferromagnetism, charge and spin density waves. The study of chemically induced many-body effects in graphene would add a new chapter to the short but fascinating history of this material. Only future will tell but the potential for more amazement is lurking on the horizon.

# VII. ACKNOWLEDGMENTS

We have benefited immensely from discussions with many colleagues and friends in the last few years but we would like to thank especially, Boris Altshuler, Eva Andrei, Alexander Balatsky, Carlo Beenakker, Sankar Das Sarma, Walt de Heer, Millie Dresselhaus, Vladimir Falko, Andrea Ferrari, Herb Fertig, Eduardo Fradkin, Ernie Hill, Mihail Katsnelson, Eun-Ah Kim, Philip Kim, Valery Kotov, Alessandra Lanzara, Leonid Levitov, Allan MacDonald, Serguey Morozov, Johan Nilsson, Vitor Pereira, Philip Phillips, Ramamurti Shankar, João Lopes dos Santos, Shan-Wen Tsai, Bruno Uchoa, and Maria Vozmediano.

N.M.R.P. acknowledges financial support from POCI 2010 via project PTDC/FIS/64404/2006. F.G. was supported by MEC (Spain) grant No. FIS2005-05478-C02-01 and EU contract 12881 (NEST). A. H. C. N was supported through NSF grant DMR-0343790. K.S.N. and A. K. G. were supported by EPSRC (UK) and the Royal Society.

# References

Abanin, D. A., P. A. Lee, and L. S. Levitov, 2006, Phys. Rev. Lett. 96, 176803.   
Abanin, D. A., P. A. Lee, and L. S. Levitov, 2007a, Solid State Comm. 143, 77.   
Abanin, D. A., and L. S. Levitov, 2007, Science 317, 641.   
Abanin, D. A., K. S. Novoselov, U. Zeitler, P. A. Lee, A. K. Geim, and L. S. Levitov, 2007b, Phys. Rev. Lett. 98, 196806.   
Abergel, D. S. L., A. Russell, and V. I. Fal'ko, 2007, Appl. Phys. Lett. 91, 063125.   
Abrikosov, A. A., 1998, Phys. Rev. B 58, 2788.   
Adam, S., E. H. Hwang, V. M. Galitski, and S. das Sarma, 2007, Proc. Natl. Acad. Sci. USA 104, 18392.   
Adebpour, N., M. Neek-Amal, R. A. a nd F. Shahbazi, N. Nafari, and M. Reza Rahimi Tabar, 2007, Phys. Rev. B 76, 195407.   
Affoune, A. M., B. L. V. Prasad, H. Saito, T. Enoki, Y. Kaburagi, and Y. Hishiyama, 2001, Chem. Phys. Lett. 348, 17.   
Akhmerov, A. R., and C. W. J. Beenakker, 2007, eprint arXiv:0710.2723.   
Aleiner, I. L., and K. B. Efetov, 2006, Phys. Rev. Lett. 97, 236801.   
Alicea, J., and M. P. A. Fisher, 2006, Phys. Rev. B 74, 075422.   
Alicea, J., and M. P. A. Fisher, 2007, Solid State Comm. 143, 504.   
Altland, A., 2006, Phys. Rev. Lett. 97, 236802.   
Ando, T., 1974a, J. Phys. Soc. Jpn. 36, 1521.   
Ando, T., 1974b, J. Phys. Soc. Jpn. 37, 622.   
Ando, T., 1974c, J. Phys. Soc. Jpn. 37, 1233.   
Ando, T., 1975, J. Phys. Soc. Jpn. 38, 989.   
Ando, T., 2000, J. Phys. Soc. Jpn. 69, 1757.   
Ando, T., 2006a, J. Phys. Soc. Jpn. 75, 124701.   
Ando, T., 2006b, J. Phys. Soc. Jpn. 75, 074716.   
Ando, T., 2007a, J. Phys. Soc. Jpn. 76, 024712.   
Ando, T., 2007b, J. Phys. Soc. Jpn. 76, 104711.   
Ando, T., T. Nakanishi, and R. Saito, 1998, J. Phys. Soc. Jpn. 67, 2857.   
Ando, T., and Y. Uemura, 1974, J. Phys. Soc. Jpn. 36, 959.   
Ando, T., Y. Zheng, and H. Suzuura, 2002, J. Phys. Soc. Jpn. 71, 1318.   
Andreoni, W., 200, The Physics of Fullerene-Based and Fullerene-related materials (Springer).   
Areshkin, D. A., and C. T. White, 2007, Nano Lett. 7, 204.   
Ashcroft, N. W., and N. D. Mermin, 1976, Solid State Physics (Saunders College, Philadelphia, PA).   
Bacon, G. E., 1950, Acta Crystalographica 4, 320.   
Baeriswyl, D., D. K. Campbell, and S. Mazumdar, 1986, Phys. Rev. Lett. 56, 1509.   
Bak, P., 1982, Rep. Prog. Phys. 45, 587.   
Balatsky, A. V., I. Vekhter, and J.-X. Zhu, 2006, Rev. Mod. Phys. 78, 373.   
Bar V. W., Y. Zhang, Y. Yayon, A. Bostwick, T. Ohta, J. L. McChesney, K. Horn, E. Rotenberg, and M. F. Crommie, 2007, Appl. Phys. Lett. 91, 122102.   
Bardarson, J. H., J. Tworzydlo, P. W. Brouwer, and C. W. J. Beenakker, 2007, Phys. Rev. Lett. 99, 106801.   
Barlas, Y., T. Pereg-Barnea, M. Polini, R. Asgari, and A. H. MacDonald, 2007, Phys. Rev. Lett. 98, 236601.   
Barone, V., O. Hod, and G. E. Scuseria, 2006, Nano Letters 6, 2748.   
Baym, G., 1969, Lectures on quantum mechanics (Benjamin, New York).   
Baym, G., and S. A. Chin, 1976, Nuclear Physics A262, 527.   
Bena, C., and S. A. Kivelson, 2005, Phys. Rev. B 72, 125432.   
Berger, C., Z. Song, X. Li, X. Wu, N. Brown, C. Naud, D. Mayou, T. Li, J. Hass, A. N. Marchenkov, E. H. Conrad, P. N. First, et al., 2006, Science 312, 1191.   
Berger, C., Z. M. Song, T. B. Li, X. B. Li, A. Y. Ogbazghi, R. Feng, Z. T. Dai, A. N. Marchenkov, E. H. Conrad, P. N. First, and W. A. de Heer, 2004, J. Phys. Chem. B 108, 19912.   
Bergman, G., 1984, Phys. Rep. 107,1.   
Bevig, B.A. T. L. Hugs, S. Rag, an D. P.Aroas, 2007, eprint cond-mat/0701436.   
Berry, M. V., and R. J. Modragon, 1987, Proc. R. Soc. Lond. A 412, 53.   
Bhattacharjee, S., and K. Sengupta, 2006, Phys. Rev. Lett. 97, 217001.   
Birrell, N. D., and P. C. W. Davies, 1982, Quantum Fields in Curved Space (Cambridge Univ. Press, Cambridge).   
Biswas, R. B., S. Sachdev, and D. T. Son, 2007, Phys. Rev. B 76, 205122.   
Blake, P., K. S. Novoselov, A. H. Castro Neto, D. Jiang, R. Yang, T. J. Booth, A. K. Geim, and E. W. Hill, 2007, Appl. Phys. Lett. 91, 063124.   
Bliokh, K. Y., 2005, Phys. Lett. A 344, 127.   
Bol, A. J. V.J. E. Crombeen, an A. V.Tooren, 1975, Surf. Sci. 48, 463.   
Bonni N. M. Lazzeri, N. Marzari, andF. Mauri, 2007, Phys. Rev. Lett. 99, 176802.   
Bostwick, A., T. Ohta, J. L. McChesney, K. V. Emtsev, T. Seyller, K. Horn, and E. Rotenberg, 2007a, New J. Phys. 9, 385.   
Bosick, A., T. Ohta, T. Seyr, K. Horn, andE. Rot, 2007b, Nature Physics 3, 36.   
Boyle, W. S., and P. Nozières, 1958, Phys. Rev. 111, 782.   
Brandt, N. B., S. M. Chudinov, and Y. G. Ponomarev, 1988, in Modern Problems in Condensed Matter Sciences, edited by V. M. Agranovich and A. A. Maradudin (North Holland (Amsterdam)), volume 20.1.   
Brey, L., and H. Fertig, 2006a, Phys. Rev. B 73, 195408.   
Brey, L., and H. Fertig, 2006b, Phys. Rev. B 73, 235411.   
Brey, L., H. A. Fertig, and S. D. Sarma, 2007, Phys. Rev. Lett. 99, 116802.   
Bunch, J. S., A. M. van der Zande, S. S. Verbridge, I. W. Frank, D. M. Tanenbaum, J. M. Parpia, H. G. Craighead, and P. L. McEuen, 2007, Science 315, 490.   
Calandra, M., and F. Mauri, 2007, Phys. Rev. B 76, 199901.   
I. . Bo F. M C. . La an .A.B, 2007, Appl. Phys. Lett. 91, 201904.   
Calogeracos, A., and N. Dombey, 1999, Contemp. Phys. 40, 313.   
Campagnoli, G., and E. Tosatti, 1989, in Progress on Electron Properties of Metals, edited by R. Girlanda et al (Kluwer Academic Publishing), p. 337.   
Cancado, L. G., M. A. Pimenta, R. B. R. Neves, G. MedeirosRibeiro, T. Enoki, Y. Kobayashi, K. Takai, K.-I. Fukui, M. S. Dresselhaus, R. Saito, and A. Jorio, 2004, Phys. Rev. Lett. 93, 047403.   
Casiraghi, C., A. Hartschuh, E. Lidorikis, H. Qian, H. Harutyunyan, T. Gokus, K. S. Novoselov, and A. C. Ferrari, 2007, Nano Lett. 7, 2711.   
Cassanello, C. R., and E. Fradkin, 1996, Phys. Rev. B 53, 15079.   
Cassanello, C. R., and E. Fradkin, 1997, Phys. Rev. B 56, 11246.   
Castillo, H.E. C.deC.Chaon, E.Fradkin, P. M.Golb, and C. Mudry, 1997, Phys. Rev. B 56, 10668.   
Castro, E. V., K. S. Novoselov, S. V. Morozov, N. M. R. Peres, J. Lopes dos Santos, J. Nilsson, F. Guinea, A. K. Geim, and A. H. Castro Neto, 2007a, Phys. Rev. Lett. 99, 216802.   
Castro, E. V., N. M. R. Peres, J. M. B. Lopes dos Santos, A. H. Castro Neto, and F. Guinea, 2007b, Phys. Rev. Lett. 100, 026802.   
Castro Neto, A. H., 2007, Nature Materials 6, 176.   
Castro Neto, A. H., and F. Guinea, 2007, Phys. Rev. B 75, 045404.   
Castro Neto, A. H., F. Guinea, and N. M. R. Peres, 2006a, Physics World 19, 33.   
Castro Neto, A. H., F. Guinea, and N. M. R. Peres, 2006b, Phys. Rev. B 73, 205408.   
Chaikn, P n T.C.Lubensky, 1995 Inri ondensed Matter Physics (Cambridge University Press).   
Chakravarty, S., and A. Schmid, 1986, Phys. Rep. 140, 193.   
Chamon, C. C., C. Mudry, and X.-G. Wen, 1996, Phys. Rev. B 53, R7638.   
Charlier, J.-C., X. Blase, and S. Roche, 2007, Rev. Mod. Phys. 79, 677.   
Charlier, J. C., J. P. Michenaud, X. Gonze, and J. P. Vigneron, 1991, Phys. Rev. B 44, 13237.   
Cheianov, V. V., V. Fal'ko, and B. L. Altshuler, 2007a, Science 315, 252.   
Cheiaov, V. V., and V. I. Fal'ko, 2006, Phys. Rev. B74, 041403.   
Cheianov, V. V., V. I. Fal'ko, B. L. Altshuler, and I. L. Aleiner, 2007b, Phys. Rev. Lett. 99, 176801.   
Chen, H.-Y., V. Apalkov, and T. Chakraborty, 2007a, Phys. Rev. Lett. 98, 186803.   
Chen, J. H., C. Jang, M. S. Fuhrer, E. D. Williams, and M. Ishigami, 2007b, eprint arXiv:0708.2408.   
Che, Z.Y.M. Lin, M. J. Rooks, and P.Avois, 2007, Physica E 40/2, 228.   
Cho, S., Y.-F. Chen, and M. S. Fuhrer, 2007, Appl. Phys. Lett. 91, 123105.   
Coey, J. M. D., M. Venkatesan, C. B. Fitzgerald, A. P. Douvalis, and I. S. Sanders, 2002, Nature 420, 156.   
Cortijo, A., and M. A. H. Vozmediano, 2007a, Nucl. Phys. B 763, 293.   
Cortijo, A., and M. A. H. Vozmediano, 2007b, Europhys. Lett. 77, 47002.   
Couraux, J., A. T. N'Diaye, C. Busse, and T. Michely, 2008, Nano Letters DOI:10.1021/nl0728874.   
Cserti, J., A. Csordás, and G. Dávid, 2007a, Phys. Rev. Lett. 99, 066802.   
Cserti, J., A. Palyi, and C. Peterfalvi, 2007b, Phys. Rev. Lett. 99, 246801.   
Dahal, H., Y. N. Joglekar, K. Bedell, and A. V. Balatsky, 2006, Phys. Rev. B 74, 233405.   
Dahal, H. P., T. O. Wehling, K. S. Bedell, J.-X. Zhu, and A. V. Balatsky, 2007, eprint cond-mat/0706.1689.   
Das, A., B. Chakraborty, and A. K. Sood, 2007, eprint arXiv:0710.4160.   
Den, R. S.K.C.ChugR. J. Nicolas, K.S.Novo, and A. K. Geim, 2007, Phys. Rev. B 76, 081406(R).   
Dell'Anna, L., 2006, Nucl.Phys.B 758, 255.   
Dharma-Wardana, M. W. C., 2007, J. Phys.:Condens. Matter 19, 386228.   
Dienwiebel, M., G. S. Verhoeven, N. Pradeep, J. W. M. Frenken, J. A. Heimberg, and H. W. Zandbergen, 2004, Phys. Rev. Lett. 92, 126101.   
Dillon, R. O., I. L. Spain, and J. W. McClure, 1977, J. Phys. Chem. Sol. 38, 635.   
DiVincenzo, D. P., and E. J. Mele, 1984, Phys.Rev.B 29, 1685.   
Dong, S.-H., X.-W. Hou, and Z.-Q. Ma, 1998, Phys. Rev. A 58, 2160.   
Doretto, R. L., and C. Morais Smith, 2007, eprint condmat/0704.3671.   
Dresselhaus, G., and M. S. Dresselhaus, 1965, Phys. Rev. 140, A401.   
Dresselhaus, M. S., and G. Dresselhaus, 2002, Advances in Physics 51, 1.   
Dresselhaus, M. S., G. Dresselhaus, J. E. Fischer, and M. J. Moran, 1983, Intercalated graphite (North-Holland, New York).   
Dresselhaus, M. S., and J. G. Mavroides, 1964, IBM J. Res. Dev 8, 262.   
Dugaev, V. K., V. I. Litvinov, and J. Barnas, 2006, Phys. Rev. B 74, 224438.   
Duke, C. B., 1968, Phys. Rev. 168, 816.   
Eisenstein, J. P., L. N. Pfeiffer, and K. W. West, 1994, Phys. Rev. B 50, 1760.   
Eizenberg, M., and J. M. Blakely, 1979, Surf. Sci. 82, 228.   
Esquinazi, P., A. Setzer, R. Höhne, C. Semmelhack, Y. Kopelevich, D. Spemann, T. Butz, B. Kohlstrunk, and M. Lösch, 2002, Phys. Rev. B 66, 024429.   
Esquinazi, P. D. Spemann, R. Höhne, A. Setzer, K.H. Han, and T. Butz, 2003, Phys. Rev. Lett. 91, 227201.   
Eun-Ah Kim, and A. H. Castro Neto, 2007, eprint condmat/0702562.   
Fasolino, A., J. H. Los, and M. I. Katsnelson, 2007, Nat. Mat. 6, 858.   
Faugeras, C., A. Nerriere, M. Potemski, A. Mahmood, E. Dujardin, C. Berger, and W. A. de Heer, 2007, eprint arXiv:0709.2538.   
Fr,BJ.Tolksor an Zeier, 2007 Quity (Birkhäuser).   
Ferrari, A. C., J. C. Meyer, V. Scardaci, C. Casiraghi, M. Lazzeri, F. Mauri, S. Piscanec, D. Jiang, K. S. Novoselov, and A. K. G. S. Roth, 2006, Phys. Rev. Lett. 97, 187401.   
Fertig, H. A., and L. Brey, 2006, Phys. Rev. Lett. 97, 116805.   
Fogler, M. M., L. I. Glazman, D. S. Novikov, and B. I. Shklovskii, 2007a, eprint arXiv:0710.2150.   
Fogler, M. M., D. S. Novikov, and B. I. Shklovskii, 2007b, Phys. Rev. B 76, 233402.   
Forbeaux, I., J.-M. Themlin, and J.-M. Debever, 1998, Phys. Rev. B 58, 16396.   
Foster, M. S., and A. W. W. Ludwig, 2006a, Phys. Rev. B 73, 155104.   
Foster, M. S., and A. W. W. Ludwig, 2006b, Phys. Rev. B 74, 241102.   
Fradkin, E., 1986a, Phys. Rev. B 33, 3257.   
Fradkin, E., 1986b, Phys. Rev. B 33, 3263.   
Fritz, L., S. Florens, and M. Vojta, 2006, Phys. Rev. B 74, 144410.   
Fuchs, J.-N., and P. Lederer, 2007, Phys. Rev. Lett. 98, 016803.   
Fujita, M., K. Wakabayashi, K. Nakada, and K. Kusakabe, 1996, J. Phys. Soc. Jpn. 65, 1920.   
Fukuyama, H., 1971, Prog. Theor. Phys. 45, 704.   
Gasparoux, H., 1967, Carbon 5, 441.   
Geim, A. K., and A. H. MacDonald, 2007, Physics Today 60, 35.   
Geim, A. K., and K. S. Novoselov, 2007, Nature Materials 6, 183.   
de Gennes, P. G., 1964, Rev. Mod. Phys. 36, 225.   
Ghosal, A., P. Goswami, and S. Chakravarty, 2007, Phys. Rev. 75, 115123.   
Giesbers, A. J. M., U. Zeitler, M. I. Katsnelson, L. A. Ponomarenko, T. M. G. Mohiuddin, and J. C. Maan, 2007, Phys. Rev. Lett. 99, 206803.   
Giovannetti, G., P. A. Khomyakov, G. Brocks, P. J. Kelly, and J. van der Brink, 2007, Phys. Rev. B 76, 073103.   
Giuliani, G. F., and G. Vignale, 2005, Quantum theory of the electron liquid (Cambridge Press).   
Goerbig, M. O., J.-N. Fuchs, K. Kechedzhi, and V. I. Fal'ko, 2007, Phys. Rev. Lett. 99, 087402.   
Goerbig, M. O., R. Moessner, and B. Douot, 2006, Phys. Rev. B 74, 161407.   
González, J., F. Guinea, and M. A. H. Vozmediano, 1992, Phys. Rev. Lett. 69, 172.   
González, J., F. Guinea, and M. A. H. Vozmediano, 1993a, Mod. Phys. Lett. B7, 1593.   
González, J., F. Guinea, and M. A. H. Vozmediano, 1993b, Nucl. Phys. B 406 [FS], 771.   
González, J., F. Guinea, and M. A. H. Vozmediano, 1994, Nucl. Phys. B 424, 596.   
González, J., F. Guinea, and M. A. H. Vozmediano, 1996, Phys. Rev. Lett. 77, 3589.   
González, J., F. Guinea, and M. A. H. Vozmediano, 1999, Phys. Rev. B 59, R2474.   
González, J., F. Guinea, and M. A. H. Vozmediano, 2001, Phys. Rev. B 63, 134421.   
Gorbachev, R. V., F. V. Tikhonenko, A. S. Mayorov, D. W. Horsell, and A. K. Savchenko, 2007, Phys. Rev. Lett. 98, 176805.   
Gorbar, E. V., V. P. Gusynin, V. A. Miransky, and I. A. Shovkovv. 2002. Phvs. R.ev. B 66. 045108.   
Goswami, P., X. Jia, and S. Chakravarty, 2007, Phys. Rev. B 76, 205408.   
Graf, D., F. Molitor, K. Ensslin, C. Stampfer, A. Jungen, C. Hierold, and L. Wirtz, 2007, Nano Letters 7, 238.   
Guinea, F., 2007, Phys. Rev. B 75, 235433.   
Guinea, F., A. H. Castro Neto, and N. M. R. Peres, 2006, Phys. Rev. B 73, 245426.   
Guinea, F., M. I. Katsnelson, and M. A. H. Vozmediano, 2007, Phys. Rev. B 76, 235309.   
Gumbs, G., and P. Fekete, 1997, Phy. Rev. B 56, 3787.   
Gunlycke, D., D. A. Areshkin, and C. T. White, 2007, Appl. Phys. Lett. 90, 12104.   
Gupta, A., G. Chen, P. Joshi, S. Tadigadapa, and P. C. Eklund, 2006, Nano Letters 12, 2667.   
Gusynin, V. P., V. A. Miransky, S. G. Sharapov, and I. A. Shovkovy, 2006, Phys. Rev. B 74, 195429.   
Gusynin, V. P., and S. G. Sharapov, 2005, Phys. Rev. Lett. 95, 146801.   
Gusynin, V. P., S. G. Sharapov, and J. P. Carbotte, 2007, J. Phys.: Condens.Matter 19, 026222.   
Haldane, F. D. M., 1988, Phys. Rev. Lett. 61, 2015.   
Han, M. Y., B. Özyilmaz, Y. Zhang, and P. Kim, 2007, Phys. Rev. Lett. 98, 206805.   
Harper, P. G., 1955, Proc. Phys. Soc. London A 68, 874.   
Harrison, W. A., 1980, Solid State Theory (Dover).   
Hass, J., R. Feng, J. E. Millán-Otoya, X. Li, M. Sprinkle, P. N. First, W. A. de Heer, and E. H. Conrad, 2007a, Phys. Rev. B 75, 214109.   
Hass, J., F. Varchon, J. E. Millan-Otoya, M. Sprinkle, W. A. de Heer, C. Berger, P. N. First, L. Magaud, and E. H. Conrad, 2007b, eprint cond-mat/0706.2134.   
de Heer, W. A., C. Berger, X. Wu, P. N. First, E. H. Conrad, X. Li, T. Li, M. Sprinkle, J. Hass, M. L. Sadowski, M. Potemski, and G. Martinez, 2007, Sol. State Comm. 143, 92.   
Heersche, H. B., P. Jarillo-Herrero, J. B. Oostinga, L. M. K. Vandersypen, and A. Morpurgo, 2007, Nature 446, 56.   
Hentschel, M., and F. Guinea, 2007, eprint condmat/0705.0522.   
Herbut, I. F., 2006, Phys. Rev. Lett. 97, 146401.   
Herbut, I. F., 2007, Phys. Rev. B 75, 165411.   
Heremans, J., C. H. Olk, and D. T. Morelli, 1994, Phys. Rev. B 49, 15122.   
Hill, E. W., A. K. Geim, K. Novoselov, F. Schedin, and P. Blake, 2007, IEEE Trans. Magn. 42, 2694.   
Himsl, F. J.K.Chrian, P. Hen, D.E. Easn, and P. J. Feibelman, 1982, Surf. Sci. 115, L159.   
Hobson, J. P., and W. A. Nierenberg, 1953, Phys. Rev. 89, 662.   
Hod, O., V. Barone, J. E. Peralta, and G. E. Scuseria, 2007, Nano Letters 7, 2295.   
Horovitz, B., and P. L. Doussal, 2002, Phys. Rev. B 65, 125323.   
Hou, C.-Y., C. Chamon, and C. Mudry, 2007, Phys. Rev. Lett. 98, 186809.   
Huar, B.J. A. Sulpizio, .Staner, K.To, B.an a D. Goldhaber-Gordon, 2007, Phys. Rev. Lett. 98, 236803.   
Huertas-Hernando, D., F. Guinea, and A. Brataas, 2006, Phys. Rev. B 74, 155426.   
Hwang, E. H., and S. Das Sarma, 2007, Phys. Rev. B 75, 205418.   
Hwang, E. H., B. Y.-K. Hu, and S. D. Sarma, 2007, Phys. Rev. B 76, 115434.   
Ishigami, M., J. H. Chen, D. W. G. Cullan, M. S. Fuhrer, and

E. D. Williams, 2007, Nano Letters 7, 1643. Itzykson, C., and J.-B. Zuber, 2006, Quantum Field Theory (Dover). Jackiw, R., and C. Rebbi, 1976, Phys. Rev. D 13, 3398. Jiang, Z., E. A. Henriksen, L. C. Tung, Y.-J. Wang, M. E. Schwartz, M. Y. Han, P. Kim, and H. L. Stormer, 2007a, Phys. Rev. Lett. 98, 197403. Jiang, Z., Y. Zhang, H. L. Stormer, and P. Kim, 2007b, Phys. Rev. Lett. 99, 106802. de Juan, F., A. Cortijo, and M. A. H. Vozmediano, 2007, Phys. Rev. B 76, 165409. Kane, C. L., and E. J. Mele, 1997, Phys. Rev. Lett. 78, 1932. Kane, C. L., and E. J. Mele, 2004, Phys. Rev. Lett. 93, 197402. Kane, C. L., and E. J. Mele, 2005, Phys. Rev. Lett. 95, 226801. Katsnelson, M. I., 2006a, Phys. Rev. B 74, 201401. Katsnelson, M. I., 2006b, Eur. J. Phys. B 51, 157. Katsnelson, M. I., 2007a, Eur. Phys. J B 57, 225. Katsnelson, M. I., 2007b, Materials Today, 10, 20. Katsnelson, M. I., 2007c, Phys. Rev. B 76, 073411. Katsnelson, M. I., and A. K. Geim, 2008, Phil. Trans. R. Soc. A 366, 195. Katsnelson, M. I., and K. S. Novoselov, 2007, Sol. Stat. Comm. 143, 3. Katsnelson, M. I., K. S. Novoselov, and A. K. Geim, 2006, Nature Physics 2, 620. Kechedzhi, K., V. I. Fal'ko, E. McCann, and B. L. Altshuler, 2007, Phys. Rev. Lett. 98, 176806. Khveshchenko, D. V., 2001, Phys. Rev. Lett. 87, 246802. Khveshchenko, D.V., 2007, eprint cond-mat/0705.4105. Khveshchenko, D. V., and W. F. Shively, 2006, Phys. Rev. B 73, 115104. Kobayashi, Y., K. Fukui, T. Enoki, K. Kusakabe, and Y. Kaburagi, 2005, Phys. Rev. B 71, 193406. Kolesnikov, D. V., and V. A. Osipov, 2004, JETP Letters 79, 532. Kolesnikov, D. V., and V. A. Osipov, 2006, Eur. J. Phys. B 49, 465. Kolezhuk, A., S. Sachdev, R. B. Biswas, and P. Chen, 2006, Phys. Rev. B 74, 165114. Kopelevich, Y., and P. Esquinazi, 2006, eprint condmat/0609497. Kopelevich, Y., J. C. Medina Pantoja, R. R. da Silva, F. Mrowka, and P. Esquinazi, 2006, Physics Letters A 355, 233. Kopelevich, Y., J. H. S. Torres, R. R. da Silva, F. Mrowka, H. Kempa, and P. Esquinazi, 2003, Phys. Rev. Lett. 90, 156402. Koshino, M., and T. Ando, 2006, Phys. Rev. B 73, 245403. Koshino, M., and T. Ando, 2007, Phys. Rev. B 75, 235333. Kumazaki, H., and D. S. Hirashima, 2006, J. Phys. Soc. Jpn. 75, 053707. Kumazaki, H., and D. S. Hirashima, 2007, J. Phys. Soc. Jpn. 76, 064713. Kusminskiy, S. V., J. Nilsson, D. K. Campbell, and A. H. Castro Neto, 2007, eprint cond-mat/0706.2359. Lammert, P. E., and V. H. Crespi, 2004, Phys. Rev. B 69, 035406. Lanau, L. D. an E. M.Lihiz, 1959, Theory f Elasicy (Pergamon Press (London)): Landau, L. D., and E. M. Lifshitz, 1981, Quantum Mechanics: Non-Relativistic Theory (Pergamon Press (London)). Laughlin, R. B.. 1981, Phvs. Rev. B 23. 5632.

Lazzeri, M., and F. Mauri, 2006, Phys. Rev. Lett. 97, 266407.   
LeClair, A., 2000, Phys. Rev. Lett. 84, 1292.   
Lee, P. A., 1993, Phys. Rev. Lett. 71, 1887.   
Lee, P. A., and T. V. Ramakrishnan, 1985, Rev. Mod. Phys. 57, 287.   
Lee, S.-S., and P. A. Lee, 2005, Phys. Rev. Lett. 95, 036403.   
Leenaerts, O., B. Partoens, and F. M. Peeters, 2007, eprint arXiv:0710.1757.   
Lemme, M. C., T. J. Echtermeyer, M. Baus, and H. Kurz, 2007, IEEE Electron Device Letters 28, 282.   
Lenosky, T., X. Gonze, M. Teter, and V. Elser, 1992, Nature 355, 333.   
Lewenkopf, C. H., E. R. Mucciolo, and A. H. Castro Neto, 2007, eprint arXiv:0711.3202.   
Li, G., and E. Y. Andrei, 2007, Nature Physics 3, 623.   
Li, Z. Q., S.-W. Tsai, W. J. Padilla, S. V. Dordevic, K. S. Burch, Y. J. Wang, and D. N. Basov, 2006, Phys. Rev. B 74, 195404.   
Lopes dos Santos, J. M. B., N. M. R. Peres, and A. H. Castro Neto, 2007, Phys. Rev. Lett. 99, 256802.   
Louis, E., J. A. Vergés, F. Guinea, and G. Chiappe, 2007, Phys. Rev. B 75, 085440.   
Ludwig, A. W. W., M. P. A. Fisher, R. Shankar, and G. Grinstein, 1994, Phys. Rev. B 50, 7526.   
Lukose, V., R. Shankar, and G. Baskaran, 2007, Phys. Rev. Lett. 98, 16802.   
Luk'yanchuk, I. A., and Y. Kopelevich, 2004, Phys. Rev. Lett. 93, 166402.   
Mañes, J. L., 2007, Phys. Rev. B 76, 045430.   
Mañes, J. L., F. Guinea, and M. A. H. Vozmediano, 2007, Phys. Rev. B 75, 155424.   
Maiti, M., and K. Sengupta, 2007, Phys. Rev. B 76, 054513.   
Malard, L. M., J. Nilsson, D. C. Elias, J. C. Brant, F. Plentz, E.S.Alves, A. H. Castro Neto, and M. A. Pimenta, 2007, Phys. Rev. B 76, 201401.   
Malle, P. F.Varcon, C. Nau L. Mag, C. Berger, and J.-Y. Veuillen, 2007, Phys. Rev. B 76, 041403(R).   
Maple, M. B., 1998, Jou. Mag. Mag. Mat. 177, 18.   
Mariani, E., L. Glazman, A. Kamenev, and F. von Oppen, 2007, Phys. Rev. B 76, 165402.   
Mariani, E., and F. von Oppen, 2007, eprint arXiv:condmat/0707.4350.   
Martelo, L. M., M. Dzierzawa, L. Siffert, and D. Baeriswyl, 1997, Z. Physik B 103, 335.   
Martin, I., and Y. M. Blanter, 2007, eprint condmat/0705.0532.   
Martin, J., N. Akerman, G. Ulbricht, T. Lohmann, J. H. Smet, K. von Klitzing, and A. Yacoby, 2007, eprint condmat/0705.2180.   
MarT.B.R. H.Mi A. J.R. Siva an A. zi, 2007, Phys. Rev. Lett. 98, 196803.   
Matsui, T., H. Kambara, Y. Niimi, K. Tagami, M. Tsukada, and H. Fukuyama, 2005, Phys. Rev. Lett. 94, 227201.   
Mattausch, A., and O. Pankratov, 2007, eprint condmat/0704.0216.   
McCann, E., 2006, Phys. Rev. B 74, 161403.   
McCann, E., and V. I. Fal'ko, 2006, Phys. Rev. Lett. 96, 086805.   
McCann, E., K. Kechedzhi, V. I. Fal'ko, H. Suzuura, T. Ando, and B. L. Altshuler, 2006, Phys. Rev. Lett. 97, 146805.   
McClure, J. W., 1956, Phys. Rev. 104, 666.   
McClure, J. W., 1957, Phys. Rev. 108, 612.   
McClure, J. W., 1958, Phys. Rev. 112, 715.   
McClure, J. W., 1960, Phys. Rev. 119, 606.   
McClure, J. W., 1964, IBM J. Res. Dev. 8, 255.   
McClure, J 91, Physc  Smimeals n arr semiconductors (D. L. Carter and R. T. Bate, Pergamon Press, New York).   
McEuen, P. L., M. Bockrath, D. H. Cobden, Y.-G. Yoon, and S. G. Louie, 1999, Phys. Rev. Lett. 83, 5098.   
Meyer, J. C., A. K. Geim, M. I. Katsnelson, K. S. Novoselov, T. J. Booth, and S. Roth, 2007a, Nature 446, 60.   
Meyer, J. C., A. K. Geim, M. I. Katsnelson, K. S. Novoselov, D. Obergfell, S. Roth, C. Girit, and A. Zettl, 2007b, Solid State Commun. 143, 101.   
Miao, F., S. Wijeratne, U. Coskun, Y. Zhang, and C. N. Lau, 2007, Science 317, 1530.   
Milton Pereira Junior, J., P. Vasilopoulos, and F. M. Peeters, 2007, Nano Letters 7, 946.   
Min, H., J. Hill, N. Sinitsyn, B. Sahu, L. Kleinman, and A. MacDonald, 2006, Phys. Rev. B 74, 165310.   
Mishchenko, E. G., 2007, Phys. Rev. Lett. 98, 216801.   
Morozov, S., K. Novoselov, M. Katsnelson, F. Schedin, D. Jiang, and A. K. Geim, 2006, Phys. Rev. Lett. 97, 016801.   
Morozov, S. V., K. S. Novoselov, F. Schedin, D. Jiang, A. A. Firsov, and A. K. Geim, 2005, Phys. Rev. B 72, 201401.   
Morpurgo, A. F., and F. Guinea, 2006, Phys. Rev. Lett. 97, 196804.   
Mrozowski, S., 1952, Phys. Rev. 85, 609.   
Muñoz-Rojas, F., D. Jacob, J. Fernández-Rossier, and J. J. Palacios, 2006, Phys. Rev. B 74, 195417.   
Nair, R. R., U. Bangert, M. H. Gass, K. S. Novoselov, A. K. Geim, and A. L. Bleloch, 2008, eprint unpublished.   
Nakada, K., M. Fujita, G. Dresselhaus, and M. S. Dresselhaus, 1996, Phys. Rev. B 54, 17954   
Nakamura, M., 2007, Phys. Rev. B 76, 113301.   
Nakao, K., 1976, J. Phys. Soc. Jpn. 40, 761.   
Nelson, D., D. R. Piran, and S. Weinberg, 2004, Statistical mechanics of membranes and surfaces (World Scientific, Singapore).   
Nelson, D. R., and L. Peliti, 1987, J. Physique 48, 1085.   
Nersesyan, A. A., A. M. Tsvelik, and F. Wenger, 1994, Phys. Rev. Lett. 72, 2628.   
Niimi, Y., H. Kambara, T. Matsui, D. Yoshioka, and H. Fukuyama, 2006, Phys. Rev. Lett. 97, 236804.   
Nilsson, J., and A. H. Castro Neto, 2007, Phys. Rev. Lett. 98, 126801.   
Nilsson, J., A. H. Castro Neto, F. Guinea, and N. M. R. Peres, 2006a, Phys. Rev. Lett. 97, 266801.   
Nilsson, J., Å. H. Castro Neto, F. Guinea, and N. M. R. Peres, 2007a, eprint arXiv:0712.3259.   
Nilsson, J., A. H. Castro Neto, F. Guinea, and N. M. R. Peres, 2007b, Phys. Rev. B 76, 165416.   
Nilsson, J., A. H. Castro Neto, N. M. R. Peres, and F. Guinea, 2006b, Phys. Rev. B 73, 214418.   
Nilsson, J., A. H. Castro Neto, N. M. R. Peres, and F. Guinea, 2006c, Phys. Rev. B 73, 214418.   
Nomura, K., M. Koshino, and S. Ryu, 2007, Phys. Rev. Lett. 99, 146806.   
Nomura, K., and A. H. MacDonald, 2006, Phys. Rev. Lett. 96, 256602.   
Nomura, K., and A. H. MacDonald, 2007, Phys. Rev. Lett. 98, 076602.   
Novikov, D. S., 2007a, Phys. Rev. B 76, 245435.   
Novikov, D. S., 2007b, Appl. Phys. Lett. 91, 102102.   
Novikov, D. S., 2007c, Phys. Rev. Lett. 99, 056802.   
Novoselov, K. S., A. K. Geim, S. V. Morozov, D. Jiang, M. I. Katsnelson, I. V. Grigorieva, S. V. Dubonos, and A. A. Firsov, 2005a, Nature 438, 197.   
Novoselov, K. S., A. K. Geim, S. V. Morozov, D. Jiang, Y. Zhang, S. V. Dubonos, I. V. Gregorieva, and A. A. Firsov, 2004, Science 306, 666.   
Novoselov, K. S., D. Jiang, F. Schedin, T. J. Booth, V. V. Khotkevich, S. V. Morozov, and A. K. Geim, 2005b, Natl. Acad. Sci. USA 102, 10451.   
Novoselov, K. S., Z. Jiang, Y. Zhang, S. V. Morozov, H. L. U.Z J.   . r ., and A. K. Geim, 2007, Science 315, 1379.   
Novoselov, K. S., E. McCann, S. V. Mozorov, V. I. Fal'ko, M. I. Katsnelson, U. Zeitler, D. Jiang, F. Schedin, and A. K. Geim, 2006, Nature Physics 1, 177.   
Nozières, P., 1958, Phys. Rev. 109, 1510.   
Ohishi, M., M. Shiraishi, R. Nouchi, T. Nozaki, T. Shinjo, and Y. Suzuki, 2007, Jap. J. Appl. Phys. 46, L605.   
Ohta, K., 1968, Jpn. J. Appl. Phys. 10, 850.   
Ohta, K., 1971, J. Phys. Soc. Jpn. 31, 1627.   
Ohta, T., A. Bostwick, J. L. McChesney, T. Seyller, K. Horn, and E. Rotenberg, 2007, Phys Rev. Lett. 98, 206802.   
T. BckT. y .Hor anRo, 2006, Science 313, 951.   
Ostinga, J. B. H. B. Heerce, X. Li, A. Morugo, and L. M. K. Vandersypen, 2007, Nature Materials DOI: 10.1038/nmat2082.   
Oshima, C., and A. Nagashima, 1997, J. Phys.: Condens. Matter 9, 1.   
Osipov, V. A., E. A. Kochetov, and M. Pudlak, 2003, JETP 96, 140.   
Ossipov, A., M. Titov, and C. W. J. Beenakker, 2007, Phys. Rev. B 75, 241401.   
Ostrovsky, P. M., I. V. Gornyi, and A. D. Mirlin, 2006, Phys. Rev. B 74, 235443.   
Ostrovsky, P. M., I. V. Gornyi, and A. D. Mirlin, 2007, Phys. Rev. Lett. 98, 256801.   
Öilmaz, B. P.JarioHerreo,D.Efetov, D.A.Aban, L. S. Levitov, and P. Kim, 2007, Phys. Rev. Lett. 99, 166804.   
Paiva, T., R. T. Scalettar, W. Zheng, R. R. P. Singh, and J. Oitmaa, 2005, Phys. Rev. B 72, 085123.   
Parr, R. G., D. P. Craig, and I. G. Ross, 1950, J. Chem. Phys. 18, 1561.   
Pauling, L., 1972, The Nature of the Chemical Bond (Cornell U. P., Ithaca, NY).   
Peliti, L., and S. Leibler, 1985, Phys. Rev. B 54, 1690.   
Pereira, V. M., F. Guinea, J. M. B. L. dos Santos, N. M. R. Peres, and Å. H. Castro Neto, 2006, Phys. Rev. Lett. 96, 036801.   
Pereira, V. M., J. M. B. Lopes dos Santos, and A. H. Castro Neto, 2007a, eprint arXiv:0712.0806.   
Pereira, V. M., J. Nilsson, and A. H. Castro Neto, 2007b, Phys. Rev. Lett. 99, 166802.   
Peres, N. M. R., M. A. N. Araújo, and D. Bozi, 2004, Phys. Rev. B 70, 195122.   
Peres, N. M. R., and E. V. Castro, 2007, Jou. Phys. Cond. Mat. 19, 406231.   
Peres, N. M. R., A. H. Castro Neto, and F. Guinea, 2006a, Phys. Rev. B 73, 195411.   
Peres, N. M. R., A. H. Castro Neto, and F. Guinea, 2006b, Phys. Rev. B 73, 241403.   
Peres, N. M. R., F. Guinea, and A. H. Castro Neto, 2005, Phys. Rev. B 72, 174406.   
Peres, N. M. R., F. Guinea, and A. H. Castro Neto, 2006c, Phys. Rev. B 73, 125411.   
Peres, N. M. R., F. Guinea, and A. H. Castro Neto, 2006d, Annals of Physics 321, 1559.   
Peres, N. M. R., F. D. Klronomos, S.-W. Tsai, J. R. Santos, J. M. B. Lopes dos Santos, and A. H. Castro Neto, 2007a, Europhys. Lett. 80, 67007.   
Peres, N. M. R., J. M. Lopes dos Santos, and T. Stauber, 2007b, Phys. Rev. B 76, 073412.   
Petroski, H., 1989, The Pencil: A History of Design and Circumstan ce (Alfred Knopf, New York).   
Phillips, P., 2006, Annals of Physics 321, 1634.   
Pisana, S., M. Lazzeri, C. Casiraghi, K. S. Novoselov, A. K. Geim, A. C. Ferrari, and F. Mauri, 2007, Nature Materials 6, 198.   
Polini, M., R. Asgari, Y. Barlas, T. Pereg-Barnea, and A. H. MacDonald, 2007, Solid State Commun. 143, 58.   
Polkovnikov, A., 2002, Phys. Rev. B 65, 064503.   
Polkovnikov, A., S. Sachdev, and M. Vojta, 2001, Phys. Rev. Lett. 86, 296.   
Radzihovsky, L., and P. Le Doussal, 1992, Phys. Rev. Lett. 69, 1209.   
Rammal, R., 1985, J. Physique 46, 1345.   
Recher, P., B. Trauzettel, Y. M. Blaner, C. W. J. Beenakker, and A. F. Morpurgo, 2007, Phys. Rev. B 76, 235404.   
Reich, S., J. Maultzsch, C. Thomsen, and P. Ordejón, 2002, Phys. Rev. B 66, 035412.   
Robinson, J. P., and H. Schomerus, 2007, Phys. Rev. B 76, 115430.   
Rollings, E., G.-H. Gweon, S. Y. Zhou, B. S. Mun, J. L. McChesney, B. S. Hussain, A. V. Fedorov, P. N. First, W. A. de Heer, and A. Lanzara, 2005, J. Phys. Chem. Sol. 67, 2172.   
Rong, Z. Y., and P. Kuiper, 1993, Phys. Rev. B 48, 17427.   
Rosenstein, B., B. J. Warr, and S. H. Park, 1989, Phys. Rev. Lett. 62, 1433.   
Rosenstein, B., B. J. Warr, and S. H. Park, 1991, Phys. Rep. 205, 59.   
Russo, S., J. B. Oostinga, D. Wehenkel, H. B. Heersche, S. S. Sobhani, L. M. K. Vandersypen, and A. F. Morpurgo, 2007, eprint arXiv:0711.1508.   
Rutter, G. M., J. N. Crain, N. P. Guisinger, T. Li, P. N. First, and J. A. Stroscio, 2007, Science 317, 219.   
Rycerz, A., J. Tworzydlo, and C. W. J. Beenakker, 2007, Nature Physics 3, 172.   
Rydberg, H., M. Dion, N. Jacobson, E. Schröder, P. Hyldgaard, S. I. Simak, D. C. Langreth, and B. I. Lundqvist, 2003, Phys. Rev. Lett. 91, 126402.   
Ryu, S., C. Mudry, H. Obuse, and A. Furusaki, 2007, Phys. Rev. Lett. 99, 116601.   
Sabio, J., C. Seoanez, S. Fratini, F. Guinea, A. H. Castro Neto, and F. Sols, 2007, eprint arXiv:0712.222.   
Sadowski, M. L., G. Martinez, M. Potemski, C. Berger, and W. A. de Heer, 2006, Phys. Rev. Lett. 97, 266405.   
Safran, S. A., 1984, Phys. Rev. B 30, 421.   
Safran, S. A., and F. J. DiSalvo, 1979, Phys. Rev. B 20, 4889.   
Saha, S. K., U. V. Waghmare, H. R. Krishnamurth, and A. K. Sood, 2007, eprint cond-mat/0702627.   
Saito, R., G. Dresselhaus, and M. S. Dresselhaus, 1998, Physical properties of carbon nanotubes (Imperial College Press, London).   
Saito, R., M. Fujita, G. Dresselhaus, and M. S. Dresselhaus, 1992a, Appl. Phys. Lett. 60, 2204.   
Saito, R., M. Fujita, G. Dresselhaus, and M. S. Dresselhaus, 1992b, Phys. Rev. B 46, 1804.   
San-Jose, P., E. Prada, and D. Golubev, 2007, Phys. Rev. B 76, 195445.   
Saremi, S., 2007, Phys. Rev. B 76, 184430.   
Sarma, S. D., E. H. Hwang, and W. K. Tse, 2007, Phys. Rev. B 75, 121406.   
Schakel, A. M. J., 1991, Phys. Rev. D 43, 1428.   
Schedin, F., A. K. Geim, S. V. Morozov, D. Jiang, E. H. Hill, P. Blake, and K. S. Novoselov, 2007, Nature Materials 6, 652.   
Schomerus, H., 2007, Phys. Rev. B 76, 045433.   
Schroeder, P. R., M. S. Dresselhaus, and A. Javan, 1968, Phys. Rev. Lett. 20, 1292.   
Semenoff, G. W., 1984, Phys. Rev. Lett. 53, 2449.   
Sengupta, K., and G. Baskaran, 2007, Phys. Rev. B 77, 05417   
Seoanez, C., F. Guinea, and A. H. Castro Neto, 2007, Phys. Rev. B 76, 125427.   
Shankar, R., 1994, Rev. Mod. Phys. 66, 129.   
Sharma, M. P., L. G. Johnson, and J. W. McClure, 1974, Phys. Rev. B 9, 2467.   
Shelton, J. C., H. R. Patil, and J. M. Blakely, 1974, Surf. Sci. 43, 493.   
Sheng, D. N., L. Sheng, and Z. Y. Wen, 2006, Phys. Rev. B 73, 2333406.   
Sheng, L., D. N. Sheng, F. D. M. Haldane, and L. Balents, 2007, Phys. Rev. Lett. 99, 196802.   
Shklovskii, B. I., 2007, eprint arXiv:0706.4425.   
Shon, N. H., and T. Ando, 1998, J. Phys. Soc. Jpn. 67, 2421.   
Shung, K. W. K., 1986a, Phys. Rev. B 34, 979.   
Shug, K. W.K. 198, Phys. Rev. B3,   
Shytov, A. V., M. I. Katsnelson, and L. S. Levitov, 2007, Phys. Rev. Lett. 99, 236801.   
Silvestrov, P. G., and K. B. Efetov, 2007, Phys. Rev. Lett. 98, 016802.   
Sinitsyna, O. V., and I. V. Yaminsky, 2006, Russ. Chem. Rev. 75, 23.   
Skrypnyk, Y. V., and V. M. Loktev, 2006, Phys. Rev. B 73, 241402 (R).   
Skrypnyk, Y. V., and V. M. Loktev, 2007, Phys. Rev. B 75, 245401.   
Slonczewski, J. C., and P. R. Weiss, 1958, Phys. Rev. 109, 272.   
Snyman, I., and C. W. J. Beenakker, 2007, Phys. Rev. B 75, 045322.   
Sols, F., F. Guinea, and A. H. Castro Neto, 2007, Phys. Rev. Lett. 99, 166803.   
Son, D. T., 2007, Phys. Rev. B 75, 235423.   
Son, Y.-W., M. L. Cohen, and S. G. Louie, 2006a, Phys. Rev. Lett. 97, 216803.   
Son, Y.-W., M. L. Cohen, and S. G. Louie, 2006b, Nature 444, 347.   
Sorella, S., and E. Tosatti, 1992, Europhys. Lett. 19, 699.   
Soule, D. E., J. W. McClure, and L. B. Smith, 1964, Phys. Rev. 134, A453.   
Spataru, C. D., M. A. Cazalilla, A. Rubio, L. X. Benedict, P. M. Echenique, and S. G. Louie, 2001, Phys. Rev. Lett. 87, 246405.   
Spry, W. J., and P. M. Scherer, 1960, Phys. Rev. 120, 826.   
Stauber, T., F. Guinea, and M. A. H. Vozmediano, 2005, Phys. Rev. B 71, 041406.   
Stauber, T., N. M. R. Peres, F. Guinea, and A. H. Castro Neto, 2007, Phys. Rev. B 75, 115425.   
Stephan, O., P. M. Ajayan, C. Colliex, P. Redlich, J. M. Lambert, P. Bernier, and P. Lefin, 1994, Science 266, 1683.   
Stolyarova, E., K. T. Rim, S. Ryu, J. Maultzsch, P. Kim, L. E. Brus, T. F. Heinz, M. S. Hybertsen, and G. W. Flynn, 2007, Proc. Natl. Acad. Sci. USA 104, 9209.   
Stone, M., 1992, Quantum Hall Effect (World Scientific, Singapore).   
Su, W. P., J. R. Schrieffer, and A. J. Heeger, 1979, Phys. Rev. Lett. 42, 1698.   
Su, W. P., J. R. Schrieffer, and A. J. Heeger, 1980, Phys. Rev. B 22, 2099.   
Sugawara, K., T. Sato, S. Souma, T. Takahashi, and H. Suematsu, 2007, Phys. Rev. Lett. 98, 036801.   
Suzuura, H., and T. Ando, 2002a, Phys. Rev. Lett. 89, 266603.   
Suzuura, H., and T. Ando, 2002b, Phys. Rev. B 65, 235412.   
Swain, P. S., and D. Andelman, 1999, Langmuir 15, 8902.   
Tanuma, S., and H. Kamimura, 1985, Graphite intercalation compounds: Progress of Research in Japan (World Scientific, PA).   
Tersoff, J., 1992, Phys. Rev. B 46, 15546.   
Tikhonenko, F. V., D. W. Horsell, R. V. Gorbachev, and A. K. Savchenko, 2007, eprint cond-mat/0707.0140.   
Titov, M., 2007, Europhys. Lett. 79, 17004.   
Titov, M., and C. Beenakker, 2006, Phys. Rev. B 74, 041401.   
Tománek, D., S. G. Louie, H. J. Mamin, D. W. Abraham, R. E. Thomson, E. Ganz, and J. Clarke, 1987, Phys. Rev. B 35, 7790.   
Tombros, N., C. Jozsa, M. Popinciuc, H. T. Jonkman, and J. van Wees, 2007, Nature 448, 571.   
Trushin, M., and J. Schliemann, 2007, Phys. Rev. Lett. 99, 216602.   
Tu, Z.-C., and Z.-C. Ou-Yang, 2002, Phys. Rev. B 65, 233407.   
Tworzydlo, J., I. Snyman, A. R. Akhmerov, and C. W. J. Beenakker, 2007, Phys. Rev. B 76, 035411.   
Tworzydlo, J., B. Trauzettel, M. Titov, A. Rycerz, and C. W. J. Beenakker, 2006, Phys. Rev. Lett. 96, 246802.   
Uchoa, B., and A. H. Castro Neto, 2007, Phys. Rev. Lett. 98, 146801.   
Uchoa, B., C.-Y. Lin, and A. H. Castro Neto, 2007, Phys. Rev. B 77, 035420.   
Vafek, O., 2006, Phys. Rev. Lett. 97, 266406.   
Vafek, O., 2007, Phys. Rev. Lett. 98, 216401.   
Varchon, F., R. Feng, J. Hass, X. Li, B. N. Nguyen, C. Naud, P. Mallet, J. Y. Veuillen, C. Berger, E. H. Conrad, and L. Magaud, 2007, Phys. Rev. Lett. 99, 126805.   
Vázquez de Parga, A. L., F. Calleja, B. Borca, M. C. G. Passeggi Jr, J. J. Hinarejo, F. Guinea, and R. Miranda, 2007, Periodically rippled graphene: growth and spatially resolved electronic structure, eprint arXiv:0709.0360.   
Visscher, P. B., and L. M. Falicov, 1971, Phys. Rev. B 3, 2541.   
Vozmediano, M. A. H., M. P. López-Sancho, and F. Guinea, 2002, Phys. Rev. Lett. 89, 166401.   
Vozmediano, M. A. H., M. P. López-Sancho, and F. Guinea, 2003, Phys. Rev. B 68, 195122.   
Vozmediano, M. A. H., M. P. López-Sancho, T. Stauber, and F. Guinea, 2005, Phys. Rev. B 72, 155121.   
Wakabayashi, K., M. Fujita, H. Ajiki, and M. Sigrist, 1999, Phys. Rev. B 59, 8271.   
Wakayabashi, K., and M. Sigrist, 2000, Phys. Rev. Lett. 84, 3390.   
Wallace, P. R., 1947, Phys. Rev. 71, 622.   
Wang, H., D. N. Sheng, L. Sheng, and F. D. M. Haldane, 2007, eprint cond-mat/0708.0382.   
Wang, X.-F., and T. Chakraborty, 2007a, Phys. Rev. B 75, 033408.   
Wang, X.-F., and T. Chakraborty, 2007b, Phys. Rev. B 75, 041404(R).   
Wehling, T. O., K. S. Novoselov, S. V. Morozov, E. E. Vdovin, M. I. Katsnelson, A. K. Geim, and A. I. Lichtenstein, 2007, eprint cond-mat/0703390.   
Williams, J. R., L. DiCarlo, and C. M. Marcus, 2007, Science 317, 638.   
Williamson, S. J., S. Foner, and M. S. Dresselhaus, 1965, Phys. Rev. 140, A1429.   
Wirtz, L., and A. Rubio, 2004, Sol. Stat. Comm. 131, 141.   
Wu, X., X. Li, Z. Song, C. Berger, and W. A. de Heer, 2007, Phys. Rev. Lett. 98, 136801.   
Wunsch, B., T. Stauber, F. Sols, and F. Guinea, 2006, New J. Phys. 8, 318.   
Xin, Z., Z. Jianjun, and O.-Y. Zhong-can, 2000, Phys. Rev. B 62, 13692.   
Xu, S., J. Cao, C. C. Miller, D. A. Mantell, R. J. D. Miller, and Y. Gao, 1996, Phys. Rev. Lett. 76, 483.   
Yan, J., Y. Zhang, P. Kim, and A. Pinczuk, 2007, Phys. Rev. Lett. 98, 166802.   
Yang, L., M. V. Cohen, and S. G. Louie, 2007a, Nano Lett. 7, 3112.   
Yang, L., C.-H. Park, Y.-W. Son, M. L. Cohen, and S. G. Louie, 2007b, Phys. Rev. Lett. 99, 186801.   
Yang, X., and C. Nayak, 2002, Phys. Rev. B 65, 064523.   
Yao, Y., F. Ye, X.-L. Qi, S.-C. Zhang, and Z. Fang, 2007, Phys. Rev. B 75, 041401.   
Zarea, M., and N. Sandler, 2007, Phys. Rev. Lett. 99, 256804.   
Zhang, C.-H., and Y. N. Joglekar, 2007, Phys. Rev. B 75, 245414.   
Zhang, L. M., and M. M. Fogler, 2007, eprint arXiv:0708.0892.   
Zhang, Y., Z. Jiang, J. P. Small, M. S. Purewal, Y.-W. Tan, M. Fazlollahi, J. D. Chudow, J. A. Jaszczak, H. L. Stormer, and P. Kim, 2006, Phys. Rev. Lett. 96, 136806.   
Zhang, Y., Y.-W. Tan, H. L. Stormer, and P. Kim, 2005, Nature 438, 201.   
Zhong-can, O.-Y., Z.-B. Su, and C.-L. Wang, 1997w, Phys. Rev. Lett. 78, 4055.   
Zhou, S., G.-H. Gweon, and A. Lanzara, 2006a, Annals of Physics 321, 1730.   
Zhou, S. Y., G.-H. Gweon, A. V. Fedorov, P. N. First, W. A. de Heer, D.-H. Lee, F. Guinea, A. H. Castro Neto, and A. Lanzara, 2007, Nature Materials 6, 770.   
Zhou, S. Y., G.-H. Gweon, J. Graf, A. V. Fedorov, C. D. Spataru, R. D. Diehl, Y. Kopelevich, D.-H. Lee, S. G. Louie, and A. Lanzara, 2006b, Nature Physics 2, 595.   
Zhou, S. Y., G.-H. Gweon, J. Graf, A. V. Fedorov, C. D. Spataru, R. D. Diehl, Y. Kopelevich, D.-H. Lee, S. G. Louie, and A. Lanzara, 2006c, Nature Physics 2, 595.   
Ziegler, K., 1998, Phys. Rev. Lett. 80, 3113.   
Ziman, J. M., 1972, Principles of the Theory of Solids (Cambridge).