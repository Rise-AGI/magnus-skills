# Surface plasmon polariton propagation around bends at a metal-dielectric interface

Keisuke Hasegawa, Jens U. Nöckel and Miriam Deutsch Oregon Center for Optics,

1274 University of Oregon, Eugene, OR 97403-1274 http://oco.uoregon.edu/

Published in Appl. Phys. Lett. 84, 1835 (2004)

We analyze theoretically the propagation of surface plasmon polaritons about a metallic corner with a finite bend radius, using a one-dimensional model analogous to the scattering from a finite-depth potential well. We obtain expressions for the energy reflection and transmission coefficients in the short wavelength limit, as well as an upper bound for the transmittance. In certain cases we find that propagation on non-planar interfaces may result in lower losses than on flat surfaces, contrary to expectation. In addition, we also find that the maximum transmittance depends non-monotonously on the bend radius, allowing increased transmission with decreasing radius.

Structured materials which allow nanoscale control of light are necessary for achieving compact integrated photonic devices. While the size of standard optical components and beams is typically set by the diffraction limit, low dimensional excitations such as surfaceplasmon polaritons may be confined to dimensions much smaller than the wavelength of light. Surface-plasmon polaritons (SPPs), coupled modes of plasmons and photons, are excited when visible electromagnetic (EM) radiation couples into surface guided modes at metal-dielectric interfaces [1, 2]. When propagating along flat interfaces, these are essentially two-dimensional (2D) waves, with an EM field intensity which peaks at the interface and decays exponentially into the two adjoining media.

Recently, SPP waveguiding and bending in nano-patterned metallic films were studied [3]. Alternately, it was shown that EM energy may be efficiently transported by near field coupling in plasmon waveguides comprised of ordered arrays of metal nanoparticles [4]. Optical elements such as linear waveguides [5], mirrors, beamsplitters and interferometers [6] were recently demonstrated.

Interestingly, while significant progress has been made in understanding SPP propagation in nano-structures, certain fundamental issues pertaining to their guiding on smooth metallic films remain unknown. In particular, quantifying guiding and energy losses in SPPs propagating around bends in metal-dielectric interfaces is of great importance, as it should set a limit on feature size in certain plasmonic-circuit devices. Previously, the problems of refraction [7] and reflection of SPPs [8] at interfaces have been addressed in this context. In this Letter we present a study of the efficiency of SPP propagation at a curved metal-dielectric interface in the short wavelength limit.

![](images/d171df3e9531a9ba9264f12d4aff589bd7d239bbc21dae7b4877cc9951c54c0a.jpg)  
Figure 1: Cross-sectional (a) and top view (b) of the SPP intensity; (b) also shows a drawing of the geometry. A metallic corner characterized by a dielectric constant $\epsilon _ { i }$ and $\mathrm { R e } [ \epsilon _ { i } ] < 0$ has a bend angle $\theta > 9 0 ^ { \circ }$ and a finite bend radius $R$ . The bend is confined to the region of space shown, with the center of curvature at the origin. The rest of space is occupied by a dielectric with $\epsilon _ { o }$ . Axes $x _ { 1 }$ and $x _ { 2 }$ extend along the boundaries between Regions I and II and Regions II and III, respectively. In the $x _ { 1 }$ $x _ { 2 }$ plane, Regions I and III are semiinfinite. The system is also infinite in extent along the entire $z$ axis. In (a) we illustrate the single-mode approximation developed in the text: the field profiles in Regions I (solid line) and II (dashed), calculated for $\omega R / c = 8 0 0$ , are well matched. In (b), the intensity is overlayed in grayscale, showing the overlap with the metal (dielectric constant $\epsilon _ { i }$ ) and the outer dielectric $\left( \epsilon _ { o } \right)$ . Arrows indicate incident and reflected fields in Region I, and transmitted field in Region III.

The geometry we study is that of propagation about a rounded edge, as shown in Fig. 1. SPPs are incident along the interface from Region I, and propagate counterclockwise through the bend in Region II, in the direction of Region III. We calculate the energyreflection and transmission coefficients, as well as bend-induced radiation losses. For lossy metals absorption losses are also evaluated.

Our approach exploits known expressions for the SPP fields in each region and matches them at the two ends of the bend. The procedure differs from related numerical techniques [9] in that we consider the SPP itself as the incident wave, with the goal of manipulating it as a well-defined quasi-particle in a non-trivial geometry. Favorable conditions for this will be seen to emerge in the short-wavelength limit, and we use this to arrive at analytic expressions.

The solutions for SPP propagation on an infinite flat surface and on cylindrical surfaces are known analytically [2]. On a flat interface, SPPs at frequency $\omega$ are two dimensional waves, decaying exponentially into the two adjoining media, with decay constants $\gamma _ { i } =$ $- \omega \epsilon _ { i } \sqrt { - 1 / ( \epsilon _ { i } + \epsilon _ { o } ) } / c$ in the metal and $\gamma _ { o } = \omega \epsilon _ { o } \sqrt { - 1 / ( \epsilon _ { i } + \epsilon _ { o } ) } / c$ in the dielectric. In the limit $\mathrm { R e } [ \gamma _ { i } ] R \gg 1$ the interference of SPPs in Regions I and II is negligible, allowing us to use the infinite flat surface 2D solution in these regions.

We construct the solution in Region II using the known solutions for SPPs propagating around the perimeter of an infinitely long cylindrical metal rod of radius $R$ [2]. Here, the magnetic field is given by

$$
\mathbf { B } = - i \hat { \mathbf { z } } \sum _ { \{ n \} } \Bigl [ A _ { n } ^ { + } e ^ { + i n \phi } - A _ { n } ^ { - } e ^ { - i n \phi } \Bigr ] \sqrt { \epsilon _ { i } } J _ { n } ( k _ { i } r ) e ^ { - i \omega t }
$$

where $k _ { i } = \omega \sqrt { \epsilon _ { i } } / c$ and $J _ { n }$ is the Bessel function. The set $\{ n \}$ is determined by the metal boundary matching equation,

$$
\frac { 1 } { k _ { i } } \frac { J _ { n } ^ { \prime } ( k _ { i } R ) } { J _ { n } ( k _ { i } R ) } = \frac { 1 } { k _ { o } } \frac { H _ { n } ^ { ( 1 ) ^ { \prime } } ( k _ { o } R ) } { H _ { n } ^ { ( 1 ) } ( k _ { o } R ) }
$$

where $k _ { o } = \omega \sqrt { \epsilon _ { o } } / c$ , $H _ { n } ^ { ( 1 ) }$ is the Hankel function of the first kind, and the prime denotes differentiation with respect to the argument. Assuming $\omega$ real, one finds $n$ to be complex, as a consequence of radiation loss and absorption in the bend. Since the wave depends on $\phi$ as $\exp [ \pm i n \phi ] = \exp \left\lfloor \pm \left( i \mathrm { R e } [ n ] - \mathrm { I m } [ n ] \right) \phi \right\rfloor$ , where $\phi$ is measured from the $x _ { 1 }$ axis, only solutions with $\mathrm { I m } [ n ] \ge 0$ are physical for damped propagation.

Solving exactly for the transmission and reflection coefficients requires matching an infinite number of solutions at the boundaries along the $x _ { 1 }$ and $x _ { 2 }$ -axes separately. However, it is possible to render this problem tractable by a few simple approximations. Noting that the incident SPP carries momentum proportional to $k = \omega \sqrt { \epsilon _ { o } \epsilon _ { i } / ( \epsilon _ { o } + \epsilon _ { i } ) } / c$ , in the short wavelength limit $\omega R / c \gg 1$ its angular momentum with respect to the origin is approximately equal to $\mathrm { R e } [ k ] R$ . In Region II the solution has angular momentum proportional to the various $n$ -values. Conservation of angular momentum dictates that the incident SPP couple predominantly to that cylindrical mode with $n$ closest in value to $k R$ . Therefore it is necessary to consider only a single term of the expansion. Formally, this is shown by examining the set $\{ n \}$ and noting that it contains an element $m$ which minimizes the mismatch between the field profiles perpendicular to the surface. The role of angular momentum conservation in this matching problem is analogous to that of tangential momentum conservation in refraction at a dielectric interface. We call the clockwise and counterclockwise modes corresponding to $m$ the fundamental modes. In the short wavelength limit of a fundamental mode $n = m \approx k R$ , and the decay rate is $\sqrt { m ^ { 2 } / R ^ { 2 } - k _ { i } ^ { 2 } } \approx \gamma _ { i }$ .Thus, $J _ { m } ( k _ { i } r ) \sim \exp [ \gamma _ { i } r ]$ near the interface, identical to the fields in the metal in Regions I and III. Similarly, in the dielectric $H _ { m } ^ { ( 1 ) } ( k _ { o } r ) \sim \exp [ - \gamma _ { o } r ]$ [10].

The modes $n \neq m$ have decay rates not as close to $\mathrm { R e } [ \gamma _ { i } ]$ as that of the fundamental mode's. For this reason, it is possible to assume that in the short wavelength limit the incident SPPs couple predominantly to the fundamental modes and ignore other mode coupling. In order to satisfy the standard Maxwell boundary conditions it is therefore necessary to match only a small number of solutions at a single point on each axis, at a distance $R$ from the origin. The boundary conditions are thus also satisfied approximately over the entire extent of the axes. As can be seen from Fig. 1, the mode mismatch at the boundaries may be very slight. The problem has now essentially become one dimensional (1D), analogous to scattering from a 1D finite potential well [11]. Since the allowed $m$ -values are always complex, bound-state solutions in the well do not exist. This distinguishes the SPP on a bent surface from waveguide bends enclosed on all sides by infinite potential walls [12].

Applying the appropriate boundary conditions to the fields at the two boundaries results in the familiar expression for the transmittance

$$
\mathsf { T } = \left| \frac { 4 m k R } { - e ^ { i m \theta } ( m - k R ) ^ { 2 } + e ^ { - i m \theta } ( m + k R ) ^ { 2 } } \right| ^ { 2 } .
$$

When the losses in the metal are accounted for, $\mathrm { I m } [ m ]$ increases with $R$ , such that when $\mathrm { I m } [ m ] \theta \gg 1$ the transmittance becomes

$$
\mathsf { T } \approx 1 6 \frac { | m k R | ^ { 2 } } { | m + k R | ^ { 4 } } e ^ { - 2 \mathrm { I m } [ m ] \theta } .
$$

The reflectance R may be obtained in a similar manner. In the limit $\omega R / c \to \infty$ these expressions become exact.

For lossless metals, the bend-induced radiation losses are simply given by $\mathsf { A } \equiv 1 - \mathsf { T } - \mathsf { R }$ . Accounting for absorption in the metal we find that A now includes both radiation and absorption losses. We extract the radiation losses by integrating the Poynting vector $\mathbf { s }$ for unit incident flux in Region II at $r  \infty$ :

$$
\mathsf { P } \equiv \int _ { \phi _ { 0 } } ^ { \theta + \phi _ { 0 } } \mathbf { S } \cdot \hat { \mathbf { r } } r d \phi .
$$

Since the radiation carries angular momentum with respect to the origin, the energy radiated into the far-field from the surface at $\phi = 0$ propagates at an angle $\phi _ { 0 }$ with respect to the $x _ { 1 }$ axis, setting the lower integration limit in (4). In the short-wavelength limit only the amplitude of the forward-propagating mode is significant, therefore the radiation losses are well approximated by integrating only this mode. To obtain $\phi _ { 0 }$ we use a stationary phase approximation. The position-dependent phase is $\Phi = k _ { o } r + \mathrm { R e } \lfloor m \rfloor \phi$ , and the vector normal to a surface of constant phase is $\mathbf { v } ( \mathbf { r } ) = \nabla \Phi = k _ { o } \hat { \mathbf { r } } + \mathrm { R e } [ m ] / r \hat { \mathbf { e } } _ { \phi }$ .

![](images/8820f465ac6eab0aae2a24e7cc5c46d1c9e72ba99df55a85d042670e1dc94ef9.jpg)  
Figure 2: The upper bound for the transmittance, $\mathsf { T } _ { u }$ , plotted for a silver-air interface with $\theta = 9 0 ^ { \circ }$ , as a function of bend radius $R$ for wavelengths of $\lambda = 5 0 0 \mathrm { { n m } }$ (dashed-dotted), $\lambda = 6 0 0 \mathrm { { n m } }$ (dashed) and $\lambda = 7 0 0 \mathrm { { n m } }$ (solid). Inset: $\mathsf { T } _ { u }$ in grayscale as a function of $R$ and $\lambda$ .

The change in angle as the wave propagates a radial distance $\delta r$ is $\delta \phi = \mathrm { R e } [ m ] / ( k _ { o } r ^ { 2 } ) \delta r$ giving $\begin{array} { r } { \phi _ { 0 } = \int _ { R } ^ { \infty } \mathrm { R e } [ m ] / ( k _ { o } r ^ { 2 } ) d r = \mathrm { R e } [ m ] / ( k _ { o } R ) } \end{array}$ .

We have carried out calculations for typical values of silver ( $\epsilon _ { i } = - 1 5 + i 0 . 5 ,$ in air $\epsilon _ { o } ~ = ~ 1$ ) when $\omega R / c = 8 0 0$ and $\theta \ : = \ : 9 0 ^ { \circ }$ . Ignoring the losses in the metal we find ${ \sf T } = 0 . 9 9 7$ , $\mathsf { R } = 1 . 1 9 \times 1 0 ^ { - 8 }$ and $\mathsf { P } \approx 0 . 0 0 3$ . When the losses are included the results change drastically to ${ \sf T } = 0 . 0 5 1 6$ , $\mathsf { R } = 1 . 1 8 \times 1 0 ^ { - 6 }$ and $\mathsf { P } \approx 0 . 0 0 2 8 2$ , indicating that most of the energy is lost to absorption in the metal. Comparing the latter overall absorption and radiation losses to the energy absorbed when SPPs propagate the equivalent arc distance on a flat surface, we find that propagation on a non-planar interface may result in lower losses. We explain this counterintuitive result using an analogous picture of semi-classical motion under an effective potential in a central potential field. In the short wavelength and large angular momentum limit the SPP fields propagating on the curved interface sample less of the metal volume than that available when propagating on a flat interface, hence the reduced absorption.

We evaluate the accuracy of our result by examining the coupling efficiency $\Delta$ of a single mode on a flat interface to a fundamental mode $m$ . We define this by

$$
\Delta ^ { 2 } \equiv \frac { \int _ { R } ^ { R + \eta \gamma _ { o } ^ { - 1 } } \left| \frac { H _ { m } ^ { ( 1 ) } ( k _ { o } r ) } { H _ { m } ^ { ( 1 ) } ( k _ { o } R ) } - \exp [ - \gamma _ { o } ( r - R ) ] \right| ^ { 2 } d r } { \int _ { R } ^ { R + \eta \gamma _ { o } ^ { - 1 } } \big | \exp [ - \gamma _ { o } ( r - R ) ] \big | ^ { 2 } d r }
$$

where $\eta = O ( 1 )$ . We find that the condition $\Delta \ll 1$ constitutes a stricter criterion for the validity of our approximation. When the latter holds, the incident SPP couples predominantly to the fundamental modes, making the approach described above selfconsistent. For example, when $\eta = 3$ , $\epsilon _ { i } = - 1 5$ , and $\epsilon _ { o } = 1$ , $\Delta ^ { 2 } = 0 . 0 0 2$ for $\omega R / c = 8 0 0$ , rendering our result applicable. On the other hand, for $\omega R / c = 1 0 0$ we obtain $\Delta ^ { 2 } = 0 . 3$ , signifying that the expression is not reliable because the coupling to modes $n$ other than the fundamental can no longer be neglected. In this regime a more physical quantity is the upper bound for the transmittance, given from (3) by

$$
\mathsf { T } _ { u } = \exp \bigl ( - 2 \mathrm { I m } [ m ] \theta \bigr ) .
$$

Here we neglect reflections at the boundaries between the different regions, thus excluding interference with the counter-propagating mode in Region II. To understand why this is an upper bound, recall that in the wavelength range of interest, where the metal is not very lossy, $\mathrm { I m } [ n ] > \mathrm { I m } [ m ]$ . Since the wave depends on $n$ as $\exp [ \pm i n \phi ]$ , modes with large $\mathrm { I m } [ n ]$ decay rapidly. Thus, the transmission in the presence of coupling to non-fundamental modes does not exceed $\textsf { T } _ { u }$ , and the latter is a true upper bound. Fig. 2 is a plot of $\mathsf { T } _ { u }$ . A peak is clearly visible, moving to higher values of $R$ as the wavelength increases. To the right of the peak, at large radii of curvature absorption losses in the metal dominate, and the maximum transmittance decreases with increasing radius. To the left of the peak radiation due to the high curvature is the dominant loss mechanism, leading to a rapid drop in $\mathsf { T } _ { u }$ . At very high curvature ( $R \leq 1 0 \mu \mathrm { m }$ ) there is a change in trend, and ${ \sf T } _ { u }$ starts to increase with decreasing $R$ (see Inset.) When calculating the radiation loss per arclength, we find that for this range of radii it increases slower than elsewhere, allowing $\mathsf { T } _ { u }$ to increase even as $R$ attains very small values. This anomalous behavior can be observed for all wavelengths, and is independent of the dispersion in the metal.

The formalism developed here may also be used to analyze the complementary reversed geometry, where the metal occupies three quadrants in space, and the SPPs propagate around a dielectric void in it. Surprisingly, in this case we find that in the single mode approximation SPP propagation around the bend is non-radiative. Separate work will address radiation and absorption processes in this system in greater detail.

In summary, we have analyzed the scattering of SPPs at a curved metal-dielectric interface in the short wavelength limit. Utilizing an analogy to a quantum mechanical 1D finite square well we obtained the energy transmission and reflection coefficients. Interestingly, propagation on a curved interface may result in lower losses than at a flat metallic surface, due to the unique field distributions which arise in our system. An expression for an upper bound on the transmittance was also obtained, showing that at high curvature radiation is the main loss mechanism, while at low curvature material losses dominate. An unexpected behavior where the maximum transmittance increases with curvature was also observed. We explain this as an interplay between various loss rates in the system. These results shed new light on the mesoscopic behavior of SPPs, and should play an important role in the design and optimization of SPP devices. Future work will address SPP propagation in waveguides, splitters and interferometers.

J.U.N. acknowledges support from NSF Grant ECS-02-39332; K.H. and M.D. acknowl. edge support from NSF Grant DMR-02-39273 and ARO Grant DAAD19-02-1-0286.

# References

[1] V. M. Agranovich, D. L. Mills, eds., Surface Polaritons (North Holland, Amsterdam, 1982).   
[2] B. E. Sernelius, Surface Modes in Physics (Wiley, Berlin, 2001).   
[3] S. I. Bozhevolnyi, J. Erland, K. Leosson, P. M. W. Skovgaard, J. M. Hvam, Phys. Rev. Lett. 86, 3008 (2001), S. I. Bozhevolnyi, V. S. Volkov, K. Leosson, A. Boltasseva, Appl. Phys. Lett. 79, 1076 (2001).   
[4] M. Quinten, A. Leitner, J. R. Krenn, F. R. Aussenegg, Opt. Lett. 23, 1331 (1998).   
[5] S. A. Maier, P. G. Kik, H. A. Atwater, S. Meltzer, E. Harel, B. E. Koel, A. A. G. Requicha, Nature Materials 2, 229 (2003).   
[6] H. Ditlbacher, J. R. Krenn, G. Schider, A. Leitner, F. R. Aussenegg, Appl. Phys. Lett. 81, 1762 (2002).   
[7] G. I. Stegeman, A. A. Maradudin, T. S. Rahman, Phys. Rev. B 23, 2576 (1981).   
[8] R. F. Wallis and A. A. Maradudin, Appl. Phys. Lett. 42, 764 (1983); P. Dawson, F. de Fornel, J-P. Goudonnet, Phys. Rev. Lett. 72, 2927 (1994).   
[9] E. Moreno, D. Erni, C. Hafner, R. Vahldieck, J. Opt. Soc. Am. A 19, 101 (2002).   
[10] To show this we use the Bessel equation $r ^ { 2 } d ^ { 2 } f / d r ^ { 2 } + r d f / d r + ( k _ { i , o } ^ { 2 } r ^ { 2 } - n ^ { 2 } ) f = 0$ , to which $J _ { m } ( k _ { i } r )$ and $H _ { m } ^ { ( 1 ) } ( k _ { o } r )$ are solutions, respectively. With $r = R + x$ , in the limit $x \ll R$ this reduces to $R ^ { 2 } d ^ { 2 } f / d x ^ { 2 } + ( k _ { i , o } ^ { 2 } R ^ { 2 } - n ^ { 2 } ) f = 0$ , to which the solutions are exponentials.   
[11] A. Mekis, J. C. Chen, I. Kurland, S. Fan, P. R. Villeneuve, J. D. Joannopoulos, Phys. Rev. Lett. 77, 3787 (1996).   
[12] F. Sols and M. Macucci, Phys. Rev. B 41, 11887 (1990)