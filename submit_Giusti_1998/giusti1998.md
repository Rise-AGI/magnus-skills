# The QCD Chiral Condensate from the Lattice

L. Giusti

Scuola Normale Superiore, P.zza dei Cavalieri 7, I-56100 Pisa, Italy INFN-Sezione di Pisa, I-56100 Pisa, Italy

F. Rapuano

INFN-Sezione di Roma, c/o Dipartimento di Fisica, Università di Roma 'La Sapienza', P.le A. Moro 2, I-00185 Roma, Italy.

M. Talevi

Department of Physics & Astronomy, University of Edinburgh The King's Buildings, Edinburgh EH9 3JZ, UK

A. Vladikas

INFN-Sezione di Roma II, c/o Dipartimento di Fisica, Università di Roma 'Tor Vergata', Via della Ricerca Scientifica 1, I-00133 Roma, Italy.

# Abstract

We determine the chiral condensate from simulations of quenched lattice QCD with Wilson fermions. Our measurements have been obtained with high statistics at three values of the gauge coupling, corresponding to UV cutoffs in the range $2 - 4$ GeV. Several improvements have been made with respect to earlier lattice computations. The most important are the non-perturbative renormalization of the condensate, the use of the tree-level improved Clover action and the reduction of the systematic error due to uncertainties in the lattice calibration. Our result for the chiral condensate in the $\overline { { \mathrm { M S } } }$ scheme is

$\langle \bar { \psi } \psi \rangle ^ { \overline { { { \mathrm { M S } } } } } ( \mu = \ 2 \mathrm { G e V } ) \ = \ - \ 0 . 0 1 4 7 ( 8 ) ( 1 6 ) ( 1 2 ) \ \mathrm { G e V } ^ { 3 } = \ - \ [ 2 4 5 ( 4 ) ( 9 ) ( 7 ) ] ^ { 8 } ,$ MeV]3 where the first error is statistical, the second is due to the non-perturbative renormalization and the third due to the lattice calibration.

# 1 Introduction

QCD has a remarkably small mass gap. Its pattern is described by a nearly degenerate isospin multiplet containing pseudoscalar mesons, but not their partners of opposite parity (vector mesons). This implies that chiral symmetry is spontaneously broken in massless QCD. This should be signalled by the non-vanishing of some order parameter (for a given order parameter this is a sufficient but not necessary condition for symmetry breaking). The chiral condensate is such an order parameter. The standard expectation is that it does not vanish in the chiral limit. However, it cannot be a priori excluded (see ref. [1] and, more recently, [2]) that it might tend to zero in this limit, while chiral symmetry breaking is nevertheless realized through the non-vanishing of the pion decay constant. A more general option is that the chiral condensate has a small non-zero value [2].

The intrinsically non-perturbative nature of the chiral condensate implies that its determination from first principles is not straightforward. QCD sum rules and lattice computations have been used in order to estimate its value. The lattice regularization is particularly well suited for the computation of low-energy non-perturbative quantities from first principles. Although lattice measurements are affected by a series of approximations (finite configuration ensemble, finite volume and lattice spacing, non-zero quark mass, quenched approximation), they can be controlled and systematically improved, at least in principle. Two of these sources of error have the most negative impact on the credibility of the lattice determination of the chiral condensate, namely the quenched approximation and the chiral extrapolation of data computed at non-zero quark masses. Quenching is destined to stay with lattice computations for some time, since realistic simulations with light sea quark masses are still unaccessible to present-day computers. Quenched simulations with fairly light valence quarks are more accessible, but the zero quark limit is elusive for hadronic quantities 1.

Interestingly, the two approximations (quenching and chiral extrapolations) are interwoven by the presence of quenched chiral logarithms in many physical quantities, such as the chiral condensate and the pion mass (but not in the pion decay constant, in the isospin limit) [6]. There has been some debate recently on the evidence for quenched chiral logarithms in lattice results obtained with light quark masses [7]. However, for the bulk of quenched lattice results, obtained at the strange quark mass region, the power dependence on the quark mass dominates the logarithms. In such simulations many physical quantities known from experiment, including those related to the chiral condensate (e.g meson masses, decay constants), have been determined with a 5-15% accuracy (see for example ref. [8]). Also, the recent lattice determinations of the quark masses [9] agree with the sum rules predictions. The good agreement of many quenched results with experiment 2 , obtained in a quark mass region where chiral singularities can be ignored, implies that the regular dependence of these observables on the quark mass is adequately reproduced in the quenched approximation.

Given the above observations, we have performed an extensive quenched lattice computation of the chiral condensate, with Wilson fermions, using the data accumulated by the APE collaboration in the last few years. The valence quark masses used in these simulations are in the range of the strange flavour. Several important improvements have been possible compared to earlier computations:

• We have improved our statistics by computing the condensate on ensembles of $\mathcal { O } ( 1 0 0 )$ gauge field configurations. Most previous computations with Wilson fermions were performed on $\mathcal { O } ( 1 0 )$ configurations. • Reliable measurements have been performed at two values of the lattice gauge coupling $g _ { 0 } ^ { 2 }$ , corresponding to UV cutoffs of $\sim 2$ GeV and $\sim 3$ GeV, using both the Wilson and the tree-level SW-Clover[10,11] actions. Results obtained with the Wilson action at fixed UV cutoff (i.e. fixed lattice spacing $a$ ) suffer from $\mathcal O ( a )$ discretization errors; those obtained in the tree-level Clover formalism suffer from $\mathcal { O } ( a g _ { 0 } ^ { 2 } )$ errors. Thus we can study the influence of finite lattice spacing effects. We also present results at a third value of the gauge coupling (corresponding to an UV cutoff of $\sim 4$ GeV), which are however somewhat unreliable, due to the small physical size of the lattice. Different determinations of the chiral condensate (expressed in lattice units) have been implemented. They all derive from the same Ward Identity(WI), obtained from the lattice action, and are equivalent only in the chiral and continuum limits. They are implemented, however, at non-zero lattice spacing and quark masses and subsequently extrapolated to the chiral limit. Thus, by comparing the chiral condensate from different determinations, we can cross-check the reliability of the extrapolations. √ The systematic errors arising from the chiral extrapolation and the uncertainties in the determination of the lattice spacing are amplified in the standard computation of the chiral condensate. The chiral extrapolation amplifies the error because we extrapolate from data obtained in the strange quark region. The uncertainty in the determination of the lattice spacing is amplified because the condensate is a dimension-3 matrix element. We have obtained the lattice condensate in physical units in a way that avoids these

problems.

The chiral condensate measured on the lattice at a given fixed lattice spacing is a bare quantity which must be renormalized. In previous studies, its renormalization constant was calculated in lattice Perturbation Theory (PT) at 1-loop. We will show that the perturbatively renormalized chiral condensate suffers from large $\mathcal { O } ( g _ { 0 } ^ { 4 } )$ corrections, which render its determination unreliable (an identical conclusion was drawn for the quark mass in ref. [9]). In the present work, we have used non-perturbative (NP) estimates of the required renormalization constant, which are known in the Regularization Independent (RI) renormalization scheme (also known as the MOM scheme). This removes lattice PT completely from the renormalization procedure. We only use continuum PT in order to transform the results from the RI to the intrinsically perturbative $\mathrm { \overline { { M S } } }$ scheme. Also, RG running has been done in PT in order to express all results in the conventional scale of 2 GeV. All PT coefficients used in this work are at NNLO.

Our best estimate for the chiral condensate agrees with the earlier lattice determinations of refs. [12,13] with greatly improved statistical accuracy. For the first time we are also able to estimate the systematic uncertainties in a reliable way. Our result is also in excellent agreement with the one obtained from sum rules (see [14] and references therein).

We now outline the contents of this paper. In sect. 2 we fix our notation and give the basic definitions. In sect. 3 we discuss issues related to our renormalization strategy. In particular we review the most important properties of the NP renormalization of the lattice condensate (see Appendix A for details). We also review the connection between the RI and $\mathrm { \overline { { M S } } }$ renormalization schemes (details can be found in Appendix B). In sect. 4 we discuss how the chiral condensate can be determined from the lattice regularization with Wilson fermions. Chiral symmetry, explicitly broken by the Wilson fermionic action, is recovered in the continuum and chiral limits through the implementation of lattice WIs; see ref. [15]. With the aid of such a WI we obtain the proper lattice definition of the chiral condensate which, for Wilson fermions, is subject to a power subtraction. Moreover, from this WI we derive the well known Gell-Mann-Oakes-Renner (GMOR) relation [16] for the condensate. Some technicalities related to the WIs with Clover fermions are discussed in Appendix C. In sect. 5 we present our results and compare them to those of previous computations of the chiral condensate with Wilson fermions 3 . Our conclusions are summarized in sect. 6.

In this section we define the quantities of interest to this paper, in order to fix our notation. We work in the lattice regularization scheme proposed by Wilson [18]. The partition function of the theory is

$$
Z = \int { \mathcal D } \psi { \mathcal D } \bar { \psi } { \mathcal D } U _ { \mu } \exp \left( - S _ { g } - S _ { f } \right) ,
$$

where the gluonic action is given by

$$
S _ { g } = \frac { 6 } { g _ { 0 } ^ { 2 } } \sum _ { P } \left\{ 1 - \frac { 1 } { 6 } T r [ U _ { P } + U _ { P } ^ { \dagger } ] \right\}
$$

and the fermionic action by

$$
\begin{array} { l } { { S _ { f } = - a ^ { 4 } \displaystyle \sum _ { x , \mu } \displaystyle \frac { 1 } { 2 a } [ { \bar { \psi } } ( x ) ( 1 - \gamma _ { \mu } ) U _ { \mu } ( x ) \psi ( x + \mu ) + { \bar { \psi } } ( x + \mu ) ( 1 + \gamma _ { \mu } ) U _ { \mu } ^ { \dagger } ( x ) \psi ( x ) ] } } \\ { { \displaystyle ~ + a ^ { 4 } \sum _ { x } { \bar { \psi } } ( x ) ( M _ { 0 } + \frac { 4 } { a } ) \psi ( x ) ~ . } } \end{array}
$$

In standard notation, $a$ is the lattice spacing, $\psi ( x )$ is the quark field (flavour indices are implied), $U _ { P }$ is the Wilson plaquette, $U _ { \mu } ( x )$ is the lattice gauge link and $g _ { 0 }$ is the bare coupling constant. The diagonal bare mass matrix is denoted by $M _ { 0 }$ . In most of this work we will assume mass degeneracy; i.e. all elements of $M _ { 0 }$ are equal to $m _ { 0 }$ . As is well known, for Wilson fermions, the quark mass, besides a multiplicative renormalization, is also subject to an additive one. Defining $m = m _ { 0 } - m _ { C }$ , the renormalized is

$$
m _ { R } = Z _ { m } m = Z _ { m } [ m _ { 0 } - m _ { C } ]
$$

with $Z _ { m }$ the multiplicative renormalization constant. The chiral limit is then $m \to 0$ ; i.e. $m _ { 0 }  m _ { C }$ . At tree-level, $m _ { C } = - 4 / a$ .

We define the scalar and pseudoscalar densities as:

$$
\begin{array} { l } { { S ^ { f } ( x ) = \bar { \psi } ( x ) \displaystyle \frac { \lambda ^ { f } } { 2 } \psi ( x ) } } \\ { { P ^ { f } ( x ) = \bar { \psi } ( x ) \gamma _ { 5 } \displaystyle \frac { \lambda ^ { f } } { 2 } \psi ( x ) } } \end{array}
$$

and the vector and axial currents as

$$
\begin{array} { l } { { \displaystyle V _ { \mu } ^ { f } ( x ) = \bar { \psi } ( x ) \gamma _ { \mu } \frac { \lambda ^ { f } } { 2 } \psi ( x ) } } \\ { { \displaystyle A _ { \mu } ^ { f } ( x ) = \bar { \psi } ( x ) \gamma _ { \mu } \gamma _ { 5 } \frac { \lambda ^ { f } } { 2 } \psi ( x ) \ : , } } \end{array}
$$

where $\lambda ^ { f } / 2$ are the generators of the $S U ( N _ { f } )$ flavour group in the fundamental representation (e.g. Pauli matrices for $N _ { f } = 2$ , Gell-Mann matrices for $N _ { f } = 3$ with $f = 1 , \dots , N _ { f } ^ { 2 } - 1 )$ . These generators satisfy:

$$
\begin{array} { l } { { \displaystyle T r [ \lambda ^ { f } \lambda ^ { g } ] = 2 \delta ^ { f g } } } \\ { { \displaystyle \left[ \frac { \lambda ^ { f } } { 2 } , \frac { \lambda ^ { g } } { 2 } \right] = i f ^ { f g h } \frac { \lambda ^ { h } } { 2 } } } \\ { { \displaystyle \left\{ \frac { \lambda ^ { f } } { 2 } , \frac { \lambda ^ { g } } { 2 } \right\} = d ^ { f g h } \frac { \lambda ^ { h } } { 2 } + \frac { \delta ^ { f g } } { N _ { f } } \lambda ^ { 0 } ~ , } } \end{array}
$$

where $f ^ { f g h }$ are the $S U ( N _ { f } )$ structure constants, $d ^ { f g h }$ are symmetric coefficients and $\lambda ^ { 0 }$ represents the $N _ { f } \times N _ { f }$ identity matrix. We will also consider the corresponding singlet bilinear scalar operator

$$
S ^ { 0 } ( x ) = \bar { \psi } ( x ) \lambda ^ { 0 } \psi ( x ) = \bar { \psi } ( x ) \psi ( x ) .
$$

The operators defined above are to be understood as bare lattice quantities 4 . They are subject to multiplicative renormalization. The renormalization of a generic bilinear operator $O _ { \Gamma } ^ { f } = \bar { \psi } \Gamma ( \lambda ^ { f } / 2 ) \psi$ is given by ( $\Gamma$ stands for any Dirac matrix)

$$
\hat { \cal O } _ { \Gamma } ( g _ { R } ^ { 2 } , m _ { R } , \mu ) = \operatorname * { l i m } _ { a \to 0 } [ Z _ { \cal O } ( g _ { 0 } ^ { 2 } , a \mu ) { \cal O } _ { \Gamma } ( g _ { 0 } ^ { 2 } , m , a ) ] ,
$$

where $g _ { R }$ is the renormalized gauge coupling and $\mu$ the renormalization scale. Note that in the bare operator the bare mass $m _ { 0 }$ has been traded off for the more convenient (from the point of view of the chiral limit) subtracted mass $m$ . Renormalized quantities will be denoted in general by "hats"; e.g. the renormalized quark propagator is

$$
\hat { S } ( p ) = \operatorname * { l i m } _ { a  0 } [ Z _ { q } S ( p ) ] \ ,
$$

where $Z _ { q } ^ { 1 / 2 }$ is the quark field renormalization.

The properties of the renormalization constants of the operators of eqs. (5) and (6), obtained in refs. [15] (for a recent treatment see [19]), are as follows: The renormalization constants $Z _ { S }$ and $Z _ { P }$ of the scalar and pseudoscalar densities are logarithmically divergent; thus they depend on the coupling $g _ { 0 } ^ { 2 }$ and the renormalization scale in lattice units $a \mu$ . Those of the vector and axial currents ( $Z _ { V }$ and $Z _ { A }$ respectively) are finite; thus they depend only on the coupling $g _ { 0 } ^ { 2 }$ . A mass independent renormalization scheme is implied.

# 3 Renormalization scheme dependence of the chiral condensate

In this section we discuss some issues related to the choice of renormalization scheme for the chiral condensate. The renormalization procedure which gives a finite condensate starting from its bare (lattice) counterpart is postponed until the next section. Here we will concentrate on how to connect the condensate in the RI renormalization scheme, obtained non-perturbatively, to the one in the $\mathrm { \overline { { M S } } }$ scheme. This is necessary in view of the fact that the RI scheme is most convenient for lattice computations, as discussed below, whereas the $\mathrm { \overline { { M S } } }$ scheme is conventionally preferred when expressing "final" physical results.

Let us sketch step by step the renormalization procedure for the chiral condensate. The breaking of chiral symmetry by the Wilson term of the fermionic action $S _ { f }$ induces some subtleties in the renormalization of the lattice condensate. These have been dealt with in ref. [15] and will be fully exposed in sec. 4. Here we simply anticipate the basic result: by using lattice WIs and by requiring that in the continuum and chiral limit the nominal WIs are recovered, we obtain that the correct renormalization of the lattice chiral condensate is given by

$$
\langle \bar { \psi } \psi \rangle ( \mu ) = Z _ { P } ( \mu a ) \left[ \langle \bar { \psi } \psi \rangle ( a ) - b _ { 0 } \right] ,
$$

where $\langle \psi \psi \rangle ( a ) = \langle S ^ { 0 } \rangle$ stands for the bare chiral condensate, computed nonperturbatively on the lattice, and $b _ { 0 }$ is a cubicly diverging subtraction. Note that $Z _ { P }$ is the renormalization of the non-singlet pseudoscalar density. Its determination from 1-loop lattice perturbation theory (see refs. [2022] for the perturbative renormalization of quark bilinears) suffers from large uncertainties due to the presence of "tadpole" diagrams [23]. In order to avoid this problem non-perturbative renormalization techniques have been developed [5,12,24,25,40] consisting in computing the renormalization constants at fixed UV cutoff (ranging, in present-day simulations, in the region $a ^ { - 1 } \simeq 2 - 4$ GeV). This is a tradeoff between higher order corrections (present in PT) and lattice artifacts (present in all non-perturbative methods). The latter, however, can be systematically controlled by improving the action and operators in the spirit of

Symanzik's work [5,10,11,26]. Thus in principle non-perturbative evaluations of the renormalization constants are preferable. The finite renormalization constants ( $Z _ { V }$ , $Z _ { A }$ and the ratio $Z _ { P } / Z _ { S }$ ) can be obtained non-perturbatively with the help of WIs (see [5,12,15,19,24,40]). All renormalization constants (including the diverging ones $Z _ { S }$ and $Z _ { P }$ ) can be obtained non-perturbatively with the help of the more recent "Non-Perturbative" (NP) method of ref. [25]. Results of the NP renormalization constants of quark bilinear operators can be found in refs. [25,27,28].

In practice, the necessity of non-perturbative renormalization was first demonstrated in the measurement of several B-parameters of the $\Delta S \ = \ 2$ fourfermion operators; see refs. [29]. More recently, it has also been shown in ref. [9] that the discrepancy, observed in the past between quark masses extracted from the Vector and Axial WIs, was due to the very poor determination of $Z _ { P }$ in PT. Use of the NP value led to an excellent agreement between the two determinations. The same conclusions are drawn in the present work for the chiral condensate.

We will therefore use the NP value of $Z _ { P }$ in the present work. The underlying theory and method of computation has been fully exposed in [25]; it consists in applying the RI (or MOM) scheme on the lattice. The amputated Green function of the insertion of the operator $P$ in the quark propagator is computed at fixed external momentum in the deep Euclidean region. The Green function is subsequently projected by a "suitable projector" and the RI renormalization condition fixes the renormalized projected Green function be equal to its treelevel value at a given scale $\mu$ . This condition can then be solved for $Z _ { P }$ ; the result is a fully non-perturbative estimate of the renormalization constant at fixed UV cutoff. For this procedure to be reliable $\mu$ must satisfy the conditions $\mu \ll \mathcal { O } ( a )$ (in order to avoid discretization errors) and $\mu \gg \Lambda _ { Q C D }$ . The former bound is approximately satisfied in practical simulations; we will be using the renormalization constants of ref. [28] at the scale $\mu \simeq a ^ { - 1 }$ . The latter bound is necessary for two reasons. (i) The RI scheme imposes that the renormalized operator transform as an irreducible representation of the chiral group. This is true only at large scales, where chiral symmetry breaking effects are negligible. (ii) At large scales we avoid higher order corrections in the continuum perturbative expansion in which the RI- $\mathrm { \overline { { M S } } }$ matching coefficient is calculated; see below. Unlike $\mathrm { \overline { { M S } } }$ , the RI renormalization condition does not depend on the regularization chosen, hence its name, Regularization Independent (RI). The explicit formulae of the RI renormalization condition and a discussion of some subtleties related to wavefunction renormalization can be found in Appendix A.

Thus the chiral condensate of eq. (11), renormalized with the NP method outlined above at scale $\mu \simeq a ^ { - 1 }$ , is expressed in the RI scheme. In order to enable comparisons with other methods, we must convert our result to the $\mathrm { \overline { { M S } } }$ scheme 5 . This is expressed in terms of a perturbative finite matching coefficient

$$
\langle \bar { \psi } \psi \rangle ^ { \overline { { \mathrm { M S } } } } ( \mu ) = \Delta Z ^ { \overline { { \mathrm { M S , R I } } } } \langle \bar { \psi } \psi \rangle ^ { R I } ( \mu ) .
$$

The matching coefficient $\Delta Z ^ { \overline { { \mathrm { M S } } } , \mathrm { R I } }$ , known in PT to NNLO (see refs. [28] and [30]), is given in Appendix B. Note that the renormalization constants obtained in the RI scheme are in general gauge dependent (we work in the Landau gauge). This dependence is cancelled by the gauge dependence of ∆ZMS,RI.

In order to express our result at the conventional scale of $\mu ^ { \prime } = 2$ GeV, we use perturbative RG running

$$
\langle \bar { \psi } \psi \rangle ^ { \overline { { \mathrm { M S } } } } ( \mu ^ { \prime } ) = \frac { c _ { S } ^ { \overline { { \mathrm { M S } } } } ( \mu ^ { \prime } ) } { c _ { S } ^ { \overline { { \mathrm { M S } } } } ( \mu ) } \langle \bar { \psi } \psi \rangle ^ { \overline { { \mathrm { M S } } } } ( \mu ) ~ .
$$

Finally, it is convenient to express the chiral condensate in a scale independent way. We define the RGI chiral condensate as

$$
\langle \bar { \psi } \psi \rangle ^ { \mathrm { R G I } } = \frac { 1 } { c _ { S } ^ { \overline { { \mathrm { M S } } } } ( \mu ) } \langle \bar { \psi } \psi \rangle ^ { \overline { { \mathrm { M S } } } } ( \mu ) ~ .
$$

to h This quantity does not depend neither on the scale $c _ { S } ^ { \mathrm { M S } } ( \mu )$ is known to NNLO; see Ap- $\mu$ nor on the renormalpendix B.

# 4 The chiral condensate from lattice WIs with Wilson fermions

In this section we obtain the proper definition of the chiral condensate with Wilson fermions from lattice WIs. The subject of chiral symmetry breaking by the Wilson term and its restoration in the continuum limit with the aid of WIs is well understood (see the original works [15], a recent review [19] and a recent theoretical treatment of the subject [31]). Thus we only present in some detail those results essential to the understanding of the lattice determination of the chiral condensate. First we shall briefly digress to vector WIs in order to demonstrate that the quark mass renormalization is the inverse of the renormalization of the scalar density. Then we will obtain the proper definition of the lattice chiral condensate and its renormalization using axial WIs from which the Gell-MannOakes-Renner relation [16] will be derived on the lattice.

The local $S U _ { L } ( N _ { f } ) \times S U _ { R } ( N _ { f } )$ chiral transformations of the fermionic fields are defined as follows:

$$
\begin{array} { l } { { \delta \psi ( x ) = i \left[ \alpha _ { V } ^ { f } ( x ) \frac { \lambda ^ { f } } { 2 } + \alpha _ { A } ^ { f } ( x ) \frac { \lambda ^ { f } } { 2 } \gamma _ { 5 } \right] \psi ( x ) } } \\ { { \delta \bar { \psi } ( x ) = - i \bar { \psi } ( x ) \left[ \alpha _ { V } ^ { f } ( x ) \frac { \lambda ^ { f } } { 2 } - \alpha _ { A } ^ { f } ( x ) \frac { \lambda ^ { f } } { 2 } \gamma _ { 5 } \right] ~ . } } \end{array}
$$

With degenerate bare quark masses $m _ { 0 }$ , the global vector transformations are a symmetry of the action. From the corresponding local transformations (eqs. (15) with $\alpha _ { A } ^ { f } = 0$ ) vector WIs can be derived, which indicate that the conserved lattice vector current is a point-split operator. The no-renormalization theorem applies to this operator. On the other hand, the lattice vector local current of eq. (6) is subject of a finite renormalization $Z _ { V } ( g _ { 0 } ^ { 2 } ) \ne 1$ , due to its non-conservation on the lattice. A detailed treatment of these results can be found in refs. [15,19]. Here we illustrate the point by writing down a typical lattice vector WI. For simplicity, the flavour $f$ of the vector variation of eq. (15) is chosen to be non-singlet and non-diagonal (e.g. the lowering operator $\lambda ^ { f } = \lambda ^ { 1 } - i \lambda ^ { 2 }$ ). The two quark masses related to this flavour are denoted by $m _ { 1 }$ and $m _ { 2 }$ (i.e. they are non-degenerate). Then we obtain the WI

$$
\begin{array} { r } { \displaystyle \frac { \delta } { \delta \alpha _ { V } ^ { f } ( x ) } \langle V _ { \nu } ^ { g } ( 0 ) \rangle = 0 \Longleftrightarrow } \\ { \displaystyle Z _ { V } \sum _ { \mu } \nabla _ { x } ^ { \mu } \langle V _ { \mu } ^ { f } ( x ) V _ { \nu } ^ { g } ( 0 ) \rangle = ( m _ { 2 } - m _ { 1 } ) \langle S ^ { f } ( x ) V _ { \nu } ^ { g } ( 0 ) \rangle \ , } \end{array}
$$

where $a \nabla _ { x } ^ { \mu } f ( x ) = f ( x ) - f ( x - \mu )$ is an asymmetric lattice derivative. We have imposed $x \neq 0$ , in order to avoid contact terms. We can now multiply both sides of the WI by $Z _ { V }$ and require that the renormalized quantities obey the nominal vector WI. This implies that the product of the quark mass and the scalar density is renormalization group invariant. In other words, the renormalization of the scalar density $S ( x )$ is the inverse of the multiplicative renormalization of the quark mass:

$$
Z _ { S } = Z _ { m } ^ { - 1 } .
$$

Far less straightforward is the implementation of the axial symmetry with Wilson fermions, because of the presence of the chiral symmetry breaking Wilson term in the action. The basic idea is that, by imposing suitable renormalization conditions, PCAC is recovered in the continuum. Consider the axial

WI, obtained from eq. (15) with $\alpha _ { V } ^ { f } = 0$ , arising from the variation of the non-singlet pseudoscalar density. For $f , g \neq 0$ we have:

$$
\begin{array} { r l r } & { } & { \frac { \delta } { \delta \alpha _ { A } ^ { f } ( x ) } \langle P ^ { g } ( 0 ) \rangle = 0 \Longleftrightarrow } \\ & { } & { - \delta ( x ) \delta ^ { f g } \frac { 1 } { N _ { f } } \langle S ^ { 0 } ( 0 ) \rangle = \displaystyle \sum _ { \mu } \nabla _ { x } ^ { \mu } \langle A _ { \mu } ^ { f } ( x ) P ^ { g } ( 0 ) \rangle - 2 m _ { 0 } \langle P ^ { f } ( x ) P ^ { g } ( 0 ) \rangle } \\ & { } & { - \langle X ^ { f } ( x ) P ^ { g } ( 0 ) \rangle \ . } \end{array}
$$

The l.h.s. of this equation is the variation of the operator which vanishes when $x \neq 0$ . The r.h.s. is the variation of the action from which the lattice version of the standard PCAC relation is recovered. The operator $X ^ { f }$ in the above WI arises from the variation of the Wilson term under the axial transformation 6. Unlike the vector current case, $X ^ { f }$ cannot be cast in the form of a four-divergence. It is a dimension-4 operator which, in the naive continuum limit vanishes, being of the form ( $a \times$ dimension-5 operator). However, it has divergent matrix elements beyond tree-level. Its mixing with operators of equal and lower dimensions can be expressed as follows [15]:

$$
\overline { { { X } } } ^ { f } ( x ) = X ^ { f } ( x ) + 2 \overline { { { m } } } P ^ { f } ( x ) + ( Z _ { A } - 1 ) \nabla _ { x } ^ { \mu } A _ { \mu } ^ { f } ( x ) \ ,
$$

where naive dimensional arguments tell us that the mixing constant $Z _ { A } ( g _ { 0 } ^ { 2 } )$ is finite, whereas $\overline { { m } } ( g _ { 0 } ^ { 2 } , m _ { 0 } )$ diverges linearly like $a ^ { - 1 }$ . Logarithmic divergences have been shown to be absent at all orders in PT [32]. The two mixing constants $Z _ { A }$ and $\overline { { m } }$ , and consequently $\overline { { X } } ^ { f } ( x )$ , are defined by requiring that the renormalized axial current satisfies the nominal (continuum) axial WI. This requirement is satisfied upon demanding that, in the continuum limit and for vanishing quark mass, the correlation functions of ${ \overline { { X } } } ^ { g }$ vanishes except for localized contact terms. For the insertion of ${ \overline { { X } } } ^ { f }$ in correlation functions such as $\langle \overline { { X } } ^ { f } ( x ) P ^ { g } ( 0 ) \rangle$ the contact terms, by flavour symmetry, must have the form

$$
\langle \overline { { { X } } } ^ { f } ( x ) P ^ { g } ( 0 ) \rangle = \frac { 1 } { N _ { f } } b _ { 0 } \delta ( x ) \delta ^ { f g } + b _ { 1 } \delta ^ { f g } \sqsupset \delta ( x ) ~ ,
$$

where the last term on the r.h.s. is a Schwinger term (it vanishes upon integration over $\textstyle \int d ^ { 4 } x )$ . Combining eqs. (18), (19) and (20) we arrive at

$$
Z _ { A } \nabla _ { x } ^ { \mu } \langle A _ { \mu } ^ { f } ( x ) P ^ { g } ( 0 ) \rangle - 2 ( m _ { 0 } - \overline { { { m } } } ) \langle P ^ { f } ( x ) P ^ { g } ( 0 ) \rangle - b _ { 1 } \delta ^ { f g } \sqsubset \partial \delta ( x )
$$

$$
= - \delta ( x ) \delta ^ { f g } \frac { 1 } { N _ { f } } \left[ \langle S ^ { 0 } \rangle - b _ { 0 } \right] .
$$

By multiplying both sides of the above by $Z _ { P }$ , we obtain a completely finite expression, which is the renormalized WI in the continuum. The $b _ { 1 }$ -subtraction on the l.h.s. compensates the contact terms arising in the correlation functions $\nabla _ { x } ^ { \mu } \langle \hat { A } _ { \mu } ^ { f } ( x ) \hat { P } ^ { g } ( 0 ) \rangle$ and $\langle \hat { P } ^ { f } ( x ) \hat { P } ^ { g } ( 0 ) \rangle$ , when $x$ is close to the origin. Otherwise, we have recovered the nominal WI, which is now satisfied by renormalized operators. This implies the following useful results:

(1) The product $\hat { A } _ { \mu } ^ { f } = Z _ { A } A _ { \mu } ^ { f }$ must be interpreted as the renormalized axial current ( $Z _ { A } ( g _ { 0 } ^ { 2 } ) \neq 1$ is a finite renormalization).

(2) The chiral limit is to be defined as the value $m _ { C }$ of $m _ { 0 }$ for which the difference $m _ { C } - \overline { { m } } ( g _ { 0 } ^ { 2 } , m _ { C } )$ vanishes. We see from eq. (21) that the product $( m _ { 0 } - { \overline { { m } } } ) P ^ { f }$ is renormalization group invariant. Thus, we can write the renormalized quark mass as:

$$
m _ { R } = \overline { { { Z } } } _ { m } [ m _ { 0 } - \overline { { { m } } } ( m _ { 0 } ) ] ,
$$

and obtain for the renormalization of the pseudoscalar density

$$
Z _ { P } = \overline { { Z } } _ { m } ^ { - 1 } \ .
$$

Combining eqs. (4), (17), (22) and (23), we obtain for the ratio $Z _ { P } / Z _ { S }$

$$
\frac { Z _ { P } } { Z _ { S } } = \frac { m _ { 0 } - \overline { { m } } ( m _ { 0 } ) } { m _ { 0 } - m _ { C } } .
$$

(3) Two useful expressions are derived from eq. (21), when $x \neq 0$ (i.e. all terms proportional to $\delta ( x )$ vanish). The first expression consists simply in the definition of the ratio (for $f = g$ ):

$$
2 \rho = \frac { 2 ( m _ { 0 } - \overline { { { m } } } ) } { Z _ { A } } = \frac { \int d ^ { 3 } x \nabla _ { x } ^ { 0 } \langle A _ { 0 } ^ { f } ( x ) P ^ { f } ( 0 ) \rangle } { \int d ^ { 3 } x \langle P ^ { f } ( x ) P ^ { f } ( 0 ) \rangle } ,
$$

where the spatial derivatives of the divergence of the axial current vanish under the integration. This quantity, once $Z _ { A }$ is known, can be used for the determination of $( m _ { 0 } - { \overline { { m } } } )$ . The second useful expression arises upon applying the LSZ reduction formula:

$$
Z _ { A } \nabla _ { 0 } ^ { \mu } \langle 0 | A _ { \mu } ^ { f } ( 0 ) | P \rangle = 2 ( m _ { 0 } - { \overline { { { m } } } } ) \langle 0 | P ^ { f } ( 0 ) | P \rangle ,
$$

where $| P \rangle$ is the lowest lying pseudoscalar state (in the chiral limit this is the Goldostone boson). Both equations are lattice versions of the standard PCAC

relation.

(4) Last but not least, the $b _ { 0 }$ -subtraction on the r.h.s. of eq. (21) implies that the proper definition of the lattice chiral condensate with Wilson fermions, determined by axial WIs, is given by

$$
\langle \bar { \psi } \psi \rangle _ { s u b } = \left[ \langle { \cal S } ^ { 0 } \rangle - b _ { 0 } \right] .
$$

The above "subtracted chiral condensate" is a logarithmically divergent quantity, which is multiplicatively renormalized by $Z _ { P }$ :

$$
\begin{array} { r } { \langle \bar { \psi } \psi \rangle = Z _ { P } \left[ \langle S ^ { 0 } \rangle - b _ { 0 } \right] . } \end{array}
$$

To obtain an expression from which the chiral condensate can be computed in a convenient way, we need to integrate both sides of WI (21) over all spacetime. Upon integration the Schwinger $b _ { 1 }$ -term on the l.h.s. vanishes. For the remaining terms on the l.h.s. two possibilities arise: either work directly in the chiral limit or work with a non-zero quark mass and subsequently extrapolate to the chiral limit. If we work directly in the chiral limit, the $( m _ { 0 } - { \overline { { m } } } )$ term vanishes. The integrated $\nabla _ { x } ^ { \mu } \hat { A } _ { \mu } ^ { f }$ would also vanish, being the integral of a four-divergence, if chiral symmetry were not broken (absence of Goldstone bosons). In this case the WI implies a vanishing chiral condensate. If the symmetry is broken, the presence of a Goldstone boson guarantees a nonvanishing surface term upon integrating $\nabla _ { x } ^ { \mu } \hat { A } _ { \mu } ^ { f }$ . Thus, also the r.h.s. of the WI (i.e. the chiral condensate) is non-zero. In the present work, in accordance with standard numerical computations, we will be working with a finite quark mass subsequently extrapolated to the chiral limit. Then the surface term of the divergence of the axial current vanishes, and the integrated WI becomes (for $f = g$ )

$$
{ \frac { 1 } { N _ { f } } } \langle \bar { \psi } \psi \rangle _ { s u b } = \operatorname * { l i m } _ { m _ { 0 }  m _ { C } } 2 ( m _ { 0 } - \overline { { { m } } } ) \int d ^ { 4 } x \langle P ^ { f } ( x ) P ^ { f } ( 0 ) \rangle \ .
$$

From the above expression, we can also derive two other useful identities for the chiral condensate. They are equivalent, in theory, to the one above, in the chiral limit and up to discretization errors. Since in practice they are implemented at finite quark mass and lattice spacing, they may reveal the importance of these systematic effects. Both expressions are obtained by writing the correlation function $\left. P ^ { f } ( x ) P ^ { f } ( 0 ) \right.$ as a time-ordered product and by inserting a complete set of states in standard fashion. The spatial integration $\int d ^ { 3 } x$ projects zeromomentum states. Contributions from higher mass states vanish in the chiral limit. Upon performing the time integration $\textit { j } d t$ , we find

$$
{ \frac { 1 } { N _ { f } } } \langle \bar { \psi } \psi \rangle _ { s u b } = - \operatorname * { l i m } _ { m _ { 0 }  m _ { C } } { \frac { ( m _ { 0 } - { \overline { { m } } } ) } { m _ { P } ^ { 2 } } } { | \langle 0 | P ( 0 ) | P \rangle | } ^ { 2 } ,
$$

where $m p$ is the mass of the pseudoscalar state $| P \rangle$ . Note that the pseudoscalar operator without a flavour index is defined above as, say, $P = \bar { u } \gamma _ { 5 } d$ (and the pion state has a $d u$ flavour content). The same is true for the axial current (without flavour index) $A _ { \mu } = \bar { u } \gamma _ { \mu } \gamma _ { 5 } d$ .

Another expression for the chiral condensate is obtained by substituting in the above equation the matrix element $\langle 0 | P ( 0 ) | P \rangle$ by $\nabla _ { 0 } ^ { \mu } \langle 0 | \hat { A } _ { \mu } ( 0 ) | P \rangle$ (see eq.(26)). The matrix element of the axial current is parameterized in terms of the pion decay constant:

$$
\langle 0 | \hat { A } _ { \mu } ( x ) | P ( \vec { p } ) \rangle = i f _ { P } p _ { \mu } \exp ( - i p x ) ,
$$

( $f _ { P }$ is the decay constant of the pseudoscalar state $| P \rangle$ ) and thus

$$
\nabla _ { 0 } ^ { \mu } \langle 0 | \hat { A } _ { \mu } ( 0 ) | P ( \vec { 0 } ) \rangle = f _ { P } m _ { P } ^ { 2 } \ .
$$

Combining eqs. (26), (30) and (32) we find

$$
\frac { 1 } { N _ { f } } \langle \bar { \psi } \psi \rangle _ { s u b } = - \operatorname * { l i m } _ { m _ { 0 }  m _ { C } } \frac { f _ { P } ^ { 2 } m _ { P } ^ { 2 } } { 4 ( m _ { 0 } - \overline { { { m } } } ) } ,
$$

which is the familiar Gell-Mann-Oakes-Renner (GMOR) relation [16]. Note that the non-vanishing of the chiral condensate in the above equation implies the well known linear dependence of the pseudoscalar mass squared from the quark mass.

All determinations of the chiral condensate discussed above are subject to taking the chiral limit ( $m  m _ { C }$ ) and the continuum limit ( $a \ \to \ 0$ ). In practice we work at non-zero quark mass $m _ { 0 }$ (fixed hopping parameter $\kappa$ )and lattice spacing $a$ (fixed inverse square coupling $\beta$ ). At fixed $\beta$ the chiral limit is normally taken by extrapolating the lattice data to $\kappa _ { C }$ (see comments in the introduction on the reliability of this extrapolation and the presence of chiral logarithms). Assuming that this extrapolation is reliable, we obtain a result which is subject to systematic errors primarily due to the finiteness of the lattice spacing. These effects are $\mathcal O ( a )$ for Wilson fermions. To reduce these errors, improved actions in the spirit of Symanzik's programme [26] must be used. We have used the tree-level improved Clover action and operators of refs. [10,11], which leaves us with $\mathcal { O } ( a g _ { 0 } ^ { 2 } )$ discretization errors. Comparison of results obtained with the Wilson and Clover actions enable us to estimate finite lattice spacing effects. Some technical details on the implementation of chiral WIs with the tree-level Clover action are given in Appendix C.

# 5 Lattice measurements of the chiral condensate

Having shown how lattice WIs uniquely determine the subtracted chiral condensate with Wilson fermions, in this section we will explore several ways of measuring it from simulations. We will first discuss some lattice details in subsect. 5.1. In subsect. 5.2 we will explore several standard methods, previously used in the literature, for extracting the chiral condensate from lattice simulations. They are based on eqs. (29), (30) and (33). Their main shortcoming is that, upon passing from lattice to physical units, the systematic error of the condensate due to the uncertainty in the lattice calibration is amplified. Another (minor) error amplification is due to the extrapolation of the result from the strange quark region to the chiral limit. In subsect. 5.3 we will present what we consider the best method of measuring the lattice condensate, which sidesteps these problems. The reader uninterested in lattice technicalities may skip these subsections. Finally, in subsect. 5.4 we discuss our final results.

# 5.1 Computational details

In this work we have used the results of the best runs performed by the APE group in the last few years. Both the Wilson (W) and the tree level Clover (C) actions have been used at three values of the coupling $\beta = 6 / g _ { 0 } ^ { 2 }$ and several quark masses in the strange quark region (the range of the lowest-lying pseudoscalar meson in these runs is in the range $m _ { K } - 2 m _ { K }$ ). Table 1 shows the parameters of the runs from which we have extracted the necessary masses and matrix elements. A discussion of the technical details can be found in [8] and [9]. Here we just state that all statistical errors have been estimated with the jacknife method by blocking the data in groups of 10 configurations. We have checked that our error estimates are roughly independent of the blocking size. The APE collaboration has also performed dedicated runs at the same $\beta$ values with both the Wilson and the tree-level Clover actions, in order to obtain NP estimates of the renormalization constants in the RI scheme [28]. The parameters of these runs and the results for the renormalization constants of interest to the present work are reported in table 2. More results and details on the error quoted can be found in ref. [28]. The main qualitative conclusion of ref. [28] (as seen in Table 2) is that the NP method gives excellent estimates of the renormalization constants $Z _ { A }$ and $Z _ { S }$ , whereas $Z _ { P }$ (and subsequently the ratio $Z _ { P } / Z _ { S }$ ) are subject to a bigger systematic error.

In lattice simulations we use dimensionless quantities; we thus define the operators (in lattice units) 7

<table><tr><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>C60a</td><td rowspan=1 colspan=1>C60b</td><td rowspan=1 colspan=1>C62</td><td rowspan=1 colspan=1>C64</td><td rowspan=1 colspan=1>W60</td><td rowspan=1 colspan=1>W62</td><td rowspan=1 colspan=1>W64</td></tr><tr><td rowspan=1 colspan=1>β# ConfsVolume</td><td rowspan=1 colspan=1>6.0490183 × 64</td><td rowspan=1 colspan=1>6.06004|243 × 40</td><td rowspan=1 colspan=1>6.2250243 × 64|243 × 64|</td><td rowspan=1 colspan=1>6.4400243 × 64|243 × 64|</td><td rowspan=1 colspan=1>6.0320183 × 64 </td><td rowspan=1 colspan=1>6.2250243 × 642</td><td rowspan=1 colspan=1>6.4400243 × 64</td></tr><tr><td rowspan=1 colspan=1>κ</td><td rowspan=1 colspan=1>0.14250.14320.1440</td><td rowspan=1 colspan=1>0.14250.14320.1440</td><td rowspan=1 colspan=1>0.141440.141840.142240.14264</td><td rowspan=1 colspan=1>0.14000.14030.14060.1409</td><td rowspan=1 colspan=1>0.15300.15400.1550</td><td rowspan=1 colspan=1>0.15100.15150.15200.1526</td><td rowspan=1 colspan=1>0.14880.14920.14960.1500</td></tr><tr><td rowspan=1 colspan=1>t1 − t2</td><td rowspan=1 colspan=1>15-28</td><td rowspan=1 colspan=1>15-19</td><td rowspan=1 colspan=1>18-28</td><td rowspan=1 colspan=1>24-30</td><td rowspan=1 colspan=1>15-28</td><td rowspan=1 colspan=1>18-28</td><td rowspan=1 colspan=1>24-30</td></tr><tr><td rowspan=1 colspan=1>[a−1(K*)|</td><td rowspan=1 colspan=1>2.12(6)</td><td rowspan=1 colspan=1>2.16(4)</td><td rowspan=1 colspan=1>2.7(1)</td><td rowspan=1 colspan=1>4.0(2)</td><td rowspan=1 colspan=1>2.26(5)</td><td rowspan=1 colspan=1>3.00(9)</td><td rowspan=1 colspan=1>4.1(2)</td></tr></table>

Table 1   
Summary of the parameters of the runs with the Clover (C) and Wilson (W) fermion actions used in this work. We also show the time intervals in which correlation functions were fitted in order to extract meson masses and decay constants. The bottom line is our preferred value for the lattice spacing, obtained by fixing the mass of the vector $K ^ { * }$ -meson. The error is statistical.

$$
\begin{array} { l } { { \displaystyle { \mathcal P } ( x ) = a ^ { 3 } P ( x ) } } \\ { { \displaystyle \mathcal A _ { \mu } ( x ) = a ^ { 3 } A _ { \mu } ( x ) } } \\ { { \displaystyle \quad \langle \bar { \chi } \chi \rangle = a ^ { 3 } \frac 1 { N _ { f } } \langle \bar { \psi } \psi \rangle _ { s u b } ~ . } } \end{array}
$$

The mass and decay constant of the pseudoscalar meson in lattice units are denoted by $M _ { P } = a m _ { P }$ and $F _ { P } = a f _ { P }$ The extraction of quantities of interest in lattice units at fixed quark mass (e.g. $M _ { P }$ , $F _ { P } ^ { \prime }$ , $a \rho$ ) consists in fitting the large time behaviour of correlation functions of suitable operators. Since this is a standard procedure we will not present it here; the interested reader can find all relevant details of our analysis in ref. [8]. We recall that in terms of the Wilson hopping parameter $\kappa$ , the bare quark mass is given by8

$$
a ( m _ { 0 } - m _ { C } ) = \frac { 1 } { 2 \kappa } - \frac { 1 } { 2 \kappa _ { C } } + \mathcal { O } ( a ) .
$$

To pass from lattice units to physical ones, the experimental value of several

<table><tr><td rowspan=1 colspan=1>β</td><td rowspan=1 colspan=1>6.0</td><td rowspan=1 colspan=1>6.0</td><td rowspan=1 colspan=1>6.2</td><td rowspan=1 colspan=1>6.2</td><td rowspan=1 colspan=2>6.4</td><td rowspan=1 colspan=1>6.4</td></tr><tr><td rowspan=2 colspan=1>Action# ConfsVolume</td><td rowspan=2 colspan=1>C100163 × 32</td><td rowspan=2 colspan=1>W100163 × 32</td><td rowspan=2 colspan=1>C180163 × 32</td><td rowspan=1 colspan=1>W100</td><td rowspan=2 colspan=2>C602243 × 32[</td><td rowspan=2 colspan=1>W60243 × 32</td></tr><tr><td rowspan=1 colspan=1>163 × 32</td></tr><tr><td rowspan=2 colspan=1>κ</td><td rowspan=2 colspan=1>0.14250.14320.1440</td><td rowspan=2 colspan=1>0.15300.15400.1550</td><td rowspan=2 colspan=1>0.141440.141840.142240.14264</td><td rowspan=1 colspan=1>0.15100.1515</td><td rowspan=2 colspan=2>0.14000.14030.14060.1409</td><td rowspan=2 colspan=1>0.14880.14920.14960.1500</td></tr><tr><td rowspan=1 colspan=1>0.15200.1526</td></tr><tr><td rowspan=4 colspan=1>ZSZPZP /ZsZA</td><td rowspan=4 colspan=1>0.83(2)0.41(6)0.49(6)1.05(3)</td><td rowspan=4 colspan=1>0.68(1)0.45(6)0.66(8)0.81(1)</td><td rowspan=3 colspan=1>0.85(2)0.47(5)0.56(5)</td><td rowspan=1 colspan=1>0.72(1)0.50(5)</td><td rowspan=1 colspan=2>0.85(2)0.55(3)</td><td rowspan=4 colspan=1>0.74(1)0.57(4)0.77(5)0.82(1)</td></tr><tr><td rowspan=2 colspan=1>0.69(7)</td><td rowspan=2 colspan=2>0.65(2)</td><td rowspan=1 colspan=1>−</td></tr><tr><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>1.02(2)</td><td rowspan=1 colspan=1>0.81(1)</td><td rowspan=1 colspan=2>1.01(1)</td></tr></table>

Table 2 Parameters of the runs used for the NP calculation of the renormalization constants and values of $Z _ { S }$ , $Z _ { P }$ , $Z _ { P } / Z _ { S }$ and $Z _ { A }$ at a scale $\mu a \simeq 1$ . All results are in the chiral limit. The error includes both statistical and systematic effects.

physical quantities is required (one for setting the scale and one per quark flavour for setting the corresponding quark mass). In practice more than one physical quantity is used for the determination of the lattice spacing, in order to estimate the systematic error of the lattice calibration. In the present work we have set the scale $a ^ { - 1 }$ from the vector meson mass $M _ { K ^ { * } } = 0 . 8 9 3 1$ GeV. As suggested in ref. [8], this gives the most reliable estimate, since it corresponds to the strange quark mass region, where our simulations have been carried out, and it is free of any uncertainties due to renormalization. The result (obtained with the methods of ref. [8] on the dataset of ref. [9]) is shown in Table 1. We see that the statistical error is about 2%-5%. The systematic error on the lattice spacing has been estimated as the spread of the $a ^ { - 1 }$ values obtained from 10 other physical quantities 9 , as listed in Table 11 of ref. [8]. This error is about 7%. Thus, we quote an overall 8% error due to the uncertainties in the lattice calibration. The discrepancies between the various lattice calibrations is believed to reflect, at least in part, the effect of quenching. We will also be using some physical quantities in order to determine the chiral limit, quark masses etc. In the region of the up and down quark masses we have used the experimental values $m _ { \pi } = 0 . 1 3 8$ GeV and $f _ { \pi } = 0 . 1 3 0 7$ GeV. In the region of the strange quark mass we have taken $m _ { K } = 0 . 4 9 5$ GeV and $f _ { K } = 0 . 1 5 9 8$ GeV. In the chiral limit ( $m _ { P } = 0$ ) we use the "experimental"

value $f _ { \chi } = 0 . 1 2 8 2$ GeV obtained by considering $f _ { K }$ and $f _ { \pi }$ , as functions of $m _ { K } ^ { 2 }$ and $m _ { \pi } ^ { 2 }$ respectively(PCAC) and by extrapolating them naively to vanishing $m p$ .

5.2 Standard lattice measurements of the dimensionless chiral condensate

We will now present several "standard" measurements of the subtracted chiral condensate based on eqs. (29), (30) and (33). In theory all three equations are equivalent in the continuum limit. However, since they are implemented in practice at fixed lattice spacing (fixed $\beta$ ) and non-zero quark mass, they are subject to systematic errors. Thus, by obtaining the chiral condensate in several ways, we are able to study the systematic errors introduced by the finiteness of the lattice spacing, the extrapolations to the chiral limit etc.

Once the hadronic correlation functions have been fitted and the lattice masses and matrix elements have been extracted, we must extrapolate the results to the chiral limit. The standard way consists in determining the chiral point $\kappa _ { C }$ by fitting the pseudoscalar mass $M _ { P } ^ { 2 }$ as a linear function of the quark masses

$$
M _ { P } ^ { 2 } = C ^ { H S } \left( \frac { 1 } { \kappa } - \frac { 1 } { \kappa _ { C } } \right) ~ .
$$

The slope $C ^ { H S }$ related to the hadron spectrum $( H S )$ is also determined from this fit. The problem is that this relation is based on eq. (35) which, for our Clover dataset is only satisfied to $\mathcal O ( a )$ and therefore it should not be implemented for this action (see Appendix C). In order to avoid this problem, we fit instead $a \rho$ as a linear function of $M _ { P } ^ { 2 }$

$$
2 a \rho = { \frac { 1 } { C ^ { A W I } } } M _ { P } ^ { 2 } \ ,
$$

and determine the slope $C ^ { A W I }$ related to the Axial WIs $( A W I )$ . We have verified for the Wilson data that $a \rho$ vanishes at $\kappa _ { C }$ within our statistical errors.

The most straightforward way of computing the subtracted condensate is that of eq. (29), where the term $( m _ { 0 } - { \overline { { m } } } )$ has been computed from $a \rho$ , with the aid of eq. (25):

$$
\frac { \langle \bar { \chi } \chi \rangle } { Z _ { A } } = \operatorname * { l i m } _ { M _ { P } ^ { 2 } \to 0 } a \rho \sum _ { x } \langle { \mathcal P } ( x ) { \mathcal P } ^ { \dagger } ( 0 ) \rangle \ .
$$

Unfortunately this method is unreliable: at $\beta ~ = ~ 6 . 0$ and 6.2 we obtained results which are clearly incompatible with those obtained from other methods; at $\beta = 6 . 4$ the statistical error overwhelms the signal. The reason is that the correlation function $\left. \mathcal { P } ( x ) \mathcal { P } ^ { \dagger } ( 0 ) \right.$ , integrated over all space-time, has a contact term at $x \approx 0$ which from dimensional counting is seen to behave as $1 / a ^ { 2 }$ . Although its contribution to the chiral condensate vanishes in the chiral limit (it is to be multiplied by $( m _ { 0 } - { \overline { { m } } } ) )$ , it is the dominant contribution at the values of lattice spacing and quark mass of the simulation. Therefore the extrapolation to $m _ { C }$ becomes unreliable and with big statistical errors. Thus, this computation of the condensate must be discarded. This effect was first observed (with much poorer statistics and control of the systematics) in ref. [12]. However we note, that in ref. [13] a variant of the above determination, performed with smeared operators, gave satisfactory results in the chiral limit. Presumably, smearing reduces the influence of the contact term on the data.

Two robust ways of computing the subtracted chiral condensate, which avoid contact terms, are based on eqs. (30) and (33). Once more eq. (25) is used to compute the term $( m _ { 0 } - { \overline { { m } } } )$ . We find for $\langle \bar { \chi } \chi \rangle$ :

$$
\begin{array} { r l } & { \frac { \left. \bar { \chi } \chi \right. _ { 1 } } { Z _ { A } } { = - \frac 1 { 2 C ^ { A W I } } \displaystyle \operatorname* { l i m } _ { M _ { P } ^ { 2 } \to 0 } \left| \left. 0 | \mathcal P ( 0 ) | P \right. \right| ^ { 2 } } } \\ & { \frac { \left. \bar { \chi } \chi \right. _ { 2 } } { Z _ { A } } { = - \frac 1 2 C ^ { A W I } \displaystyle \operatorname* { l i m } _ { M _ { P } ^ { 2 } \to 0 } F _ { P } ^ { 2 } . } } \end{array}
$$

Other computations of the subtracted chiral condensate are also possible by using eqs. (24),(35) and (36) in order to express $( m _ { 0 } - \overline { { m } } )$ in terms of $\left( m _ { 0 } - m _ { C } \right)$ . We avoid this option for two reasons. First, it requires the ratio $Z _ { P } / Z _ { S }$ which, compared to $Z _ { A }$ , has a larger error (see results in Table 2). Second, it involves eq. (35) which cannot be implemented with our Clover data.

<table><tr><td rowspan=1 colspan=1>Run</td><td rowspan=1 colspan=1>¯xχ1/ZA</td><td rowspan=1 colspan=1>hχχ2/ZA</td><td rowspan=1 colspan=1>αa−3h¯χχ1</td><td rowspan=1 colspan=1>a−3h¯χ2</td></tr><tr><td rowspan=2 colspan=1>C60aC60bW60</td><td rowspan=2 colspan=1>0.0040(3)0.0039(2)0.0045(3)</td><td rowspan=2 colspan=1>0.0040(3)0.0039(2)0.0051(4)</td><td rowspan=1 colspan=1>0.044(4)0.045(3)</td><td rowspan=2 colspan=1>0.044(4)0.046(3)0.053(6)</td></tr><tr><td rowspan=1 colspan=1>0.048(4)</td></tr><tr><td rowspan=1 colspan=1>C62W62</td><td rowspan=1 colspan=1>0.0018(2)0.0019(1)</td><td rowspan=1 colspan=1>0.0015(2)0.0016(1)</td><td rowspan=1 colspan=1>0.041(7)0.046(5)</td><td rowspan=1 colspan=1>0.034(5)0.038(4)</td></tr><tr><td rowspan=1 colspan=1>C64W64</td><td rowspan=1 colspan=1>0.0007(1)0.0008(1)</td><td rowspan=1 colspan=1>0.0007(1)0.0007(1)</td><td rowspan=1 colspan=1>0.047(12)0.053(11)</td><td rowspan=1 colspan=1>0.049(8)0.047(7)</td></tr></table>

Table 3 The subtracted chiral condensate determined from eqs. (39) in lattice units (second and third columns) and in physical, GeV $^ { 3 }$ , units (fourth and fifth columns). An overall minus sign has been omitted from the data. All errors are statistical.

In Table 3 we present our results. Our comments are the following:

•The results from runs C60a and C60b, which only differ by a volume factor of about 1.5, are perfectly compatible. We see no finite volume effects.   
All $\langle \bar { \chi } \chi \rangle _ { 1 } / Z _ { A }$ and $\langle \bar { \chi } \chi \rangle _ { 2 } / Z _ { A }$ results, obtained with the Wilson action at the same $\beta$ , are compatible within two standard deviations. Those from the Clover action are compatible within one standard deviation. This suggests that our chiral extrapolations are robust, especially in the Clover case where the discretization error is reduced.   
The $\beta = 6 . 4$ results have a larger relative error for both actions. This is probably due to the small physical volume at this value of $\beta$ . We recall that the quark mass measured at this coupling on the same dataset in ref. [9] also suffered from finite lattice effects. Thus the $\beta = 6 . 4$ results will not be taken into account.

Earlier works also computed the subtracted condensate at $\beta = 6 . 0$ with the Wilson action. The methods used were all based on eq.(33) (GMOR), so they are comparable to our $\langle \bar { \chi } \chi \rangle _ { 2 } / Z _ { A }$ value:

$\langle \bar { \chi } \chi \rangle _ { 2 } / Z _ { A } = - 0 . 0 0 5 1 ( 4 )$ This work $\langle \bar { \chi } \chi \rangle / Z _ { A } = - 0 . 0 0 6 ( 1 )$ Ref. [12] $\langle \bar { \chi } \chi \rangle / Z _ { A } = - 0 . 0 0 6 ( 3 )$ Ref. [13] .

We see that our result is compatible with previous determinations, but the statistical accuracy is greatly improved.

The standard method for computing the chiral condensate would now involve renormalizing $\langle \bar { \chi } \chi \rangle _ { 1 }$ or $\langle \bar { \chi } \chi \rangle _ { 2 }$ by multiplying it with $Z _ { P }$ (c.f. eqs.(27) and (28)) and expressing it in physical units by multiplying it with $a ^ { - 3 }$ (c.f. eq.(34)). The latter step is the source of error amplification, due to the cubic power of the lattice spacing. Our 8% overall estimate of the error due to the lattice calibration would be amplified to $2 4 \%$ . This increased error is a serious shortcoming of all standard methods of computing the condensate on the lattice.

# 5.3 Lattice measurements of the physical chiral condensate

We will now present a better method of computing the chiral condensate on the lattice. Rather than the standard two-step approach (i.e. first measure the dimensionless condensate in terms of dimensionless quantities and then multiply it by $a ^ { - 3 }$ ) we will obtain it in physical units, by using the physical value of the pseudoscalar decay constant multiplied by $a ^ { - 1 }$ . Thus we write the

GMOR relation (33) for the renormalized condensate as follows

$$
{ \frac { 1 } { N _ { f } } } \langle \bar { \psi } \psi \rangle _ { 1 } = - { \frac { 1 } { 2 } } a ^ { - 1 } f _ { \chi } ^ { 2 } Z _ { S } C ^ { H S } \ .
$$

Use of eqs.(24) and (36) has been made. Note that $f _ { \chi } = 0 . 1 2 8 2$ GeV is the "experimental" value in physical units (see subsect. 5.1). The computation of the condensate based on the above formula has several advantages. The most important is that, by expressing it in terms of $f _ { \chi }$ , we are left with only one power of the UV cutoff $a ^ { - 1 }$ . Thus, we avoid the error amplification of the standard methods. Another advantage is that there is no extrapolation to the chiral limit. Instead, we only need to determine the slope $C ^ { H S }$ of the squared mass of the pseudoscalar meson with respect to the quark mass. This is reminiscent of similar procedures used in spectroscopic studies in refs. [8,33]. Note that since we work at $\kappa$ values typical of the strange quark mass, we are implicitly assuming that the slope will not change in the chiral region. There is no way round this assumption, besides simulating at lighter quark masses. Finally, another advantage of eq. (41) is that the renormalization constant is now $Z _ { S }$ , which is much better determined than $Z _ { P }$ with the NP method.

The only drawback of the above determination is that it cannot be used on our Clover data as explained in Appendix C. The alternative is to write eq. (33), in terms of eq.(37), as follows

$$
{ \frac { 1 } { N _ { f } } } \langle \bar { \psi } \psi \rangle _ { 2 } = - { \frac { 1 } { 2 } } a ^ { - 1 } f _ { \chi } ^ { 2 } { \frac { Z _ { P } } { Z _ { A } } } C ^ { A W I } \ .
$$

Also this determination has the advantages of having only a single power of the inverse lattice spacing and no explicit chiral extrapolation; a disadvantage is the larger systematic error of $Z _ { P }$ .

Since we use the NP estimates of the renormalization constants in both determinations of eqs. (41) and (42), our results are obtained in the RI scheme at a renormalization scale $a ^ { - 1 }$ . Subsequently, perturbation theory is used to NNLO to express them in the $\overline { { \mathrm { M S } } }$ scheme (c.f. eq. (12)) and run them to the conventional reference scale of 2 GeV (c.f. eq. (13)). In so doing, two further uncertainties are introduced. The first is the $\mathcal { O } ( \alpha _ { s } ^ { 3 } )$ error due to the truncation of the NNLO perturbative series which we consider negligible. The second uncertainty arises because the initial scale $\mu \simeq a ^ { - 1 }$ is known with an 8% precision. Since the dependence of the perturbative coefficients on the scale is logarithmic, we find that the error introduced is at most 2% and will be ignored.

In conclusion there are three important errors in our results. Firstly, we have a purely statistical eror of the lattice quantiies $C ^ { H S }$ and $C ^ { A W I }$ .Secondly, there is the 8% error of statistical and systematic origin from the lattice calibration (the $a ^ { - 1 }$ factor in eqs.(41) and (42)). Finally, we have the error in the determination of the renormalization constants (c.f. Table 2) which is also of both statistical and systematic origin.

In table 4 we present our results for the chiral condensate in the $\mathrm { \overline { { M S } } }$ scheme at 2 GeV. Our comments are as follows:

At fixed $\beta$ , the two estimates from $\langle \bar { \psi } \psi \rangle _ { 1 } ^ { \overline { { \mathrm { M S } } } }$ and $\langle \bar { \psi } \psi \rangle _ { 2 } ^ { \overline { { \mathrm { M S } } } }$ obtained with the Wilson action are perfectly compatible within the statistical errors. • At $\beta = 6 . 0$ , the $\langle \bar { \psi } \psi \rangle _ { 2 } ^ { \overline { { \mathrm { M S } } } }$ ultaiit theClover n Wilon acins are perfectly compatible within two standard deviations of the statistical error; at $\beta = 6 . 2$ they are compatible within one standard deviation of the statistical error. This implies that discretization effects are under control. • Even within our increased statistical accuracy it is not possible to see a systematic dependence of the chiral condensate in the $\beta$ range 6.0-6.2. The results at $\beta = 6 . 4$ apparently show such a dependence, but, as we have argued in subsect. 5.2, this is more likely to be a finite size effect. Thus, we will again not take into consideration results at $\beta = 6 . 4$ . • As previously stated, the systematic error due to the NP renormalization constants is small in all $\langle \bar { \psi } \psi \rangle _ { 1 } ^ { \overline { { \mathrm { M S } } } }$ results and dominant in the $\langle \bar { \psi } \psi \rangle _ { 2 } ^ { \overline { { \mathrm { M S } } } }$ results. However we cannot exclude that the very small error in $Z _ { S }$ is underestimated. A conservative attitude consists is regarding the larger error of $Z _ { P }$ in the $\langle \bar { \psi } \psi \rangle _ { 2 } ^ { \overline { { \mathrm { M S } } } }$ results as more realistic.

In order to appreciate the effects of the NP renormalization on the chiral condensate, we now present the results at $\beta = 6 . 2$ obtained from the perturbative values of $Z _ { S }$ and $Z _ { P }$ . We use the results of Boosted PT (BPT) $\langle \bar { \psi } \psi \rangle _ { 2 } ^ { \overline { { { M S } } } } = - 0 . 0 1 8 5 ( 1 0 )$ wheres fom h Wln ction $\langle \bar { \psi } \psi \rangle _ { 1 } ^ { \overline { { { M S } } } } = - 0 . 0 1 5 8 ( 5 )$ and $\langle \bar { \psi } \psi \rangle _ { 2 } ^ { \overline { { M S } } } =$ $- 0 . 0 2 0 5 ( 7 )$ . It is obvious that the excellent compatibility of the corresponding results, obtained with NP renormalization constants, has been destroyed by PT.

# 5.4 Final result

From the discussion in the previous section on the results of Table 4 we can now derive our best estimate for the chiral condensate. Our best results are those at $\beta = 6 . 2$ (high enough UV cutoff without finite volume effects). Although all results at this $\beta$ value are compatible, the Clover estimate is in principle preferable, as it has $\mathcal { O } ( g _ { 0 } ^ { 2 } a )$ discretization errors. Thus we give as our best estimate the C62 result of Table 4:

<table><tr><td rowspan=1 colspan=1>Run</td><td rowspan=1 colspan=1>CHS</td><td rowspan=1 colspan=1>$b$\$</td><td rowspan=1 colspan=1>CAWI</td><td rowspan=1 colspan=1>$b$\qr}$</td></tr><tr><td rowspan=1 colspan=1>C60aC60bW60</td><td rowspan=1 colspan=1>2.98(8)3.04(7)2.40(5)</td><td rowspan=1 colspan=1>0.0150(3)(2)</td><td rowspan=1 colspan=1>3.9(1)4.1(1)3.01(7)</td><td rowspan=1 colspan=1>0.0141(5)(21)0.0146(3)(21)0.0152(4)(20)</td></tr><tr><td rowspan=1 colspan=1>C62W62</td><td rowspan=1 colspan=1>2.9(1)2.52(8)</td><td rowspan=1 colspan=1>0.0156(5)(2)</td><td rowspan=1 colspan=1>3.7(2)2.98(9)</td><td rowspan=1 colspan=1>0.0147(8)(16)0.0158(5)(16)</td></tr><tr><td rowspan=1 colspan=1>C64W64</td><td rowspan=1 colspan=1>3.5(2)2.9(1)</td><td rowspan=1 colspan=1>0.0172(7)(2)</td><td rowspan=1 colspan=1>4.2(2)3.2(1)</td><td rowspan=1 colspan=1>0.0187(9)(10)0.0179(8)(13)</td></tr></table>

Table 4   
The pure lattice quantities $C ^ { H S }$ , $C ^ { A W I }$ and the corresponding chiral condensates (in GeV3) at a scale $\mu = 2$ GeV, in the $\mathrm { \overline { { M S } } }$ scheme. The first error is statistical, the second comes from the NP renormalization constants. An overall minus sign has been suppressed in the results of the chiral condensate.

$$
\begin{array} { r l r } {  { \frac { 1 } { N _ { f } } \langle \bar { \psi } \psi \rangle ^ { \overline { { \mathrm { M S } } } } ( \mu = 2 \mathrm { G e V } ) = - ( 0 . 0 1 4 7 \pm 0 . 0 0 0 8 \pm 0 . 0 0 1 6 \pm 0 . 0 0 1 2 ) \mathrm { G e V } ^ { 3 } } } \\ & { } & \\ & { } & { = - [ ( 2 4 5 \pm 4 \pm 9 \pm 7 \mathrm { M e V } ) ] ^ { 3 } , \phantom { \frac { 1 } { 1 } \Bigg ( } \quad \mathrm { ( 4 . } } \end{array}
$$

where the first error is statistical, the second is due to the NP renormalization and the third due to the lattice calibration. For the first time we are able to estimate the systematic uncertainties in a reliable way. We have also seen that our analysis is in accordance (with improved statistical accuracy) with earlier lattice results [12,13]; c.f. eq. (40). The results cited here disagree by a factor of 2 with the recent lattice computation of ref. [34]. Similar discrepancies were seen between the quark mass estimates of this work and, say, ref. [9]. The two effects are not independent, since the GMOR relation connects the chiral condensate to the quark mass. Presumably, the results of ref. [34] are flawed by unreliable extrapolations of the lattice data obtained at several UV cutoff values to the continuum limit.

Our result is also in agreement with the one obtained from sum rules (see [14] and references therein):

$$
\begin{array} { r l r } & { } & { \displaystyle \frac { 1 } { N _ { f } } \langle \bar { \psi } \psi \rangle ^ { \overline { { \mathrm { M S } } } } ( \mu = 2 \mathrm { G e V } ) = - \left( 0 . 0 1 4 \pm 0 . 0 0 2 \right) \mathrm { G e V } ^ { 3 } } \\ & { } & { = - \left[ \left( 2 4 2 \pm 9 \right) \mathrm { M e V } \right] ^ { 3 } . } \end{array}
$$

For completeness we also give the quenched RGI result defined in eq. (14):

$$
\begin{array} { r c l } { { } } & { { } } & { { \displaystyle \frac { 1 } { N _ { f } } \langle \bar { \psi } \psi \rangle ^ { \mathrm { R G I } } = - \left( 0 . 0 0 8 8 \pm 0 . 0 0 0 5 \pm 0 . 0 0 1 \pm 0 . 0 0 0 7 \right) \mathrm { G e V ^ { 3 } } } } \\ { { } } & { { } } & { { } } \\ { { } } & { { } } & { { = - \left[ \left( 2 0 6 \pm 4 \pm 8 \pm 5 \mathrm { M e V } \right) \right] ^ { 3 } . } } \end{array}
$$

The errors correspond to those of eq. (43).

The reference scale $\mu = 2 \mathrm { G e V }$ is the standard lattice choice. In order to facilitate comparison with other determinations of the condensate, conventionally given at the scale $\mu = 1 6 \mathrm { e V }$ , we run our results to the latter scale, using NNLO perturbation theory. We find

$$
\begin{array} { r l r } {  { \frac { 1 } { N _ { f } } \langle \bar { \psi } \psi \rangle ^ { \overline { { \mathrm { M S } } } } ( \mu = 1 \mathrm { G e V } ) = - ( 0 . 0 1 2 4 \pm 0 . 0 0 0 7 \pm 0 . 0 0 1 4 \pm 0 . 0 0 1 0 ) \mathrm { G e V } ^ { 3 } } } \\ & { } & \\ & { } & { = - [ ( 2 3 1 \pm 4 \pm 8 \pm 6 \mathrm { M e V } ) ] ^ { 3 } , \qquad ( 4 6 \pm 7 \pm 4 \mathrm { G e V } ) } \end{array}
$$

to be compared to the sum rule result of ref. [14]:

$$
\begin{array} { r } { \displaystyle \frac { 1 } { N _ { f } } \langle \bar { \psi } \psi \rangle ^ { \overline { { \mathrm { M S } } } } ( \mu = 1 \mathrm { G e V } ) = - ( 0 . 0 1 2 \pm 0 . 0 0 2 ) \mathrm { G e V } ^ { 3 } } \\ { = - \left[ ( 2 2 9 \pm 9 ) \mathrm { M e V } \right] ^ { 3 } . } \end{array}
$$

# 6 Conclusions

We have computed the QCD chiral condensate from first principles in the framework of the lattice regularization with Wilson fermions. Our result has been obtained in the quenched approximation. We have been particularly careful in understanding and controlling all important sources of error (except for quenching): (i) large field configuration ensembles have been generated in order to minimize the statistical error; (ii) two actions (Wilson and tree-level Clover) and several gauge couplings have been used to control the discretization error; (iii) two lattice volumes (at one gauge coupling) ensure some control of the finite size effects; (iv) several methods of extracting the chiral condensate, equivalent only in the chiral limit, give us some confidence about the extrapolations to the chiral limit; (v) all necessary renormalizations have been performed non-perturbatively to avoid large tadpole contributions which distort lattice perturbative renormalization at 1-loop; (vi) the RI - $\mathrm { \overline { { M S } } }$ renormalization scheme matching and the RG running from the lattice scale $a ^ { - 1 }$ to the usual scale 2 GeV have been done in continuum perturbation theory at NNLO. Our final estimate for the chiral condensate is given in eq.(43). It is compatible with previous lattice determinations with a much smaller statistical error and a careful estimate of the systematic error. It is also compatible with the sum rule result. Thus we conclude that the chiral condensate is an order parameter of the spontaneous breaking of chiral symmetry in QCD with massless quarks.

# Acknowledgements

We are extremely grateful to J. Gasser, V. Gimenez, V. Lubicz, G. Martinelli, G.C. Rossi and M. Testa for many useful discussions. M.T. acknowledges the support of PPARC through grant GR/L22744.

# APPENDICES

# A The RI renormalization scheme

We shortly review the NP method for the renormalization constants of quark bilinear operators. The full discussion of the method can be found in ref. [25]. The most up to date results, which we use in the present work, can be found in ref. [28].

Given a quark bilinear $O _ { \Gamma } ^ { f } = \bar { \psi } \Gamma ( \lambda ^ { f } / 2 ) \psi$ $\Gamma$ is any Dirac matrix), we consider the operator insertion in the quark propagator

$$
G _ { O } ( p ) = \int d x _ { 1 } d x _ { 2 } \exp [ i p ( x _ { 1 } - x _ { 2 } ) ] \langle \psi ( x _ { 1 } ) O _ { \Gamma } ^ { f } ( 0 ) \bar { \psi } ( x _ { 2 } ) \rangle ~ ,
$$

from which we obtain the amputated Green function computed between offshell quark states of momentum $p$ in the Landau gauge:

$$
\Lambda _ { \cal O } ( p a ) = { \cal S } ( p a ) ^ { - 1 } { \cal G } _ { \cal O } ( p a ) { \cal S } ( p a ) ^ { - 1 } .
$$

The above quantity is computed non-perturbatively via Monte Carlo simulations in the Landau gauge. [25]. Any effects from Gribov copies are small [24]; spurious solutions [35] have not been considered. The renormalization $Z _ { O } ( \mu a , g _ { 0 } )$ of $O _ { \Gamma }$ is determined by the RI condition

$$
Z _ { \cal O } ( \mu a ) Z _ { q } ^ { - 1 } ( \mu a ) \mathrm { T r } ~ \mathrm { P } _ { \cal O } \Lambda _ { \cal O } ( p a ) \biggr | _ { p ^ { 2 } = \mu ^ { 2 } } = 1 ~ ,
$$

where $\mathrm { P } _ { O }$ is a projector chosen so that the above condition is satisfied at tree level [25] and $Z _ { q }$ is the wave function renormalization. Since the projected

amputated Green function can be computed non-perturbatively with Monte Carlo simulations, eq.(A.3) can be solved for $Z _ { O }$ at a fixed scale $\mu$ and cutoff $a$ provided that $Z _ { q }$ is known.

In order for the RI scheme to be compatible (at large scales $\mu$ ) with the WIs, $\textstyle Z _ { q }$ must be defined as [25]

$$
Z _ { q } ( \mu a ) = \left. - i \frac { 1 } { 1 2 } T r ( \frac { \partial S ( p a ) ^ { - 1 } } { \partial \not p } ) \right| _ { p ^ { 2 } = \mu ^ { 2 } } .
$$

The above definition is not convenient, as it requires differentiation with respect to the discrete variable $p$ . Thus, following ref. [25], we have opted for the definition

$$
Z _ { q } ^ { \prime } ( \mu a ) = \left. - i \frac { 1 } { 1 2 } \frac { T r \sum _ { \mu = 1 , 4 } \gamma _ { \mu } \sin ( p _ { \mu } a ) S ( p a ) ^ { - 1 } } { 4 \sum _ { \mu = 1 , 4 } \sin ^ { 2 } ( p _ { \mu } a ) } \right| _ { p ^ { 2 } = \mu ^ { 2 } } ,
$$

which, in the Landau Gauge, differs from $Z _ { q }$ by a finite term of order $\alpha _ { s } ^ { 2 }$ . The matching coefficient can be computed using continuum perturbation theory. From refs. [28,30] we quote up to order $\alpha _ { s } ^ { 2 }$ ,

$$
{ \frac { Z _ { q } } { Z _ { q } ^ { \prime } } } = 1 - { \frac { \alpha _ { s } ^ { 2 } } { \left( 4 \pi \right) ^ { 2 } } } \Delta _ { q } ^ { ( 2 ) } + \ldots ,
$$

where

$$
\Delta _ { q } ^ { ( 2 ) } = \frac { \left( N _ { c } ^ { 2 } - 1 \right) } { 1 6 N _ { c } ^ { 2 } } \ : \left( 3 + 2 2 N _ { c } ^ { 2 } - 4 N _ { c } N _ { f } \right)
$$

( $N _ { f } = 0$ in the quenched approximation). Note that at a typical scale $\mu \sim$ 2 GeV, $\Delta _ { q }$ represents a tiny correction, smaller than $1 \%$ . The renormalization constants used in the present work have been corrected by this factor (c.f. [28]).

# B RI/ $\mathbf { \overline { { M S } } }$ matching coefficient and RG evolution coefficient

The matching coefficient $\Delta Z ^ { \mathrm { R I / \overline { { M S } } } }$ which transforms the renormalized chiral condensate from the RI to the $\mathrm { \overline { { M S } } }$ scheme is that of the renormalization constant $Z _ { P }$ (or $Z _ { S }$ ). It has been calculated in PT to $\mathcal { O } ( \alpha _ { s } ^ { 2 } )$ (see refs. [28,30] for

details). The result is

$$
\Delta Z ^ { \mathrm { R I / \overline { { { \mathrm { M S } } } } } } = 1 + \frac { \alpha _ { s } ( \mu ) } { 4 \pi } C ^ { ( 1 ) } + \frac { \alpha _ { s } ^ { 2 } ( \mu ) } { ( 4 \pi ) ^ { 2 } } C ^ { ( 2 ) } + { \mathcal O } ( \alpha _ { s } ^ { 3 } ) ,
$$

where

$$
\begin{array} { l } { { \displaystyle C ^ { ( 1 ) } = \frac { 8 \left( N _ { c } ^ { 2 } - 1 \right) } { 4 N _ { c } } } } \\ { { \displaystyle C ^ { ( 2 ) } = \frac { \left( N _ { c } ^ { 2 } - 1 \right) } { 9 6 N _ { c } ^ { 2 } } \left( - 3 0 9 + 3 0 2 9 N _ { c } ^ { 2 } \right. } } \\ { { \left. \quad \quad - 2 8 8 \zeta _ { 3 } - 5 7 6 N _ { c } ^ { 2 } \zeta _ { 3 } - 3 5 6 N _ { c } N _ { f } \right) } } \end{array}
$$

and $\zeta _ { 3 } = 1 . 2 0 2 0 6 \cdots$ .

The evolution coefficient of the chiral condensate is derived from a standard RG analysis. Its running is determined by the evolution of the strong coupling $\alpha _ { s } ( \mu )$ (i.e. the Callan-Symanzik $\beta$ function) and the scalar operator (i.e. its anomalous dimension $\gamma _ { S } ( \alpha _ { s } )$ . In any scheme which respects vector symmetry (i.e. the vector WIs are valid) the renormalization of the scalar operator is the inverse of the renormalization of the quark mass 10 (see for example eq.(17)). Thus, we can readily use the results of refs. [36] obtained in the $\overline { { \mathrm { M S } } }$ scheme (with NDR regularization) at the NNLO. The Callan -Symanzik $\beta$ function and the quark mass anomalous dimension are given by

$$
\begin{array} { r l r } & { } & { \frac { \beta ( \alpha _ { s } ) } { 4 \pi } = \mu ^ { 2 } \frac { d } { d \mu ^ { 2 } } \left( \frac { \alpha _ { s } } { 4 \pi } \right) = - \displaystyle \sum _ { i = 0 } ^ { \infty } \beta _ { i } \left( \frac { \alpha _ { s } } { 4 \pi } \right) ^ { i + 2 } } \\ & { } & { \gamma _ { m } ( \alpha _ { s } ) = - 2 Z _ { m } ^ { - 1 } \mu ^ { 2 } \frac { d Z _ { m } } { d \mu ^ { 2 } } = \displaystyle \sum _ { i = 0 } ^ { \infty } \gamma _ { m } ^ { ( i ) } \left( \frac { \alpha _ { s } } { 4 \pi } \right) ^ { i + 1 } , } \end{array}
$$

where

$$
\begin{array} { l } { { \beta _ { 0 } { = } { \displaystyle { \frac { 1 1 } { 3 } } N _ { c } } - { \frac { 2 } { 3 } } N _ { f } } } \\ { { \beta _ { 1 } { = } { \frac { 3 4 } { 3 } } N _ { c } ^ { 2 } - { \frac { 1 0 } { 3 } } N _ { c } N _ { f } - { \frac { \left( N _ { c } ^ { 2 } - 1 \right) } { N _ { c } } } N _ { f } } } \\  { \beta _ { 2 } ^ { \overline { { \mathrm { M S } } } } { = } { \displaystyle { \frac { 2 8 5 7 } { 5 4 } } N _ { c } ^ { 3 } + { \frac { \left( N _ { c } ^ { 2 } - 1 \right) ^ { 2 } } { 4 N _ { c } ^ { 2 } } } N _ { f } - { \frac { 2 0 5 } { 3 6 } } \left( N _ { c } ^ { 2 } - 1 \right) N _ { f } } } \\ { { \displaystyle ~ - { \frac { 1 4 1 5 } { 5 4 } } N _ { c } ^ { 2 } N _ { f } + { \frac { 1 1 } { 1 8 } } { \frac { \left( N _ { c } ^ { 2 } - 1 \right) } { N _ { c } } } N _ { f } ^ { 2 } + { \frac { 7 9 } { 5 4 } } N _ { c } N _ { f } ^ { 2 } } } \end{array}
$$

$$
\begin{array} { l } { { \displaystyle \gamma _ { m } ^ { ( 0 ) } = 3 \frac { N _ { c } ^ { 2 } - 1 } { N _ { c } } } } \\ { { \displaystyle \gamma _ { m } ^ { ( 1 ) } = \frac { N _ { c } ^ { 2 } - 1 } { N _ { c } ^ { 2 } } \left( - \frac { 3 } { 4 } + \frac { 2 0 3 } { 1 2 } N _ { c } ^ { 2 } - \frac { 5 } { 3 } N _ { c } N _ { f } \right) } } \\ { { \displaystyle \gamma _ { m } ^ { ( 2 ) } = \frac { N _ { c } ^ { 2 } - 1 } { N _ { c } ^ { 3 } } \left[ \frac { 1 2 9 } { 8 } - \frac { 1 2 9 } { 8 } N _ { c } ^ { 2 } + \frac { 1 1 4 1 3 } { 1 0 8 } N _ { c } ^ { 4 } \right. } } \\ { { \displaystyle \left. + N _ { f } \left( \frac { 2 3 } { 2 } N _ { c } - \frac { 1 1 7 7 } { 5 4 } N _ { c } ^ { 3 } - 1 2 N _ { c } \zeta _ { 3 } - 1 2 N _ { c } ^ { 3 } \zeta _ { 3 } \right) - \frac { 3 5 } { 2 7 } N _ { c } ^ { 2 } N _ { f } ^ { 2 } \right] } \ . }  \end{array}
$$

Note that the scheme dependence settles in only at the two-loop order for $\gamma _ { m } ( \alpha _ { s } )$ and three-loop order for $\beta ( \alpha _ { s } )$ . Finally, to NNLO the evolution coefficient of the chiral condensate is [36]

$$
\begin{array} { l } { { \displaystyle c _ { S } ^ { \overline { { \mathrm { M S } } } } \left( \mu \right) = \alpha _ { s } \left( \mu \right) ^ { \overline { { \gamma } } _ { S } ^ { \left( 0 \right) } } \left\{ 1 + \frac { \alpha _ { s } } { 4 \pi } \left( \overline { { \gamma } } _ { S } ^ { \left( 1 \right) } - \overline { { \beta } } _ { 1 } \overline { { \gamma } } _ { S } ^ { \left( 0 \right) } \right) \right. } \left. \left( \mathrm { B . 6 } \right) \right. }  \\ { { \displaystyle \qquad + \left. \frac { 1 } { 2 } \left( \frac { \alpha _ { s } \left( \mu \right) } { 4 \pi } \right) ^ { 2 } \left[ \left( \overline { { \gamma } } _ { S } ^ { \left( 1 \right) } - \overline { { \beta } } _ { 1 } \overline { { \gamma } } _ { S } ^ { \left( 0 \right) } \right) ^ { 2 } + \overline { { \gamma } } _ { S } ^ { \left( 2 \right) } + \overline { { \beta } } _ { 1 } ^ { 2 } \overline { { \gamma } } _ { S } ^ { \left( 0 \right) } - \overline { { \beta } } _ { 1 } \overline { { \gamma } } _ { S } ^ { \left( 1 \right) } - \overline { { \beta } } _ { 2 } \overline { { \gamma } } _ { S } ^ { \left( 0 \right) } \right] \right\} \ : , } } \end{array}
$$

with $\overline { { \beta } } _ { i } = \beta _ { i } / \beta _ { 0 }$ and $\overline { { \gamma } } _ { S } ^ { i } = \gamma _ { S } ^ { ( i ) } / \left( 2 \beta _ { 0 } \right) = - \gamma _ { m } ^ { ( i ) } / \left( 2 \beta _ { 0 } \right)$ The running coupling $\alpha _ { s } ( \mu )$ at NNLO in $\mathrm { \overline { { M S } } }$ scheme is given by

$$
\begin{array} { r l r } {  { \frac { \alpha _ { s } ^ { \overline { { \mathrm { M S } } } } } { 4 \pi } ( q ^ { 2 } ) = \frac { 1 } { \beta _ { 0 } \ln ( q ^ { 2 } ) } - \frac { \beta _ { 1 } } { \beta _ { 0 } ^ { 3 } } \frac { \ln \ln ( q ^ { 2 } ) } { \ln ^ { 2 } ( q ^ { 2 } ) } } } \\ & { } & { + \frac { 1 } { \beta _ { 0 } ^ { 5 } \ln ^ { 3 } ( q ^ { 2 } ) } ( \beta _ { 1 } ^ { 2 } \ln ^ { 2 } \ln ( q ^ { 2 } ) - \beta _ { 1 } ^ { 2 } \ln \ \ln ( q ^ { 2 } ) + \beta _ { 2 } ^ { \overline { { \mathrm { M S } } } } \beta _ { 0 } - \beta _ { 1 } ^ { 2 } ) \ , } \end{array}
$$

where $q ^ { 2 } = ( \mu / \Lambda _ { Q C D } ^ { \overline { { \mathrm { M S } } } } ) ^ { 2 }$ $N _ { f } = 0$ in all formulae. The QCD scale in the quenched approximation has been set to $\Lambda _ { Q C D } ^ { \overline { { \mathrm { M S } } } } = 0 . 2 5 1 \pm 0 . 0 2 1 \ \mathrm { G } \ i$ V; seeref. [37..

# C Lattice WIs with Clover fermions

In this Appendix we discuss the modifications which are necessary for a computation of the chiral condensate with the tree level Clover improved action. This involves the discussion of several technical details. At tree-level, the Clover action is defined as follows [10]:

$$
S _ { c } = S _ { f } - a ^ { 4 } \sum _ { x , \mu , \nu } \frac { a } { 4 } i g _ { 0 } \bar { \psi } ( x ) \sigma _ { \mu \nu } F _ { \mu \nu } ( x ) \psi ( x ) ,
$$

with $F _ { \mu \nu } ( x )$ the clover-leaf discretization of the field tensor. At this order all $\mathcal { O } ( a g _ { 0 } ^ { 2 n } \ln ^ { n } a )$ terms, which are effectively of $\mathcal O ( a )$ in the scaling limit ( $g _ { 0 } ^ { 2 } \sim$ $1 / \ln { a } )$ , are eliminated from correlation functions. At leading-log level the improvement of local operators can be expressed as a rotation of the fermion fields [11]:

$$
{ \cal O } _ { \Gamma } ^ { I } ( x ) = \bar { \psi } ^ { R } ( x ) \Gamma \psi ^ { R } ( x ) ,
$$

where the rotated fields are defined through

$$
\begin{array} { l } { \displaystyle \psi ^ { R } ( \boldsymbol { x } ) = [ 1 - \frac { a \stackrel {  } { \mathbb { P } } } { 2 } ] \psi ( \boldsymbol { x } ) } \\ { \displaystyle \bar { \psi } ^ { R } ( \boldsymbol { x } ) = \bar { \psi } ( \boldsymbol { x } ) [ 1 + \frac { a \stackrel {  } { \mathbb { P } } } { 2 } ] } \end{array}
$$

and $D _ { \ l }$ and $D /$ are symmetric lattice discretizations of the covariant derivatives (see ref. [11] for their definition). The improved operators $O ^ { I }$ differ from the original ones by terms proportional to the cutoff. When $O ^ { I }$ is inserted in a Green function, these extra terms combine with the $a ^ { - 1 }$ UV divergences to give finite contributions. Consequently, Clover renormalization constants differ from the Wilson ones by finite terms.

We denote by $S ^ { I } ( x - y ) = \langle \psi ( x ) \psi ( y ) \rangle$ the quark propagator obtained by solving the Dirac equation of the Clover action. The tree-level $\mathcal O ( a )$ improved quark propagator is $\left. \psi ^ { R } ( x ) \psi ^ { R } ( y ) \right.$ ; in terms of $S ^ { I } ( x - y )$ it is given by (see ref. [19] for details):

$$
\langle \psi ^ { R } ( x ) { \bar { \psi } } ^ { R } ( y ) \rangle = S ^ { e f f } ( x - y ) + \frac { a } { 2 } \delta ( x - y ) + { \mathcal O } ( a ^ { 2 } ) ,
$$

where the effective rotated propagator $S ^ { e f f } ( x - y )$ is defined as

$$
S ^ { e f f } ( x - y ) = \left[ 1 - \frac { a } { 2 } \stackrel { \right. } { \mathcal { P } } ( x ) \right] S ^ { I } ( x - y ) \left[ 1 + \frac { a } { 2 } \stackrel { \left. \right]} { \mathcal { P } } ( y )  .
$$

This is the Clover propagator we use in our computations. The reason it contains an $\scriptstyle { \mathcal { O } } ( a ^ { 2 } )$ term in its definition is that it is readily computable in numerical simulations [38]. Thus, it has been extensively used in several Clover improved lattice QCD computations of on-shell matrix elements. In computing correlation functions between external on-shell hadron states, the $\delta$ -function never contributes. To illustrate the point, consider a two-point correlation function of two bilinear improved operators $O _ { j } ^ { I } ( x ) = \psi ^ { R } ( x ) \Gamma _ { j } \psi ^ { R } ( x )$ (where $j = 1 , 2$ )given by

$$
\begin{array} { l } { { \langle { \cal O } _ { 1 } ^ { I } ( x ) { \cal O } _ { 2 } ^ { I } ( y ) \rangle = \langle { \cal O } _ { 1 } ^ { I } ( x ) { \cal O } _ { 2 } ^ { I } ( y ) \rangle _ { n f \rangle c t } \nonumber } } \\ { { \displaystyle - \frac { a } { 2 } \delta ( x - y ) \mathrm { T r } \left[ \Gamma _ { 1 } \Gamma _ { 2 } S ^ { e f f } ( x - y ) \right] - \frac { a } { 2 } \delta ( x - y ) \mathrm { T r } \left[ \Gamma _ { 2 } \Gamma _ { 1 } S ^ { e f f } ( x - y ) \right] \ : , } } \end{array}
$$

where

$$
\left. { \cal O } _ { 1 } ^ { I } ( x ) { \cal O } _ { 2 } ^ { I } ( y ) \right. _ { n \bar { p } e t } = - \left. \mathrm { T r } ~ \left[ S ^ { e f f } ( x - y ) \Gamma _ { 2 } S ^ { e f f } ( y - x ) \Gamma _ { 1 } \right] \right.
$$

(i.e. the subscript $n \not | D c t$ stands for "no $D /$ contact terms"). Note that the nΦct correlation functions are traces of the effective propagator $S ^ { e f f } ( x - y )$ .However, in the case of the chiral condensate, c.f. eq. (29), the space-time points of the two-point correlation function are allowed to coincide. The on-shell matrix elements argument is then invalidated and the complete contribution of eq. (C.6) must be used with the $\delta$ -function giving rise to contact terms.

The above improvement scheme is by no means unique. Other choices of quark field rotations are also possible, as discussed in refs. [19,22,39]. The choice made here is not necessarily the most convenient for improving lattice WIs. As previously pointed out, it is however, the one implemented in all APE simulations with tree-level Clover improvement. Thus, we need to examine how $\mathcal O ( a )$ -improved WIs can be obtained for this choice of rotation (i.e. this choice of improved operators).

Upon performing an axial variation of the improved pseudoscalar operator $P ^ { I g } ( 0 )$ with the theory defined in terms of the Clover action, we obtain

$$
\begin{array} { r l r } & { } & { \frac { \delta \langle P ^ { I g } ( 0 ) \rangle } { \delta \alpha _ { A } ^ { f } ( x ) } = 0 \Longleftrightarrow } \\ & { } & { i \langle \frac { \delta P ^ { I g } ( 0 ) } { \delta \alpha _ { A } ^ { f } ( x ) } \rangle = a ^ { 4 } \underset { \mu } { \sum } \nabla _ { x } ^ { \mu } \langle A _ { \mu } ^ { f } ( x ) P ^ { I g } ( 0 ) \rangle - a ^ { 4 } 2 m _ { 0 } \langle P ^ { f } ( x ) P ^ { I g } ( 0 ) \rangle } \\ & { } & { \qquad - a ^ { 4 } \langle X ^ { f } ( x ) P ^ { I g } ( 0 ) \rangle \ . } \end{array}
$$

This is eq. (18) with the difference that the variation of the Clover term in the action modifies the expression for $X ^ { f }$ by a term proportional to $a \bar { \psi } \sigma _ { \mu \nu } \tilde { F } _ { \mu \nu } \psi$ (c.f. eq. (C.1); $\tilde { F } _ { \mu \nu }$ is the dual field tensor). The above WI is not $\mathcal O ( a )$ improved for three reasons: (i) the correlation functions contain unimproved operators; (ii) the asymmetric derivative is not improved ( $\nabla _ { x } ^ { \mu } = \partial _ { x } ^ { \mu } + { \mathcal { O } } ( a ) )$ ; (iii) the variation of the improved operator $P ^ { I } { \boldsymbol { g } }$ (the l.h.s. of the above) does not yield the vacuum expectation value of the improved scalar density $\langle S ^ { I 0 } ( 0 ) \rangle$ (c.f. eq. (18)). The conclusion is that improved WIs cannot be obtained by the standard procedure of performing symmetry transformations on improved operators and the Clover improved action. Instead, we simply impose the validity of WIs (such as that of eq. (18)) with symmetric derivatives and correlation functions of improved operators. Since all equations of sect. 4, implemented at fixed inverse coupling $\beta$ , suffer from $\mathcal O ( a )$ corrections, we can add by hand the $\mathcal O ( a )$ terms which render the correlation functions improved. Consistency requires that the mass subtraction $\overline { { m } } ( m _ { 0 } )$ is also modified by $\mathcal O ( a )$ terms so that $( m _ { 0 } - { \overline { { m } } } ^ { I } ) P ^ { I f }$ is $\mathcal O ( a )$ improved. Consequently, the mass difference $\left( m _ { 0 } - { \overline { { m } } } ^ { I } \right)$ vanishes at a critical point $m _ { C } ^ { I }$ . In other words, we define $2 \rho ^ { I }$ as:

$$
2 \rho ^ { I } = \frac { 2 ( m _ { 0 } - \overline { { { m } } } ^ { I } ) } { Z _ { A ^ { I } } } = \frac { \int d ^ { 3 } x \overline { { { \nabla } } } _ { x } ^ { 0 } \langle A _ { 0 } ^ { I f } ( x ) P ^ { I f } ( 0 ) \rangle } { \int d ^ { 3 } x \langle P ^ { I f } ( x ) P ^ { I f } ( 0 ) \rangle } ,
$$

where $\overline { { \nabla } } _ { x } ^ { 0 }$ is the symmetric lattice derivative. This equation, with $Z _ { A ^ { I } }$ computed, say, from improved WIs [19,24,40] is a definition of the improved quantity $( m _ { 0 } - { \overline { { m } } } ^ { I } )$ . With this definition, this mass difference is renormalized by $Z _ { P ^ { I } } ^ { - 1 }$ (the improved version of eq. (23)).

We can now write down the improved WI for the chiral condensate. It is the WI (29) expressed in terms of improved operators and masses. Also, contact terms must be properly taken into account according to eq. (C.6). We obtain

$$
\begin{array} { r l r } { \displaystyle \frac { 1 } { N _ { f } } \langle \bar { \psi } \psi \rangle _ { s u b } ^ { I } = \displaystyle \operatorname* { l i m } _ { m _ { 0 }  m _ { C } ^ { I } } 2 ( m _ { 0 } - \overline { { m } } ^ { I } ) \int d ^ { 4 } x \langle P ^ { I f } ( x ) P ^ { I f } ( 0 ) \rangle } & { { \mathrm { ( C . 1 0 ) } } } & \\ { \displaystyle = \operatorname* { l i m } _ { m _ { 0 }  m _ { C } ^ { I } } 2 ( m _ { 0 } - \overline { { m } } ^ { I } ) [ \int d ^ { 4 } x \langle P ^ { I f } ( x ) P ^ { I f } ( 0 ) \rangle _ { n p \bar { j } c t } - a \mathrm { T r } \langle S ^ { e f f } ( 0 ) \rangle ] ~ . } & { } & \end{array}
$$

Apart from these modifications, all other expressions for the chiral condensate (see sect. 5) remain valid, provided all quantities are improved.

Similar arguments can be used in the case of vector WIs. It is not possible to obtain Clover-improved WIs by the standard vector variations of improved operators. Instead, we impose the validity of WI (16) with all correlation functions expressed in terms of improved operators. This implies an improved definition of the subtracted mass $( m _ { 0 } - m _ { C } ^ { I } )$ , so that the renormalization constants $Z _ { S ^ { I } }$ and $Z _ { m } ^ { I }$ obey eq. (17). Moreover, as stressed repeatedly in this work, the simple relationship between the quark mass and the hopping parameter of eq. (35), satisfied to $\mathcal O ( a )$ , needs to be improved. As this is not a straightforward task with the rotated fields of eq.(C.3), we have avoided use of eq.(35) when deriving our Clover results. This point has been explained also in ref. [9], where an improved version of eq. (35) valid for a different rotation of the quark fields is given.

A final word of caution is in place here. The relationships (17) and (23) constrain the renormalization constants involved, so that the WIs are recovered in the continuum. These constraints relate not only the anomalous dimensions, but also the finite parts of these renormalization constants. The determination of the finite part of a renormalization constant depends on the choice of lattice operator 11 and renormalization condition. There are choices for which the WIs are not satisfied, and consequently eqs. (17) and (23) are violated by finite terms. For example, in 1-loop perturbation theory with the Clover action, $Z _ { m } ^ { I }$ has been obtained by renormalizing the Clover quark propagator $\langle \psi \bar { \psi } \rangle$ (in the $\mathrm { \overline { { M S } } }$ scheme) [21]. Compared to the result for $Z _ { S ^ { I } }$ obtained in ref. [22], we see that eq. (17) is violated by the finite term arising from the quark field rotations in $S ^ { I }$ [9]. Agreement would have been reached, had the rotated propagator $\langle \psi ^ { R } \psi ^ { R } \rangle$ been used to obtain $Z _ { m } ^ { I }$ . The discrepancy is not due to the choice of scheme, but to the choice of an unrotated quark propagator. On the other hand, the renormalization constants used in this work, obtained in the RI scheme with the NP method of [25], are perfectly compatible with the Clover WIs.

# References

[1] J.F. Gunion P.C. McNamee and M.D. Scadron, Nucl. Phys. B123 (1977) 445 ; H. Sazdjian, Nucl. Phys. B129 (1977) 319 ; N.H. Fuchs and H. Sazdjian, Phys. Rev. D18 (1978) 889; N.H. Fuchs and M.D. Scadron, Phys. Rev. D20 (1979) 2421; M.D. Scadron, Rep. Prog. Phys 44 (1981) 213; J.Phys. G7 (1981) 1325; R. Dashen Phys. Rev. 183 (1969) 1245; J. Gasser and H. Leutwyler, Phys, Rep. C87 (1982) 77.   
[2] J. Stern, H. Sazdjan and N.H. Fuchs, Phys. Rev. D47 (1993) 3814; M. Knecht et al., Phys. Lett. B313 (1993) 229; J. Stern, hep-ph/9712438 and hep-ph/9801282.   
[3] T. Banks and A. Casher, Nucl. Phys. B168 (1980) 103.   
[4] I.M. Barbour et al., Phys. Lett. B136 (1984) 80; I.M. Barbour et al., Phys. Lett. B158 (1985) 61; G. Salina and A. Vladikas Phys. Letts. B249 (1990) 119; G. Salina and A. Vladikas Nucl.Phys. B348 (1991) 210.   
[5] K. Jansen et al.,Phys. Lett. B 372 (1996) 275; M. Lüscher, S. Sint, R. Sommer and P. Weisz, Nucl. Phys. B 478 (1996) 365; M. Lüscher et al., Nucl. Phys. B491 (1997) 323.   
[6] C. Bernard and M. Golterman, Phys. Rev. D46 (1992) 853; S.R. Sharpe, Phys. Rev. D46 (1992) 3146.   
[7] S. Kim and D.K Sinclair, Phys. Rev. D52 (1995) 2614; R.D. Mawhinney, Nucl Phys. B(Proc. Suppl)47 (1996) 557; S. Kim and S. Ohta, Nucl Phys. B(Proc. Suppl)63 (1998) 185.   
[8] C.R. Allton, V. Giménez, L. Giusti and F. Rapuano, Nucl. Phys. B 489, (1997) 427.   
[9] V. Giménez, L. Giusti, F. Rapuano and M. Talevi, hep-lat/9801028 to appear in Nucl. Phys. B.   
[10] B. Sheikholeslami and R. Wohlert, Nucl. Phys. B259 (1985) 572.   
[11] G.Heatlie et al, Nucl. Phys. B352 (1991) 266.   
[12] L. Maiani and G. Martinelli, Phys.Lett. B178 (1986) 265.   
[13] D. Daniel et al., Phys. Rev. D46 (1992) 3130.   
[14] H.G. Dosch and S. Narison Phys. Lett. B417 (1998) 173.   
[15] L.H. Karsten and J. Smit, Nucl.Phys. B183 (1981) 103; M. Bochicchio et al., Nucl.Phys. B262 (1985) 331.   
[16] M. Gell-Mann, R. J. Oakes and B. Renner, Phys. Rev. 175 (1968) 2195.   
[17] E. Laermann, Nucl. Phys. B(Proc.Suppl.)63 (1998) 114.   
[18] K.G.Wilson, Phys. Rev. D10 (1974) 2445; K.G.Wilson, in "New Phenomena in Subnuclear Physics" (Erice 1975), ed. A.Zichichi (New York, Plenum, 1975).   
[19] M. Crisafulli, V. Lubicz and A. Vladikas Eur.Phys.J. C4 (1998) 145.   
[20] B. Meyer and C. Smith, Phys. Lett. B123 (1983) 62; G. Martinelli and Y.C. Zhang, Phys.Lett. B123(1983) 433; Phys. Lett. B125 (1983) 77.   
[21] E. Gabrielli et al., Nucl.Phys. B362 (1991) 475.   
[22] A. Borrelli, R. Frezzotti, E. Gabrielli and C. Pittori, Nucl. Phys. B409 (1993) 382.   
[23] G. Parisi, in "High Energy Physics - 1980", Proceedings of the XXth International Conference, Madison, Wisconsin, eds. L. Durand and L. G. Pondrom (American Institute of Physics, New York, 1981); G.P. Lepage and P.B. Mackenzie, Phys. Rev. D48 (1993) 2250.   
[24] M. Paciello, S. Petrarca, B. Taglienti and A. Vladikas, Phys. Lett. B341 (1994) 187.   
[25] G.Martinelli et al, Nucl.Phys.B445(1995)81.   
[26] K. Symanzik, in Mathematical Problems in Theoretical Physics, Springer Lecture Notes vol.153 (1982) 47; eds. R. Schrader, R. Seiler and D.A. Uhlenbrock; K. Symanzik, Nucl. Phys. B226 (1983) 187; K. Symanzik, Nucl. Phys. B226 (1983) 205.   
[27] M. Göckeler et al., Nucl. Phys. B(Proc.Suppl) 47 (1996) 493.   
[28] V. Giménez, L. Giusti, F. Rapuano, M. Talevi, hep-lat/9806006 to appear in Nucl. Phys. B.   
[29] A. Donini et al., Phys. Lett. B360 (1996) 83; M. Crisafulli et al., Phys. Lett. B369 (1996) 325; A. Donini, et al., Nucl.Phys.B(Proc.Suppl.) 53 (1997) 883; JLQCD Collaboration, S. Aoki et al., Nucl.Phys.B(Proc.Suppl.) 53 (1997) 349; Phys. Rev. Lett. 81 (1998) 1778; L. Conti et al., Phys. Lett. B241 (1998) 273; C.R. Allton et al., hep-lat/9806016.   
[30] E. Franco and V. Lubicz, hep-ph/9803491.   
[31] M. Testa, JHEP 04(1998)002   
[32] G. Curci, Phys. Lett. B167 (1986) 265.   
[33] P. Lacock and C. Michael, Phys. Rev D52 (1995) 5213.   
[34] R. Gupta and T. Bhattacharya, Phys. Rev. D55 (1997) 7203.   
[35] L. Giusti, Nucl. Phys. B 498 (1997) 331; L. Giusti et al. Phys. Lett. B432 (1998) 196.   
[36] J.A.M. Vermaseren, S.A. Larin and T. van Ritbergen, Phys. Lett. B405 (1997) 327; T. van Ritbergen, J.A.M. Vermaseren and S.A. Larin, Phys. Lett. B400 (1997) 379; K.G. Chetyrkin, Phys. Lett. B404 (1997) 161.   
[37] S. Capitani et al., Nucl. Phys. B(Proc.Suppl) 63 (1998) 153.   
[38] G. Martinelli, C.T. Sachrajda, G. Salina and A. Vladikas, Nucl. Phys. B378 (1992) 591; Nucl. Phys. B397 (1993) 479.   
[39] G.Martinelli, C.T.Sachrajda and A.Vladikas, Nucl. Phys. B358 (1991) 212.   
[40] G. Martinelli, S. Petrarca, C.T. Sachrajda and A. Vladikas, Phys. Lett. B311 (1993) 241; Phys. Lett. B317 (1993) 660.