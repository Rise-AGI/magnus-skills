# A comparison between the Split Step Fourier and Finite-Difference method in analysing the soliton collision of a type of Nonlinear Schrödinger equation found in the context of optical pulses

Luke Taylor University of Cape Town Department of Mathematics and Applied Mathematics Cape Town, South Africa tylchr011@uct.ac.za

# Abstract

In this report a type of Schrödinger Equation which is found in the context of optical pulses is analysed using the Split Step and Finite Difference method. The investigation shows interesting dynamics regarding certain values for parameter $S$ as well as a comparison between the two numeric schemes demonstrating the Split Step to be superior for this problem.

# 1 Introduction

Optical pulses in media with saturation nonlinearity properties are modelled by the nonlinear Schrödinger (NLS) Equation [Kato, 1989] with the following nonlinearity [Zemlyanaya and Alexeeva, 2011].

$$
i \frac { \partial \psi } { \partial t } + \frac { 1 } { 2 } \frac { \partial ^ { 2 } \psi } { \partial x ^ { 2 } } + \frac { | \psi | ^ { 2 } \psi } { 1 + S | \psi | ^ { 2 } } = 0
$$

This is a Partial Differential Equation $( P D E )$ as it describes a relation of $\psi$ in regards to change over time and space. A solution to this equation is a so called soliton of the form.

$$
\psi ( x , t ) = \frac { 2 \sqrt { 2 } e ^ { \sqrt { 2 } x } } { 1 + ( \frac { 3 } { 2 } - 2 S ) e ^ { 2 \sqrt { 2 } x } } e ^ { i t + i v }
$$

The aim of this report is to study the collision of two such solitons [Serkin and Hasegawa, 2000] with varying values for $S$ and $v$ by numerically advancing from an initial configuration (the sum of two solitons) via the Split-Step [Weideman and Herbst, 1986] and Finite Difference [Delfour et al., 1981] methods.

# 2 Split Step Method

The Split Step method is a pseudo-spectral numerical method used to solve nonlinear $P D E s$ such as $P D E 1$ The method works by splitting the equation into a nonlinear and linear part.

$$
i \frac { \partial \psi } { \partial t } + \frac { | \psi | ^ { 2 } \psi } { 1 + S | \psi | ^ { 2 } } = 0
$$

$$
i \frac { \partial \psi } { \partial t } + \frac { 1 } { 2 } \frac { \partial ^ { 2 } \psi } { \partial x ^ { 2 } } = 0
$$

Both these equations are treated separately. The soliton is advanced in time by taking a small time step $\tau$ for both solutions. For the linear solution however $\psi$ needs to be Fourier transformed with the solution being advanced in Fourier space before inverse Fourier transforming back to the time domain.

# 2.1 Solving the nonlinear part

In order to solve equation 3 we multiply both sides by $i$ which gives us.

$$
\frac { \partial \psi } { \partial t } = i \frac { | \psi | ^ { 2 } } { 1 + S | \psi | ^ { 2 } } \psi
$$

Next we can notice that $| \psi | ^ { 2 }$ is a scalar and hence equation 6 is a first order differential equation which has the following analytical solution.

$$
\psi ( x , t ) = \psi ( x , t _ { 0 } ) e ^ { \frac { | \psi | ^ { 2 } } { 1 + S | \psi | ^ { 2 } } t }
$$

# 2.2 Solving the linear part

The key insight is to write $\psi ( x , t )$ as a Fourier Series.

$$
\psi ( x , t ) = \sum _ { i = - \infty } ^ { \infty } \hat { \psi } _ { n } ( t ) e ^ { \frac { 2 \pi i n } { L } x }
$$

Now subbing expression 8 into equation 4, doing some algebraic manipulation and rearranging the terms we get.

$$
\frac { \partial \hat { \psi } _ { n } } { \partial t } = - i 2 ( \frac { \pi n } { L } ) ^ { 2 } \hat { \psi } _ { n }
$$

This is again a first order differential equation which has the following solution.

$$
\hat { \psi } _ { n } ( t ) = \hat { \psi } _ { n } ( t _ { 0 } ) e ^ { - i 2 ( \frac { \pi n } { L } ) ^ { 2 } } t
$$

# 2.3 Finite Difference

The Finite Difference method works by approximating the derivatives in the expression with finite differences. In our $P D E$ we have $\frac { { \partial } { \psi } } { { \partial } t }$ and $\frac { \partial ^ { 2 } \psi } { \partial x ^ { 2 } }$ that need to be approximated via finite differences. The way $\frac { { \partial } { \psi } } { { \partial } t }$ is approximated determines what type of Finite Difference scheme is used which has various implications with regards to accuracy, stability and implementation.

1. 1) The Forward Difference is an explicit scheme which means that the solution at each point at the latest time level can be expressed through the solutions of the previous time levels. Although this simplifies the implementation the scheme suffers from stability issues. In fact for the $P D E 1$ it can be shown that the Forward Difference has an exponential growth in error using theVon Neumann stability analysis.

2. The Backwards Difference is an implicit scheme which means that a system of equations has to be solved in order to compute the solution at the next time level which makes the implementation non-trivial. However, this method has the superior property that it does not suffer from stability issues.

3. The Central Difference method has the advantage over both the Forward Difference and Backwards Difference in regards to the accuracy as the error is of $O ( \tau ^ { 2 } )$ compared to $O ( \tau )$ of the other methods. This means that the total error of the $P D E 1$ is of $O ( \tau ^ { 2 } + h ^ { 2 } )$ where $\tau$ is the time step and $h$ is the space step. This method suffers from the same shortcoming as the Forward Difference method: Stability issues.

The Central Difference was opted to solve the $P D E 1$ as it has good accuracy and although it suffers from stability issues it is stable for certain parameters for $\tau$ and $h$ shown in section 2.5. This method was chosen over the others schemes as the Forward Difference is provably unstable and the Backwards Difference was not investigated due to the nontrivial implementation and time overhead of solving a system of equations at every time step.

# 2.4 Algorithm

The Central Difference depends on the last two time solutions and hence the first time solution was approximated via the Forward Difference method. The Finite Difference was henceforth implemented as follows.

1Define the initial solution at $t = 0$ : $\psi ( x , t _ { 0 } ) = f ( x )$ 2. Approximate the solution at $t \ : = \ : 1$ using the Forward Difference: $\psi ( x , t _ { 1 } ) = F D ( \psi ( x , t _ { 0 } ) )$ 3. Now the Central Difference can be deployed to approximate the solution at $\begin{array} { r l r l } { t } & { { } = } & { i } \end{array}$ : $\mathbf { \bar { \psi } } \psi ( x , t _ { i } ) \mathbf { \bar { \Psi } } = \mathbf { \bar { \Psi } }$ $\mathbf { \bar { \it C D } } ( \psi ( x , t _ { i - 1 } ) , \psi ( x , t _ { i - 2 } ) )$

where $F D$ and $C D$ are the Finite Difference and Central Difference schemes respectively which are defined as:

$$
\psi _ { j , k + 1 } = \tau i ( 0 . 5 \psi _ { x x } + A _ { j , k } \psi _ { j , k } ) + \psi _ { j , k }
$$

$$
\psi _ { j , k + 1 } = 2 \tau i ( 0 . 5 \psi _ { x x } + A _ { j , k } \psi _ { j , k } ) + \psi _ { j , k - 1 }
$$

where

$$
\psi _ { x x } = \frac { \psi _ { j - 1 , k } - 2 \psi _ { j , k } + \psi _ { j + 1 , k } } { h ^ { 2 } }
$$

$$
A _ { j , k } = \frac { \vert \psi _ { j , k } \vert ^ { 2 } } { 1 + S \vert \psi _ { j , k } \vert ^ { 2 } }
$$

# 2.5 Stability analysis

In this section the stability of using the Central Difference method for $P D E 1$ is analysed to be able to make good choices for the parameters $\tau$ and $h$ . The Von Neumann stability analysis [Keller and Isaacson, 1994] is used to make sense of the stability.

Let $\psi _ { j , k } ~ = ~ \alpha ^ { k } e ^ { i \beta j }$ where $\psi _ { j , k }$ is the approximation at time $k$ at point $j$ .The right hand side of the equation is an arbitrary Fourier Mode. If the coefficient $| \alpha | > 1$ this implies the solution has an exponential growing error in time and hence it is required that $| \alpha | \le 1$ .

We assume the nonlinear term of $P D E 1$ to be negligible and perform the analysis on the linearised version.

$$
i \frac { \partial \psi } { \partial t } + \frac { 1 } { 2 } \frac { \partial ^ { 2 } \psi } { \partial x ^ { 2 } } = 0
$$

which has the Central Difference formula

$$
\psi _ { j , k + 1 } = \frac { \tau i } { h ^ { 2 } } ( \psi _ { j - 1 , k } - 2 \psi _ { j , k } + \psi _ { j + 1 , k } ) + \psi _ { j , k - 1 }
$$

Substituting $\psi _ { j , k } = \alpha ^ { k } e ^ { i \beta j }$ into this expression we get

$$
\begin{array} { c } { { \alpha ^ { k + 1 } e ^ { i \beta \beta j } = \displaystyle \frac { T i } { \hbar ^ { 2 } } ( \alpha ^ { k } e ^ { i \beta ( j - 1 ) } - 2 \alpha ^ { k } e ^ { i \beta j } + \alpha ^ { k } e ^ { i \beta ( j + 1 ) } ) + \alpha ^ { k - 1 } e ^ { i \beta j } } } \\ { { \alpha = \displaystyle \frac { T i } { \hbar ^ { 2 } } ( e ^ { - i \beta } - 2 + e ^ { i \beta } ) + \displaystyle \frac { 1 } { \alpha } } } \\ { { \alpha ^ { 2 } - \displaystyle \frac { T i } { \hbar ^ { 2 } } ( 2 c o s \beta - 2 ) \alpha - 1 = 0 } } \\ { { \alpha ^ { 2 } + \displaystyle \frac { 4 \pi i } { \hbar ^ { 2 } } ( s i n ^ { 2 } ( \frac { \beta } { 2 } ) ) \alpha - 1 = 0 } } \\ { { \alpha = - \displaystyle \frac { 2 \tau i } { \hbar ^ { 2 } } s i n ^ { 2 } ( \frac { \beta } { 2 } ) \pm \sqrt { - \frac { 4 \tau ^ { 2 } } { \hbar ^ { 4 } } s i n ^ { 4 } ( \frac { \beta } { 2 } ) + 1 } } } \end{array}
$$

The discriminant can either be positive or negative. It will be sufficient to consider the discriminant to be positive as

$$
| \alpha | ^ { 2 } = ( - \frac { 2 \tau i } { h ^ { 2 } } s i n ^ { 2 } ( \frac { \beta } { 2 } ) ) ) ^ { 2 } - \frac { 4 \tau ^ { 2 } } { h ^ { 4 } } s i n ^ { 4 } ( \frac { \beta } { 2 } ) + 1 = 1
$$

For the discriminant to be positive we require

$$
\frac { 4 \tau ^ { 2 } } { h ^ { 4 } } s i n ^ { 4 } ( \frac { \beta } { 2 } ) < 1
$$

As we want this inequality to be satisfied for all $\beta$ and $s i n ^ { 4 } \big ( \frac { \beta } { 2 } \big )$ is bounded between 0 and 1 we get the inequality.

$$
\frac { 4 \tau ^ { 2 } } { h ^ { 4 } } < 1
$$

and hence we require

$$
\tau < { \frac { h ^ { 2 } } { 2 } }
$$

# 2.6 Experimental Setup

To make a comparison between both methods we require the following quantity to be conserved

$$
N = \int _ { - \infty } ^ { \infty } | \psi | ^ { 2 } d x
$$

for a given accuracy $\epsilon$ for one soliton for a short time.

In light of this the following parameters were established: The amount of mesh points $N$ were fixed to be 512 with a simulation time $T = 1$ for both methods. The spacial length $L$ was chosen to be 64 for the Split Step method and 30 for the Finite Difference method. The time step $\tau$ for the Split Step method was chosen to be 0.01 (which also satisfied the stability requirement described in section 2.5). Simulating the solution for 8 time steps yielded $N \ = \ 0 . 0 0 0 9 7$ For the Finite Difference method the $\tau$ was chosen to be 0.001 running the simulation for 8 steps yielded $N \ = \ 0 . 0 0 0 3 9 .$ . These provided parameters satisfy an $\epsilon \ : = \ : 1 0 ^ { - 3 }$ . $N$ was computed using the composite trapezoidal rule.

All implementation was completed using Python 3 using Numpy for the matrix operations, Fourier Transform and Inverse Fourier Transform. Matplotlib was used to generate $2 D$ and $3 D$ plots. A shared code base was implemented in a file called helper.py which contained methods to create initial solutions, plotting functionality and to compute the quantity $N$ .Two additional files were generated which implemented the Split Step and Finite Difference method respectively. Jupyter Notebook was used to test both methods. All the source code can be found in the Appendix Section.

In the next section the results of multiple simulations using varying values for $S$ and velocities $v _ { 1 }$ and $v _ { 2 }$ are portrayed for both methods.

# 3 Results

# 3.1 Results with a small negative $S$ of $- 0 . 1$

For these experiments $S$ was set to a small negative value of $- 0 . 1$ .Figure 1 is the simulation of the Finite Difference method with both solitons initiated with a velocity of 20 in colliding directions. The same is portrayed in figure 2 however using the Split Step method. The first observation to be made is that both methods produce colliding solitons that superimpose on each other (the red spike) and then decompose back into their original states. A difference between the plots however is that the right soliton of the first figure seems to move left while the right soliton in the next figure seems to rather move in a straight path.

![](images/7e6a80a0f577a43182ee231d3149779c6fe5bd1cbc4aad314d226d8d3a4f3c26.jpg)  
Figure 1: Finite Difference: $S = - 0 . 1$ with $v _ { 1 } = 2 0$ and $v _ { 2 } = - 2 0$ (The right hand soliton is a bit tricky to view hence: The oval shows the initial position of the soliton and the arrows attached to the oval point the direction in which it is travelling. The top arrow shows where the soliton is after the collision; Behind the big red spike (superposition))

![](images/f2dcee1c0ef89cf71d36caa68c5e95a497cbc94b25fe4fecf911550408e54038.jpg)  
Figure 2: Split Step: $S = - 0 . 1$ with $v _ { 1 } = 2 0$ and $v _ { 2 } = - 2 0$

Other velocities were also tested using the Split Step method (Further experiments using the Finite Difference methods were omitted due to time constraints). In figure 3 both solitons were initiated with velocities of 10 in colliding directions. This change seems to be apparent from the visualisation as the solitons collide at a later time step during the simulation. Another experiment was carried out with the solitons initiated with a velocity of 1 in colliding directions which is depicted in figure 4. In this instance it can be observed that the solitons seems to be travelling parallel to each other as the velocities are not set large enough.

# 3.2 Results with a small positive $S$ of 0.4

Two experiments were run, one for the Finite Difference and the other for the Split Step, however this time with value of $S \ : = \ : 0 . 4$ .The velocities were set to 20 in colliding directions. More changes in velocities were not examined due to time constraints. Figure 5 portrays the result using the Finite Difference. It can be observed that change in $S$ has caused the solitons to raise in height and breadth. Figure 6 shows the same result however using the Split Step method from which the same observation can be made.

![](images/2bfab9c24665394ca753315161c9e85face2df087c67f6935bdc529dcf7ab5f1.jpg)  
Figure 3: Split Step: $S = - 0 . 1$ with $v _ { 1 } = 1 0$ and $v _ { 2 } = - 1 0$

![](images/f35a9d7e2c65d7941706ec9f8167a50c2286e997f0a272aa52798fe990d9264e.jpg)  
Figure 4: Split Step: $S = - 0 . 1$ with $v _ { 1 } = 1$ and $v _ { 2 } = - 1$

# 3.3 Results with a large negative $S$ of $- 1 0$

Previous configurations were kept however simulations were now run using a value of $- 1 0$ for $S$ and velocities of 10 and $- 1 0$ for the two solitons. A few observations can be made: Both initial solutions have a smaller soliton height; Both simulations vary largely from each other. Figure 7 results in a rough terrain with little structure where else figure 8 results in a somewhat more structured output where the solution superposition can be seen, however afterwards small artefact waves can be observed next to the main solitons which seem to dampen after the collision.

# 3.4 Results with a large positive $S$ of 2

Figure 9 and figure 10 represent simulations run with a value of $S = 2$ using the Finite Difference and Split Step method. The plots provided are $2 D$ as full simulations were not necessary as from these solutions throughout time it can be observed that the solutions become noisy and differ from each other. Increasing $S$ resulted in the the solitons increasing in height and decreasing in width.

![](images/6d649cba2a431289837a0c92d1067b36a9fe9a82471f7ba84654c249fa160449.jpg)  
Figure 5: Finite Difference: $S = 0 . 4$ with $v _ { 1 } = 2 0$ and $v _ { 2 } =$ $- 2 0$ The right hand solution is hidden behind the big red spike after the the collision. View Figure 1 to get a better idea.

![](images/7264b2d6548ceb2ad29e9d83b9e31285d4ddf4314a1bfd841c4f1cbe87418885.jpg)  
Figure 6: Split Step: $S = 0 . 4$ with $v _ { 1 } = 2 0$ and $v _ { 2 } = - 2 0$

# 3.5 Run Times

The time step of the Split Step method ran for 5 seconds and the Finite Difference for 2.5 seconds. However the Split Step had a total of 100 time steps compared to the 1000 times steps of the Finite Difference, hence the total running time using the Split Step was 500 seconds and the running time of the Finite Difference was 2500 seconds which is considerably longer.

# 4 Discussion

In regards to the numeric methods both produced similar plots for values of $S$ close to 0. A small positive value for $S$ raised the soliton heights where else a small negative value for $S$ reduced the soliton heights. Both schemes produced different results for larger values of $S$ (negative and positive values). This increase in the absolute value of $S$ increases the nonlinearity of the $P D E$ as it is a coefficient to the absolute value of $\psi$ which has a larger effect on the dynamics of the solutions hence the possible discrepancies between the solutions produced by the two schemes.

In regards to the numeric methods the Split Step was easier to implement and did not suffer from any instability or careful consideration to choose values $\tau$ and $h$ (in fact further experiments showed that $\tau$ could be increased even more without affecting the accuracy of the outcome). The Finite Difference on the other hand was harder to implement and to debug. In addition, stability analysis had to be performed to ensure that the method would actually converge. Finding parameters $\tau$ and $h$ was more difficult compared to the former method. Preliminary experiments demonstrated that even if the stability was satisfied, if $N$ was too small the solutions did not produce desired results. Thus $N$ had to be set large enough yet this forced the step size $\tau$ to be chosen to be considerably smaller than the one chosen for the Split Step method. In addition the spacial interval had to be reduced to satisfy the stability condition.

![](images/b9156c681e84b460a07e2edc0c843714d1dcdec53f3152366186f206e5d9245e.jpg)  
Figure 7: Finite Difference: $S = - 1 0$ with $v _ { 1 } ~ = ~ 1 0$ and $v _ { 2 } = - 1 0$

![](images/e6b203620c5f855caf679f757f5abe70a3f2ed36ae5280a70df9617c336230e2.jpg)  
Figure 8: Split Step: $S = - 1 0$ with $v _ { 1 } = 1 0$ and $v _ { 2 } = - 1 0$

The Split Step is superior to the Finite Difference in this regard as the spacial interval and time step can be relatively large with no stability issues. The Finite Difference method on the other hand required careful probing of the these parameters to avoid stability issues. Lastly the Finite Difference method requires more computational time than the Split Step method due to the small time step.

![](images/891339969185c545b1db5473c299fd44277a9731586fb0ed45e6b336328832f8.jpg)  
Figure 9: Finite Difference: $S = 2$ with $v _ { 1 } = 1 0$ and $v _ { 2 } =$ −10

![](images/a0bc73e157e6afd37386e6c0bac89e29ce273e91f5dce17dc3dc17e247040cb0.jpg)  
Figure 10: Split Step: $S = 2$ with $v _ { 1 } = 1 0$ and $v _ { 2 } = - 1 0$

Difference), is faster and can solve the equation on a larger spacial interval using a larger time step. Thus it is advisable to solve problems like these using a spectral method like the Split Step over a Finite Difference method. Investigation of the equation itself also showed that the parameter $S$ has a large effect on the dynamics of the solutions with smaller values producing solitons that collide with each other, superimpose and then decompose back into their original states whereas larger absolute values destroy the dynamics of the solutions.

# References

[Delfour et al., 1981] M. Delfour, M. Fortin, and G. Payr. Finite-difference solutions of a non-linear schrödinger equation. In Journal of computational physics, pages 277 288, 1981.   
[Kato, 1989] T. Kato. Nonlinear Schrodinger Equations. Springer, Berlin, Germany, 1989.   
[Keller and Isaacson, 1994] Herbert Keller and Eugene Isaacson. Analysis of numerical methods. Courier Corporation, 1994.   
[Serkin and Hasegawa, 2000] V.N. Serkin and A. Hasegawa. Novel soliton solutions of the nonlinear schrödinger equation model. In Physical Review Letters, page 4502, 2000.   
[Weideman and Herbst, 1986] J.A.C. Weideman and B.M. Herbst. Split-step methods for the solution of the nonlinear schrödinger equation. In SIAM Journal on Numerical Analysis, pages 485507, 1986.   
[Zemlyanaya and Alexeeva, 2011] E.V. Zemlyanaya and N.V. Alexeeva. Numerical study of time-periodic solitons in the damped-driven nls. In International Journal of Numerical Analysis and Modeling, pages 248261, 2011.

# 5 Conclusion

In this report the Schrödinger Equation with a particular nonlinearity was investigated using the Split Step and Finite Difference method. In practice it was found that the Split Step method does not suffer from stability issues (like the Finite

# 6 Appendix

# 6.1 Shared code base: helper.py

from sympy import \*   
from mpl_toolkits.mplot3d import Axes3D   
import matplotlib.pyplot as plt   
from matplotlib import cm   
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def initOneSoliton(off, v, L, N, S): psi $=$ np.zeros(N, dtype $: =$ np.complex_) $ { \mathrm { ~  ~ t ~ } } =  { \mathrm { ~  ~ 0 ~ } }$ $\mathrm { ~ h ~ } = \mathrm { ~ L / N ~ }$ $\begin{array} { r l } { \exists } & { { } = } \end{array}$ np.sqrt(2) $\mathrm { ~  ~ { ~ \cal ~ B ~ } ~ } = \mathrm { ~  ~ 3 ~ } / 2 - 2 \mathrm { ~  ~ \star ~ } \mathrm { ~  ~ { ~ \cal ~ S ~ } ~ }$ for i in range(N): $\begin{array} { c c } { \mathrm {  ~ x ~ } = \mathrm {  ~ i ~ } \mathrm {  ~ \star ~ } \mathrm {  ~ h ~ } ^ { - } - \mathrm {  ~ \ o f f ~ } } \\ { \mathrm {  ~ f ~ } = \mathrm {  ~ \sigma ~ } ( 2 \mathrm {  ~ \star ~ } \mathrm {  ~ a ~ } \star \mathrm {  ~ \ e x p ~ } ( \mathrm {  ~ a ~ } \star \mathrm {  ~ \times ~ } \mathrm {  ~ x ) ~ } ) } & { / } & { ( 1 \mathrm {  ~ + ~ } \mathrm {  ~ B ~ } \star \mathrm {  ~ \times ~ } ) } \\ { \mathrm { e x p ~ } ( 2 \mathrm {  ~ \star ~ } \mathrm {  ~ a ~ } \star \mathrm {  ~ \times ~ } \mathrm {  ~ x ~ } ) ) } \\ { \mathrm { p s i ~ } [ \mathrm {  ~ i ~ } ] \mathrm {  ~ \cdot ~ } = \mathrm {  ~ f ~ } \mathrm {  ~ \star ~ } \mathrm { e x p ~ } ( \mathrm { I ~ } \star \mathrm {  ~ t ~ } + \mathrm {  ~ I ~ } \star \mathrm {  ~ v ~ } \star \mathrm {  ~ \times ~ } ) } \end{array}$ return psi

def initTwoSoliton(xloff, vl, x2off, v2, L, N, S): psi = np.zeros(N, dtype $: =$ np.complex_) $ { \mathrm { ~  ~ t ~ } } =  { \mathrm { ~  ~ 0 ~ } }$ $\mathrm { ~ h ~ } = \mathrm { ~ L / N ~ }$ a = np.sqrt(2) $\mathrm { ~  ~ { ~ \cal ~ B ~ } ~ } = \mathrm { ~  ~ 3 ~ } / 2 - 2 \mathrm { ~  ~ \star ~ } \mathrm { ~  ~ { ~ \cal ~ S ~ } ~ }$ for i in range(N): $\begin{array} { r c l } { \times 1 } & { = } & { \bot \star \mathrm {  ~ \gamma ~ h ~ } - \mathrm {  ~ \ x 1 o f ~ f ~ } } \\ { \mathrm {  ~ \gamma ~ \in ~ ( } 2 \mathrm {  ~ \gamma ~ \star ~ } \mathrm {  ~ \gamma ~ \lrcorner ~ } \ \Leftrightarrow \ \exp { ( \mathrm {  ~ a ~ } \star \mathrm {  ~ \ x 1 ~ } ) } \ ) \mathrm {  ~ \beta ~ } / \mathrm {  ~ \gamma ~ } ( 1 \mathrm {  ~ \beta ~ + ~ } \mathrm {  ~ B ~ } \star \mathrm {  ~ \gamma ~ } ) } \\ & { \displaystyle \exp { ( 2 \mathrm {  ~ \gamma ~ \star ~ } \mathrm {  ~ \alpha ~ } \widehat { \otimes } \ \mathrm {  ~ \star ~ } \mathrm {  ~ \ x 1 ~ } ) } } \\ { \times 2 } & { = \mathrm {  ~ i ~ } \star \mathrm {  ~ \ h ~ } - \mathrm {  ~ \ x 2 o f ~ f ~ } } \\ { \mathrm {  ~ \gamma ~ \geq ~ ( } 2 \mathrm {  ~ \gamma ~ \star ~ } \mathrm {  ~ \alpha ~ } \star \mathrm {  ~ \ e x p ~ } ( \mathrm {  ~ a ~ } \star \mathrm {  ~ \ x 2 ~ } ) \ ) \mathrm {  ~ \beta ~ } / \mathrm {  ~ \gamma ~ } ( 1 \mathrm {  ~ \beta ~ + ~ } \mathrm {  ~ B ~ } \star \mathrm {  ~ \gamma ~ } ) } \\ & { \displaystyle \exp { ( 2 \mathrm {  ~ \gamma ~ \star ~ } \mathrm {  ~ \alpha ~ } \widehat { \otimes } \ \mathrm {  ~ \star ~ } \mathrm {  ~ \ x 2 ~ } ) } ) } \\ { \mathrm {  ~ \ p s i ~ } [ \mathrm {  ~ i ~ } ] } & { = } & { \mathrm {  ~ f ~ } 1 \mathrm {  ~ \star ~ } \ \exp { ( \mathrm {  ~ I ~ } \star \mathrm {  ~ \biggr ~ \cdot ~ } \mathrm {  ~ t ~ } + \mathrm {  ~ I ~ } \star \mathrm {  ~ \ v 1 ~ } \star \mathrm {  ~ x 1 ~ } ) } } \\ &  \displaystyle + \mathrm {  ~ f ~ } 2 \star \mathrm {  ~ \ e x p ~ } ( \mathrm {  ~ I ~ } \star  \end{array}$ return psi

def computeN(psil, psi2, L, N): $\begin{array} { r l } { \mathbb { N } \bot } & { { } = } \end{array}$ np.trapz(abs(psil), ${ \mathrm { d } } { \mathrm { x } } { = } \mathrm { L } / \mathrm { N }$ $\begin{array} { r l } { \mathbb { N } 2 } & { { } = } \end{array}$ np.trapz(abs(psi2), ${ \mathrm { d } } { \mathrm { x } } { = } \mathrm { L } / \mathrm { N } .$ return np.abs(N1 - N2)

def plot2d(y, L, N): plt.plot(np.linspace(0, L, num $\mathrm { _ { 1 = N } }$ ),abs(y)) plt.show()

def plot3D(psiEv, L, N, T, tau): fig $=$ plt.figure(figsize $=$ (20,10)) ax $=$ fig.gca(projection $= ^ { 1 }$ 3d')

# Make data. $\begin{array} { r l } { \mathrm { ~ X ~ } } & { { } = } \end{array}$ np.arange(0, L, L/N) $\begin{array} { r l } { \mathrm { Y } } & { { } = } \end{array}$ np.arange(0, T, tau) X, Y $=$ np.meshgrid(X, Y)

# Plot the surface.   
surf $=$ ax.plot_surface(X, Y, psiEv, cmap ${ \bf \Pi } = { \bf C } { \bf m }$ .jet, linewidth ${ \tt = } 1 0$ , antialiased $\underline { { \underline { { \mathbf { \Pi } } } } } =$ True, rstride $^ { = 1 }$ , cstride $^ { = 1 }$ )

# # Customize the z axis.

ax.set_zlim(0, 2)

ax.zaxis.set_major_locator(LinearLocator(10)) ax.zaxis.set_major_formatter(FormatStrFormatter('%.1 #ax.view_init(30, 190)

# Add a color bar which maps values to colors. fig.colorbar(surf, shrink $= 0 \cdot 5$ , aspect ${ } = \mathtt { 1 0 }$ )

plt.show()

# 6.2 SplitStep.py

import numpy as np from sympy import $\star$ import helper as h

$\mathrm { ~ T ~ } = \mathrm { ~ 1 ~ }$   
$\mathtt { t a u \ = \ 0 . 0 1 }$   
$L ~ = ~ 6 4$   
$\mathrm { ~ \tt ~ N ~ } = \mathrm { ~ \tt ~ 5 1 2 ~ }$   
$\mathrm { ~  ~ { ~ S ~ } ~ } = \mathrm { ~  ~ - 0 ~ } . 1$ #3/4 - 0.1

def splitstep(psi): # Nonlinear Part for i in range(N): $\mathbf { \Sigma } _ { \subset } =$ (abs(psi[i]) $\star$ abs(psi[i])) coef $=$ $\mathrm { ~  ~ { ~ \ b ~ = ~ } ~ } ( \ I \mathrm { ~ \bf ~ \star ~ } \mathrm { ~ \bf ~ C ~ } )$ $( { \mathrm { ~ 1 ~ \ t ~ } } \ { \mathrm { ~ S ~ } } \ { \star } \ { \mathrm { ~ C ~ } } )$ psi[i] $=$ psi[i] $\star$ exp(coef $\star$ tau) # Fourier Transform $\begin{array} { r l } { \mathbf { \Sigma } } & { { } \subset \mathbf { \Sigma } } \end{array} =$ np.fft.fftshift(np.fft.fft(psi)) # Move in Fourier Space #n $=$ np.linspace(-N/2, $\mathbb { N } / 2 + 1$ , num $\mathrm { _ { 1 = N } }$ ) for i in range(N): e = (-2 \* i \* i \* pi $\star$ pi) / (L \* L); $\mathbf { C } \left[ \begin{array} { l } { \mathbf { i } } \end{array} \right] \ \mathbf { \Sigma } = \mathbf { \Sigma }$ exp(tau \* I \* e) \* c[i]; # Convert back to physical space psi $=$ np.fft.ifft(np.fft.fftshift(c)); return psi

# Run the Simulation   
psi $=$ h.initTwoSoliton(8, 20, 18, -20, L, N, S)   
psiEv $=$ np.zeros(shape $=$ (int(T/tau), int(N)))   
for i in range(int(T/tau)): psiEv[i] $=$ abs(psi) h.plot2D(psi, L, N) # Plot 2D graph at every time step psi $=$ splitstep(psi)   
h.plot3D(psiEv, L, N, T, tau) # Plot 3D graph   
def computeN(): psil $=$ h.initOneSoliton(8, 10, L, N, S) ps $\begin{array} { l l } { { \mathrm { i } 2 \ = \ } } \end{array} .$ psil for i in range(8): h.plot2D(psi2, L, N) psi2 $=$ splitstep(psi2) print(h.computeN(psil, psi2, L, N))

# 6.3 FiniteDifference.py

import numpy as np from sympy import \*

$\mathrm { ~ T ~ } = \mathrm { ~ 1 ~ }$   
tau $\mathbf { \varepsilon } = \mathbf { \varepsilon } _ { 0 } . 0 0 1$   
$\mathrm { ~ \tt ~ { ~ L ~ } ~ } = \mathrm { ~ \tt ~ { ~ 4 ~ 0 ~ } ~ }$   
$\mathrm { ~ \tt ~ N ~ } = \mathrm { ~ \tt ~ 5 1 2 ~ }$   
$\mathrm { ~  ~ { ~ S ~ } ~ } = \mathrm { ~  ~ - 0 ~ } . 1$ #3/4 - 0.1

def FD(phi1): newPhi $=$ np.zeros(N, dtype $: =$ np.complex_) $\mathrm { ~ h ~ } = \mathrm { ~ L / N ~ }$

for i in range(0, N): if $\begin{array} { l l l } { \dot { \mathrm { ~ \scriptsize ~ 1 ~ } } } & { = = } & { 0 } \end{array}$ : phixx $=$ (phil[N - 1] − 2 \* phil[i] $+ \mathrm {  ~ \ p h i 1 \left[ \dot { 1 } \dot { \mathrm {  \Phi } } + \frac { } { } 1 \right] } \mathrm {  \langle ~ \langle ~ \langle ~ \hbar ~ \star ~ \hbar \rangle ~ }$ elif $\\begin{array} { r } { \mathrm { ~  ~ \underline { ~ } { ~ } ~ } = = \mathrm { ~  ~ N - 1 ~ } } \end{array}$ : phixx $=$ (phil[i - 1] - 2 $\star$ phil[i] $+ \texttt { p h i l } [ 0 ] \ \rangle \ / \ ( \texttt { h } \star \texttt { h } )$ else: phixx $=$ (phil[i - 1] - 2 $\star$ phil[i] $^ +$ $- \mathrm {  ~ \ p h i 1 ~ } [ \mathrm {  ~ i ~ } + \mathrm {  ~ 1 ~ } ]$ ) / (h \*h) psiSquared $=$ np.abs(phil[i]) $\star$ np.abs(phil[i]) $\begin{array} { r l } { \mathrm { ~ \mathbb ~ { ~ A ~ } ~ } = } \end{array}$ psiSquared / $\begin{array}{c} (  { { \mathchoice { \mathrm { 1 } } { \mathrm { 1 } } } { \mathrm { 1 } } } +  { { \mathchoice { \mathrm { S } } { \mathrm {  ~ \star ~ } } { \mathrm { 1 } } { \mathrm { ~  ~ \star ~ } } { \mathrm { ~  ~ \times ~ } } } } \end{array} )$ psiSquared) newPhi[i] = I $\star$ tau $\star$ ((1/2) $\star$ phixx - A $\star$ phil[i]) $^ +$ phil[i]

# return newPhi

def CD(phil, phi2): newPhi $=$ np.zeros(N, dtype $=$ np.complex_) $\mathrm { ~ h ~ } = \mathrm { ~ L / N ~ }$ for i in range(0, N): if $\mathrm { ~  ~ { ~ i ~ } ~ } = = \mathrm { ~  ~ { ~ 0 ~ } ~ }$ : phixx $=$ (phil[N − 1] - 2 \* phil[i] $+ \mathrm {  ~ \ p h i 1 \left[ \dot { 1 } \dot { \mathrm {  \Phi } } + \frac { } { } 1 \right] } \mathrm {  \langle ~ \langle ~ \langle ~ \hbar ~ \star ~ \hbar \rangle ~ }$ elif( $\\\\\\mathrm { ~ i ~ \ = = ~ \ N - 1 }$ : phixx $=$ (phil[i − 1] - 2 \* phil[i] $+ \texttt { p h i l } [ 0 ] \ \rangle \ / \ ( \texttt { h } \star \texttt { h } )$ else: phixx $=$ (phil[i − 1] - 2 $\star$ phil[i] $^ +$ phil[i + 1]) / $( \mathrm { ~ h ~ \ } \star \ \mathrm { ~ h ~ } )$ psiSquared $=$ np.abs(phil[i]) $\star$ np.abs(phil[i]) $\begin{array} { r l } { \mathrm { ~ \mathbb ~ { ~ A ~ } ~ } = } \end{array}$ psiSquared / $ { \mathrm { ~  ~ { ~ \cal ~ 1 ~ } ~ } } +  { \mathrm { ~  ~ { ~ \cal ~ S ~ } ~ } } \star$ psiSquared) newPhi[i] $= \begin{array} { l l l l } { { 2 } } & { { \star } } & { { { \underline { { { \sf  ~ \cal ~ T ~ } } } } } } \end{array} ,$ \* tau $\star$ ((1/2) \* phixx $^ +$ A $\star$ phil[i]) $^ +$ phi2[i] return newPhi

# Run the Simulation   
firstPsi $= \mathrm { ~ h ~ }$ .initTwoSoliton(10, 20, 20, -20, L, N, S)   
secondPsi $=$ FD(firstPsi)   
psiEv $=$ np.zeros(shape $=$ (int(T/tau), int(N)))   
for i in range(int(T/tau)): psiEv[i] $=$ abs(firstPsi) h.plot2D(firstPsi, L, N) # Plot 2D graph at every time step temp $=$ firstPsi firstPsi $=$ secondPsi secondPsi $=$ CD(secondPsi, temp)   
h.plot3D(psiEv, L, N, T, tau) # Plot 3D graph   
firstPsi $= \mathrm { ~ h ~ }$ .initOneSoliton(8, 20, L, N, S)   
secondPsi $=$ FD(firstPsi)   
psil $=$ firstPsi   
for i in range(8): temp $=$ firstPsi firstPsi $=$ secondPsi secondPsi $=$ CD(secondPsi, temp)   
print(h.computeN(psil, secondPsi, L, N))