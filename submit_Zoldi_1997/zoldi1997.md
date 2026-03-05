# Parallel Implementations of the Split-Step Fourier Method for Solving Nonlinear Schrödinger Systems

S.M. Zoldi

zoldi@phy.duke.edu Department of Physics, Duke University Box 90305, Durham, NC 27708-0305

V. Ruban, A. Zenchuk ruban@itp.ac.ru, zenchuk@itp.ac.ru L.D. Landau Institute for Theoretical Physics Russian Academy of Sciences, 2 Kosygin St., Moscow, Russia

# S. Burtsev

burtsev@cnls.lanl.gov Theoretical Division and Center for Nonlinear Studies Los Alamos National Laboratory, Los Alamos, NM 87545

November 26, 2024

# Abstract

We present a parallel version of the well-known Split-Step Fourier method (SSF) for solving the Nonlinear Schrödinger equation, a mathematical model describing wave packet propagation in fiber optic lines. The algorithm is implemented under both distributed and shared memory programming paradigms on the Silicon Graphics/Cray Research Origin 200. The 1D Fast-Fourier Transform (FFT) is parallelized by writing the 1D FFT as a 2D matrix and performing independent 1D sequential FFTs on the rows and columns of this matrix. We can attain almost perfect speedup in SSF for small numbers of processors depending on both problem size and communication contention. The parallel algorithm is applicable to other computational problems constrained by the speed of the 1D FFT.

# 1 Introduction

The Nonlinear Schrödinger equation (NLSE)

$$
i A _ { t } + \sigma d \frac { 1 } { 2 } A _ { x x } + A ^ { * } A ^ { 2 } = G ,
$$

is a nonlinear partial differential equation that describes wave packet propagation in a medium with cubic nonlinearity. Technologically, the most important application of NLSE is in the field of nonlinear fiber optics [1, 2]. The parameter $\sigma$ specifies the fiber anomalous group velocity dispersion ( $\sigma = 1$ ) or normal group velocity dispersion $\sigma = - 1$ ), while the parameter $d$ defines the normalized absolute value of the fiber's dispersion. The perturbation $G$ is specified by details of the physical fiber being studied.

In the special case $G = 0$ , NLSE is integrable [3] and can be solved analytically. In general if $G \neq 0$ NLSE must be solved numerically. One of the most popular numerical methods to solve the perturbed NLSE is the Split-Step Fourier method (SSF) [2]. For small-scale calculations, serial implementations of SSF are adequate; however, as one includes more physics in the simulation, the need for large numbers of Fourier modes to accurately solve NLSE equation demands parallel implementations of SSF.

Many fiber optics problems demand large-scale numerical simulations based on the SSF method. One class of such problems involves accurate modeling of wave-length division multiplexed (WDM) transmission systems where many optical channels operate at their own frequencies and share the same optical fiber. WDM is technologically important as it is one of the most effective ways to increase the transmission capacity of optical lines [1, 2, 4]. To accurately model WDM one needs to include large numbers of Fourier harmonics in the numerical simulation to cover the entire transmission frequency band. Moreover, in WDM systems different channel pulses propagate at different velocities and, as a result, collide with each other. At the pulse collision, Stokes and anti-Stokes sidebands are generated; these high frequency perturbations lead to signal deterioration [5, 6]. Another fundamental nonlinear effect called four-wave mixing (FWM) [7] must be accurately simulated as the FWM components broaden the frequency domain which requires even larger numbers of Fourier modes for accurate numerical simulation.

To suppress the FWM [5, 6] and make possible the practical realization of WDM, one can use dispersion management (concatenation of fiber links with variable dispersion characteristics). The dispersion coefficient $d$ in NLSE is now no longer constant but represents a rapidly varying piecewise constant function of the distance down the fiber. As a result, one must take a small step size along the fiber to resolve dispersion variations and the corresponding pulse dynamics. A final reason to include a large number of Fourier modes in numerical simulations is to model the propagation of pseudorandom data streams over large distances.

All of the above factors make simulation of NLSE quite CPU intensive. Serial versions of the split-step Fourier method in the above cases may too be slow even on the fastest modern workstations. To address the issue of accurately simulating physical telecommunication fibers in a reasonable amount of time, we discuss the parallelization of SSF algorithm for solving NLSE. Our parallel SSF algorithm is broadly applicable to many systems and not limited to the solution of NLSE. We consider an algorithm appropriate for multiprocessor workstations. Parallel computing on multiprocessor systems raises complex issues including solving problems efficiently with small numbers of processors, limitations due to the increasingly complex memory hierarchy, and the communication characteristics of shared and distributed multiprocessor systems.

Current multiprocessors have evolved towards a generic parallel machine, which shares characteristics of both shared and distributed memory computers. Therefore most commercial multiprocessors support both shared memory and distributed memory programming paradigms. The shared memory paradigm consists of all processors being able to access some amount of shared data during the program execution. This addressing of memory on different nodes in shared memory multiprocessors causes complications in writing efficient code. Some of the most destructive complications are: cache hierarchy inefficiency (alignment and data locality), false sharing of data contained in a cache block, and cache thrashing due to true sharing of data. Most vendors provide compiler directives to share data and divide up computation (typically in the form of loop parallelism) which in conjunction with synchronization directives can be used to speed up many sequential codes. In distributed memory programming, each processor works on a piece of the computation independently and must communicate the results of the computation to the other processors. This communication must be written explicitly into the parallel code, thus requiring more costly development and debugging time. The communication is typically handled by libraries such as the message passing interface (MPI) [8] which communicates data through Ethernet channels or through the existing memory system. Our primary goal is to present a parallel split-step Fourier algorithm and implement it under these two different parallel programming paradigms on a 4-processor Silicon Graphics/Cray Research Origin 200 multiprocessor computer.

The remainder of this paper is organized as follows. In Section 2, we recall a few basics of the the split-step Fourier method. In Section 3, we introduce the parallel algorithm for SSF. Timing results and conclusions are given in Sections 4 and 5, respectively.

# 2 Split-Step Fourier Method

The Split-Step Fourier (SSF) method is commonly used to integrate many types of nonlinear partial differential equations. In simulating Nonlinear Schrödinger systems (NLS) SSF is predominantly used, rather than finite differences, as SSF is often more efficient [9]. We remind the reader of the general structure of the numerical algorithm [2].

NLS can be written in the form:

$$
\frac { \partial A } { \partial t } = ( L + N ) A ,
$$

where $L$ and $N$ are linear and nonlinear parts of the equation. The solution over a short time interval $\tau$ can be written in the form

$$
A ( t + \tau , x ) = \exp ( \tau L ) \exp ( \tau N ) A ( t , x )
$$

vhere the linear operator in NLS acting on a spatial field $B ( t , x )$ is written in Fourier ipace as,

$$
\exp ( \tau L ) B ( t , x ) = F ^ { - 1 } \exp ( - i k ^ { 2 } \tau ) F B ( t , x )
$$

where $F$ denotes Fourier transform (FT), $F ^ { - 1 }$ denotes the inverse Fourier transform, and $k$ is the spatial frequency.

We split the computation of $A$ over time interval $\tau$ into 4 steps:

Step 1. Nonlinear step: Compute $A _ { 1 } = \exp ( \tau N ) A ( t , x )$ (by finite differences);   
Step 2. Forward FT: Perform the forward FFT on $A _ { 1 }$ : $A _ { 2 } = F A _ { 1 }$ ;   
Step 3. Linear step: Compute $A _ { 3 } = \exp ( \tau L ) A _ { 2 }$ ;   
Step 4. Backward FT: Perform the backward FFT on $A _ { 3 }$ : $A ( t + \tau ) = F ^ { - 1 } A _ { 3 }$ ;

To discretize the numerical approximation of the above algorithm, the potential $A$ is discretized in the form: $A _ { l } = A ( l h ) ; l = 0 , \dots , N - 1$ , where $h$ is the space-step and $N$ is the total number of spatial mesh points.

The above algorithm of the Split-Step Fourier (SSF) method is the same for both sequential and parallel code. Parallel implementation of this algorithm involves parallelizing each of the above four steps.

# 3 The Parallel Version of the Split-Step Fourier (SSF) Method

By distributing computational work between several processors, one can often speed up many types of numerical simulations. A major prerequisite in parallel numerical algorithms is that sufficient independent computation be identified for each processor and that only small amounts of data are communicated between periods of independent computation. This can often be done trivially through loop-level parallelism (shared memory implementations) or non-trivially by storing true independent data in each processor's local memory. For example, the nonlinear transformation in the SSF algorithm involves the independent computation over subarrays of spatial elements of $A ( l )$ . Therefore $P$ processors each will work on sub-arrays of the field $A$ , e.g., the first processor updates $A _ { 0 }$ to $A _ { ( N / P - 1 ) }$ , the second processor updates $A _ { N / p }$ to $A _ { 2 ( N / p ) - 1 }$ , etc.

In the 1D-FFT, elements of $( F A ) _ { k }$ can not be computed in a straightforward parallel manner, because all elements $A _ { l }$ are used to construct any element of $( F A ) _ { k }$ . The problem of 1D FFT parallelization has been of great interest for vector [10, 11] and distributed memory computers [12]. These algorithms are highly architecture dependent, involving efficient methods to do the data re-arrangement and transposition phases of the 1D FFT. Communication issues are paramount in 1D FFT parallelization and in the past have exploited classic butterfly communication patterns to lessen communication costs [12]. However, due to a rapid change in parallel architectures, towards multiprocessor systems with highly complex memory hierarchies and communication characteristics, these algorithms are not directly applicable to many current multiprocessor systems. Shared memory multiprocessors often have eficient communication speeds, and we therefore implement the parallel 1D FFT by writing $A _ { l }$ as a two dimensional array, in which we can identify independent serial 1D FFTs of rows and columns of this matrix. The rows and columns of the matrix $A$ can be distributed to divide up the computation among several processors. Due to efficient communication speeds, independent computation stages, and the lack of the transposition stage of the 1D FFT in SSF computations, we show that this method exploits enough independent computation to result in a significant speedup using a small number of processors.

# 3.1 Algorithm of Parallel SSF

The difficulty parallelizing the split-step Fourier algorithm is in steps 2 and 4, as the other two steps can be trivially evolved due to the natural independence of the data $A$ and $A _ { 2 }$ . In Step 2 and Step 4 there are non-trivial data dependences over the entire range $0 ~ < = ~ l ~ < = ~ N$ of $A _ { 1 } ( l )$ and $A _ { 3 } ( l )$ which involve forward and backward Fourier Transforms (FFT and BFT). The discrete forward Fast Fourier Transform (FFT) is of the form

$$
F ( k ) = \sum _ { l = 0 } ^ { N - 1 } A ( l ) \exp \left( - { \frac { 2 \pi i l k } { N } } \right)
$$

which requires all elements of $A ( l )$ . Several researchers have suggested parallel 1D Fast Fourier Transform algorithms [10, 11, 12], but to date there exist no vendor-optimized parallel 1D FFT algorithms. Therefore implementations of these algorithms are highly architecture dependent. Parallel 1D FFT algorithms must deal with serious memory hierarchy and communication issues in order to achieve good speedup. This may be the reason why vendors have been slow to support the computational community with fast parallel 1D FFT algorithms. However, we show below that we can get significant parallel speedup due to the elimination of the transposition stage in 1D FFT for SSF methods and due to exploitation of independent computation by performing many sequential 1D FFTs on small subarrays of $A ( l )$ .

Our method of parallelizing the SSF algorithm requires dividing the 1D array $A ( l )$ into subarrays which are processed using vendor optimized sequential 1D FFT routines. We assume the 1D array $A$ is of the dimension $N$ of the product of two integer numbers: $N = M _ { 0 } \times M _ { 1 }$ . Therefore $A$ can be written as a 2D matrix of size $M _ { 0 } \times M _ { 1 }$ . As a result, we can reduce the expression for the Fourier transform of the array $A$ to the form

$$
F ( M _ { 1 } k _ { 1 } + k _ { 0 } ) = \sum _ { n _ { 0 } = 0 } ^ { M _ { 0 } - 1 } \sum _ { n _ { 1 } = 0 } ^ { M _ { 1 } - 1 } A ( M _ { 0 } n _ { 1 } + n _ { 0 } ) e x p \left( { \frac { - 2 \pi i } { M _ { 0 } M _ { 1 } } } ( M _ { 0 } n _ { 1 } + n _ { 0 } ) ( M _ { 1 } k _ { 1 } + k _ { 0 } ) \right)
$$

$$
= \sum _ { n _ { 0 } = 0 } ^ { M _ { 0 } - 1 } f ( k _ { 0 } , n _ { 0 } ) e x p \left( { \frac { - 2 \pi i } { M _ { 0 } M _ { 1 } } } n _ { 0 } k _ { 0 } \right) e x p \left( { \frac { - 2 \pi i } { M _ { 0 } } } n _ { 0 } k _ { 1 } \right)
$$

where $F$ is the Fourier transform of $A$ and $f$ is the result of $M _ { 1 }$ -size Fourier transform of $A ( M _ { 0 } n _ { 1 } + n _ { 0 } )$ with fixed $n _ { 0 }$

$$
f ( k _ { 0 } , n _ { 0 } ) = \sum _ { n _ { 1 } = 0 } ^ { M _ { 1 } - 1 } A ( M _ { 0 } n _ { 1 } + n _ { 0 } ) e x p \left( { \frac { - 2 \pi i } { M _ { 1 } } } n _ { 1 } k _ { 0 } \right)
$$

$$
n _ { 0 } , k _ { 1 } = 0 , . . . , M _ { 0 } - 1 \qquad n _ { 1 } , k _ { 0 } = 0 , . . . , M _ { 1 } - 1 .
$$

The reduced expression Eq. (2) demonstrates that the $N = M _ { 0 } * M _ { 1 }$ Fourier transform $F$ is obtained by making $M _ { 0 }$ size Fourier transforms of $\begin{array} { r } { f ( k _ { 0 } , n _ { 0 } ) e x p \left( \frac { - 2 \pi i } { N } n _ { 0 } k _ { 0 } \right) } \end{array}$ for fixed $k _ { 0 }$ . Therefore the 1D array $A$ is written as a 2D matrix of size $M _ { 0 } \times M _ { 1 }$ with elements $a _ { j k }$   
$( A ( 0 ) , . . . , A ( M _ { 0 } - 1 ) )$ in the first column, $( A ( M _ { 0 } ) , . . . , A ( 2 M _ { 0 } - 1 ) )$ in the second column, etc. We use this matrix $a _ { j k }$ in our parallel FFT-algorithm:

Step 1. Independent $M _ { 1 }$ -size FFTs on rows of $a _ { j k }$ . Step 2. Multiply elements $a ( j , k )$ by a factor $E _ { j k } = e x p ( - ( 2 \pi i / N ) \cdot j \cdot k )$ Step 3. Independent $M _ { 0 }$ -size FFTs on columns of $a _ { j k }$ .

The result of Step 1 - Step 3 is the $N = M _ { 0 } * M _ { 1 }$ 1D Fourier transform of $A$ stored in rows: $( F ( 0 ) , . . . , F ( M _ { 1 } - 1 ) )$ in the first row, $( F ( M _ { 1 } ) , . . . , F ( 2 M _ { 1 } - 1 ) )$ in the second row, and so on. To regain the proper ordering of $A$ (how elements were originally stored in matrix $a _ { j k }$ ) requires a transposition of the matrix which is the last step in a parallel FFT algorithm.

In the SSF method, the transposition is not necessary as we apply a linear operator $L ( k )$ and then take the steps: Step 1 - Step 3 in reverse order. This avoids the transposition because one can define a transposed linear operator array and multiply $a _ { j k }$ by this transposed linear operator. Then Step 1 - Step 3 are performed in reverse order with the conjugate of the exponential term in Step 2.

The complete SSF parallel algorithm consists of the following steps:

Step 1. nonlinear step   
Step 2. row-FFT   
Step 3. multiply by $E$   
Step 4. column-FFT   
Step 5. linear step (transposed linear operator)   
Step 6. column-BFT   
Step 7. multiply by $E ^ { * }$ (the complex conjugate of $E$ )   
Step 8. row-BFT

The parallelization is due to the natural independence of operations in steps 1, 3, 5, and 7 and by the row and column subarray FFTs in steps 2, 4, 6, and 8. The row and column subarray FFTs of size $M _ { 1 }$ and $M _ { 0 }$ are performed independently with serial optimized 1D FFT routines. Working with subarray data, many processors can be used to divide up the computation work resulting in significant speedup if communication between processors is efficient. Further, smaller subarrays allows for better data locality in the primary and secondary caches. The implementation details of the shared-memory and the distributed memory parallel SSF algorithm outlined above depend on writing Steps 1 - 8 using either shared memory directives or distributed memory communication library calls (MPI).

# 3.2 Shared Memory Approach

Much of the SSF parallel algorithm outlined above can be implemented with "\$doacross" directives to distribute independent loop iterations over many processors. The FFTs of size $M _ { 0 }$ and $M _ { 1 }$ are implemented by distributing the 1D subarray FFTs of rows and columns over the P available processors. The performance can be improved drastically by keeping the same rows and columns local in a processor's secondary cache to alleviate true sharing of data from dynamic assignments of sub-array FFTs by the "\$doacross" directive. The subarray FFTs are performed using vendor optimized sequential 1D FFT routines which are designed specifically for the architecture.

It is efficient to perform all column operations (Steps 3 - 7) in one pass: copying a column into local sub-array $S$ contained in the processor's cache and in order, multiply by the exponents in Step 3, perform the $M _ { 0 }$ -size FFT of $S$ , multiply by the transposed linear operator $\exp ( \tau L )$ , invert the $M _ { 0 }$ -size FFT, multiply by the conjugate exponents, and finally store $S$ back into the same column of $a$ . This allows for efficient use of the cache, reducing false/true sharing as we perform many operations on each subarray.

# 3.3 Distributed Memory Approach

The Massage Passing Interface (MPI) is a tool for distributed parallel computing which has become a standard used on a variety of high-end parallel computers to weakly coupled distributed networks of workstations (NOW) [8]. In distributed parallel programming, different processors work on completely independent data and explicitly use send and receive library calls to communicate data between processors.

To implement the distributed parallel SSF algorithm for the Nonlinear Schrödinger system (NLS), one needs to distribute the rows of array $A$ among all P available processors. Then Steps 1 - 3 can be executed without communication between processors. After these steps, it is necessary to endure the communication cost of redistributing the elements of $A$ among the P processors. Each processor must send a fraction $\textstyle { \frac { 1 } { P } }$ of its data to each of the other processors. Then each processor will have the correct data for Steps 4 -7 and column operations are performed independently on all P processors. Finally, there is a second redistribution prior to Step 8. To make $T$ steps of the SSF algorithm, we use the following scheme:

# subroutine Distributed SSF

distribute rows among processors Step 1. nonlinear step Step 2. row-FFT Step 3. multiply by a factor $E$

for $i = 1$ to $T - 1$ do data redistribution Step 4. column-FFT Step 5. linear step Step 6. column-BFT Step 7. multiply by a factor $E ^ { * }$ data redistribution Step 8. row-BFT Step 1. nonlinear step Step 2. row-FFT Step 3. multiply by a factor $E$   
end do   
Step 4. column-FFT   
Step 5. linear step   
Step 6. column-BFT   
Step 7. multiply by a factor $E ^ { * }$   
Step 8. row-BFT

# end

The large performance cost in this algorithm is the redistribution of data between row and column operations. If the row and column computational stages result in significant speedup compared to the communication expense of redistributing the matrix data, then this algorithm will be successful. This depends crucially on fast communication between processors which is usually the case for shared memory multiprocessors and less so for NOW computers.

# 4 Results

We performed timings of the parallel SSF algorithm on the Silicon Graphics/Cray Research Origin 200 multiprocessor. The Origin 200 was used because it allows for both shared and distributed memory parallel programming and models a generic multiprocessor. The Origin 200 is efficient at fine-grained parallelism which typically makes shared memory programming both efficient and easy. The Origin 200 workstation used in this study consisted of four MIPS R10000 64-bit processors (chip revision 2.6) with MIPS R10010 (revision 0.0) foating point units running at 180MZ. The primary cache consisted of a 32KB 2-way set-associative instruction cache and a 32KB 2-way set-associative data cache. Each processor also had a 1MB 2-way set-associative secondary cache. The machine had a sustained 1.26GB/sec memory bandwidth and 256MB of RAM.

The operating system was IRIX 6.4. We used a Mongoose f77 version 7.10 Fortran compiler. For the parallel programming we used MPI version 1.1 for the distributed computing and the native "\$doacross" and synchronization directives provided by the f77 compiler for shared memory programming. All timings are of the total wall-clock time for the code to both initialize and execute.

# 4.1 Timings

For the following timings, we use $M _ { 0 } = M _ { 1 } = 2 ^ { K }$ , so that the entire 1D array is of size $N = 2 ^ { 2 K }$ . The one-processor implementation of parallel SSF was 10% to $2 0 \%$ faster than serial SSF code using vendor optimized 1D FFTs of the entire array of size $N =$ $2 ^ { 2 K }$ . This improvement is due to better cache coherence using smaller subarrays, as an entire subarray can be contained in the L1 cache and is due to the fact that the single processor parallel SSF does not do the transposition stage of the 1D FFT. All timings are compared to the one-processor parallel code at the same optimization level (compared to sequential SSF the below speedups are even more impressive). For shared memory parallel implementations, we find over the range of $2 ^ { 1 2 } < N < 2 ^ { 1 8 }$ that two node SSF implementations have good speedup (SU) with a maximum speedup at $N = 2 ^ { 1 6 }$ . Using four nodes, for small array sizes we have $1 / 4$ less work per processor, but more contention due to the sharing of pieces of data contained in the secondary caches of four different processors. At $N = 2 ^ { 1 6 }$ , we again see the maximum speedup (now for 4 nodes), reflecting that the ratio of computational speed gain to communication contention is optimal at this problem size.

# Shared Memory

<table><tr><td rowspan=8 colspan=1>array size (N)number of steps (T )T1pr. (sec)T2pr (sec)T4pr (sec)SU = T2pr/T1pr.SU = T4pr/T1pr.</td><td></td><td></td><td></td><td></td></tr><tr><td rowspan=1 colspan=1>N = 212</td><td rowspan=1 colspan=1>N = 214</td><td rowspan=1 colspan=1>N = 216</td><td rowspan=1 colspan=1>N = 218</td></tr><tr><td rowspan=1 colspan=1>T = 8000</td><td rowspan=1 colspan=1>T = 2000</td><td rowspan=1 colspan=1>T = 500</td><td rowspan=1 colspan=1>T = 125</td></tr><tr><td rowspan=1 colspan=1>49.5</td><td rowspan=1 colspan=1>51.5</td><td rowspan=1 colspan=1>65.5</td><td rowspan=1 colspan=1>97.5</td></tr><tr><td rowspan=1 colspan=1>29.5</td><td rowspan=1 colspan=1>30.5</td><td rowspan=1 colspan=1>33.5</td><td rowspan=1 colspan=1>61.0</td></tr><tr><td rowspan=1 colspan=1>19.5</td><td rowspan=1 colspan=1>18.5</td><td rowspan=1 colspan=1>19.5</td><td rowspan=1 colspan=1>34.5</td></tr><tr><td rowspan=1 colspan=1>1.7</td><td rowspan=1 colspan=1>1.7</td><td rowspan=1 colspan=1>2.0</td><td rowspan=1 colspan=1>1.6</td></tr><tr><td rowspan=1 colspan=1>2.5</td><td rowspan=1 colspan=1>2.8</td><td rowspan=1 colspan=1>3.4</td><td rowspan=1 colspan=1>2.8</td></tr></table>

Under the shared memory programming model, subarrays are continually distributed among processors to divide up the computational work. Data in a single subarray may be contained on one or more processors requiring constant communication. The data contained in each processor's L2 cache is of size O(N/P), where P is the number of processors. Contention in the memory system is modeled as being proportional to $O ( ( N / P ) ^ { 2 } )$ which reflects the communication congestion for sharing data of large working sets. Further unlike the serial code, the parallel code endures a communication time to send data between processors proportional to ${ \cal O } ( N / P ) \tau _ { c }$ , where $\tau _ { c }$ is the time to transfer a L2 cache block between processors. Finally, the time to perform the 1D FFT is approximately $N L o g ( N ) \tau _ { M }$ , where $^ { \tau _ { M } }$ is the time to perform a floating point operation. A simple formula for the speedup (SU) of the shared memory FFT is

$$
S U = \frac { \tau _ { M } N L o g ( N ) } { \tau _ { M } N L o g ( \sqrt { N } ) / P + \tau _ { c } N / P + f ( N / P ) ^ { 2 } ) } ,
$$

where $f$ is a small number reflecting contention in the communication system. If $N = 2 ^ { K }$ we can simplify the above expression,

$$
S U = \frac { 2 P } { ( 1 + \xi / K + f 2 ^ { K } / ( P K ) ) } ,
$$

where the constants are absorbed into $f$ and $\xi$ . With $f = 0$ (no contention) one predicts for fixed P that the speedup increases for larger and larger problem size $N$ . However, for $f \neq 0$ the speedup eventually decreases with larger $N$ due to contention of communicating small pieces of subarray data between arbitrary processors. This equation reflects the trend seen in our empirical data of speedup for shared memory SSF, where speedup attains a maximum with problem size at $N = 2 ^ { 1 6 }$ .

The above SU formula must be reinterpreted for distributed SSF due to the implicit independent computational stages where no data is communicated between processors unlike shared memory SSF. Distributed SSF uses communication stages to send data between processors and does not involve contention due to sharing data between P processors during computation stages. Distributed MPI timings are compared to a one-processor MPI code at the same optimization level. The MPI one-processor code was faster than one-processor shared memory code, as it did not have synchronization steps. The parallel timings were typically faster than the shared memory parallel code, except for the $N = 2 ^ { 1 6 }$ array size for which the shared memory code did slightly better. We find that for distributed memory parallel implementations of SSF over the range of $2 ^ { 1 2 } < N < 2 ^ { 1 8 }$ two-node implementations have good speedup with maximum speedup at $N = 2 ^ { 1 4 }$ , beyond which the communication cost increases and the computation/communication ratio decreases for larger problem size. The communication cost is different in the MPI case than for shared memory, as data is communicated in "communication stages" so less than perfect speedup (SU) is due to the volume of data communicated between processors in redistribution stages. Using four nodes, we find that the speedup increases with the working set $N$ . This is due to both making the computation stages faster O(NLog(N)/8) and by communicating only O(N/16) of data between single processors in the redistribution stage. For small problem size there is not enough work to make dividing the problem among 4 processors beneficial. The speedup in the distributed SSF algorithm is attributed to the independence of data contained in a processor's local cache between data rearrangement stages, which is not true for the dynamic assignment and sharing of subarray data throughout computational stages in shared memory SSF implementations.

# Distributed Memory (MPI)

<table><tr><td rowspan=7 colspan=1>array size (N) andnumber of steps (S)T1pr. (sec)T2pr (sec)T4pr (sec)SU = T2pr/T1pr.SU = T4pr/T1pr.</td><td rowspan=1 colspan=1>N = 212</td><td rowspan=1 colspan=1>N = 214</td><td rowspan=1 colspan=1>N = 216</td><td rowspan=1 colspan=1>N = 218</td></tr><tr><td rowspan=1 colspan=1>S = 8000</td><td rowspan=1 colspan=1>S = 2000</td><td rowspan=1 colspan=1>S = 500</td><td rowspan=1 colspan=1>S = 125</td></tr><tr><td rowspan=1 colspan=1>37.9</td><td rowspan=1 colspan=1>44.5</td><td rowspan=1 colspan=1>59.4</td><td rowspan=1 colspan=1>92.4</td></tr><tr><td rowspan=1 colspan=1>24.7</td><td rowspan=1 colspan=1>25.4</td><td rowspan=1 colspan=1>34.9</td><td rowspan=4 colspan=1>65.926.81.43.4</td></tr><tr><td rowspan=1 colspan=1>18.8</td><td rowspan=1 colspan=1>16.3</td><td rowspan=1 colspan=1>20.1</td></tr><tr><td rowspan=1 colspan=1>1.5</td><td rowspan=1 colspan=1>1.8</td><td rowspan=1 colspan=1>1.7</td></tr><tr><td rowspan=1 colspan=1>2.0</td><td rowspan=1 colspan=1>2.7</td><td rowspan=1 colspan=1>3.0</td></tr></table>

These results are encouraging in that the speedup in multiprocessor SSF implementations is considerable. Speedup over sequential code using vendor optimized full array 1D FFT is even greater. We recommend implementing the parallel SSF algorithm even on sequential machines due to the $1 0 \%$ to 20% speedup over optimized 1D sequential SSF algorithms. This reflects a better use of the L1 cache and data locality by using small subarrays and removing the transposition stage of the 1D FFT in SSF. For shared memory implementations of parallel SSF, the maximum speedup requires balancing contention in the communicating data contained over more than one processor to the computation performance gain of using small subarrays. For the distributed parallel SSF there is more data locality as data is distributed statically prior to the computational stages. This division of computation and communication stages is different than for shared memory SSF which dynamically distributes subarray FFTs and shares data on more than one processor. Distributed SSF speedup is a function of the number of processors $P$ which reduces both the computational time and communication volume between single processors. The speedup of the parallel SSF is strongly dependent on reducing communication time and contention in the multiprocessor.

# 5 Conclusions

Multiprocessor systems occupy the middle ground in computing between sequential and massively parallel computation. In multiprocessor computing, one wants to write code to take advantage of between 2 and 16 processors to get good speedups over sequential code. Our parallel SSF method is designed for small numbers of tightly integrated processors to divide the 1D FFT into many subarray FFTs performed on P processors. The speedup depends on optimizing the computational speed gain to communication cost in order to speedup traditionally sequential numerical code. The shared memory parallel SSF algorithm does not scale with problem size $N$ as subarray data is distributed over more than one processor causing increases in contention due to gathering large amounts of subarray data from many processors. The distributed memory parallel SSF algorithm uses independent data during computational stages and then uses expensive data redistribution stages. The communication cost of the data redistribution stages can be reduced by using more processors, which also decreases the time for the computation stage. Our results suggest that nearly perfect speedup can be achieved over sequential SSF algorithms by tuning the number of processors and problem size. The significant speedup over sequential code is broadly applicable to many sorts of code which depend crucially on speeding up the sequential 1D FFT and should be explored for other numerical algorithms.

# ACKNOWLEDGMENTS

This work has been performed under the auspices of the US Department of Energy under contract W-7405-ENG-36 and the Applied Mathematical Sciences contract KC-07- 01-01.

V. Ruban and A. Zenchuk thank the Theoretical Division Group T-7 and Center for Nonlinear Studies at the Los Alamos National Laboratory for their hospitality and financial support during the summer of 1997. S.M Zoldi was supported by the Computational Science Fellowship Program of the Office of Scientific Computing in the Department of Energy. S.M Zoldi thanks the Center of Nonlinear Studies for their hospitality during the summer of 1997. Useful conversations with Richard Peltz are acknowledged.

# References

[1] A. Hasegawa and Y. Kodama, Solitons in Optical Communication, Oxford Univ. Press (1995).   
[2] G.P. Agrawal, Nonlinear Fiber Optics, Acad. Press, (1994) 2nd edition.   
[3] V.E. Zakharov and A.B. Shabat, Sov. Phys. JETP, 34, 62 (1972).   
[4] "Optical Communications IIIa", ed. by I.P. Kaminow and T.L. Koch, Acad. Press, (1997).   
[5] P.V. Mamyshev and L.F. Mollenauer, Opt. Lett. 21, (6), 396 (1996).   
[6] N.S. Bergano and C.R. Davidson, J. of Lightwave Technology, 14, (6), 1299 (1996).   
[7] P.N. Butcher and D. Cotter, The Elements of Nonlinear Optics, Cambridge Studies in Modern Optics, v.9, (1993).   
[8] M. Snir, S. Otto, S.H Lederman, D. Walker, J. Dongarra, MPI: The Complete Reference, MIT Press, London, (1996).   
[9] T.R Taha and M.J. Ablowitz, J. Comp Physics, 55, 203 (1984).   
10] P.N Swartztrauber, Parallel Computing, 5, 197 (1987).   
11] A. Averbuch, E. Gabber, B. Gordissky, Y. Medan, Parallel Computing, 15, 61 (1990).   
12] A. Dubey, M. Zubair, C.E. Grosch, Parallel Computing, 20, 1697 (1994).