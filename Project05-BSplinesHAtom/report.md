---
title: Computational Physics -- 5-- The $\mathrm H$-Atom with B-Splines
author:
- Ludwig Neste
header-includes: |
    \usepackage{blindtext}
    \usepackage[section, below]{placeins}
papersize: a4
fontfamily: libertinus
geometry: 
- left=1cm
- right=1cm
- top=2cm
- bottom=2cm
classoption:
- twocolumn
...

# The Problem
The stationary Schrödinger equation 
$$
\hat H(\vec x) \Psi(\vec x)  = E \Psi(\vec x)
$$
governs the behavior of non-relativistic, subatomic particles. 
One of the most common problems is one electron in the potential of a positively charged 
nucleus, Hydrogen-Like-Atoms:
$$
\hat H(\vec x) = -\frac{\hbar^2\nabla^2}{2m_{\mathrm{e}}}-\frac{Ze^2}{4\pi \epsilon_0 |\vec x|},
$$
where $Z$ is the number of protons in the nucleus.
Since this Hamiltonian is spherically symmetric, it is useful to express the wavefunction 
in terms of radial basis functions
$$
\Psi(\vec x) = \Psi(r, \theta, \varphi) =\frac{P_{nl}(r)}{r}Y_{lm}(\theta, \varphi).
$$
Where $P_{nl}/r$ is the radial part of the solution and $Y_{lm}$ are the basis functions on 
the surface of a sphere, the spherical harmonics. With this we can write the 
associated Schrödinger equation as 
$$
\left(-\frac{\hbar^2\nabla^2}{2m_{\mathrm{e}}}-\frac{Ze^2}{4\pi \epsilon_0 |\vec x|}\right)\frac{P_{nl}(r)}{r}Y_{lm}(\theta, \varphi)
= E \frac{P_{nl}(r)}{r}Y_{lm}(\theta, \varphi)
$$

$$
\Leftrightarrow
\left(-\frac{\hbar^2\partial_r^2}{2m_{\mathrm{e}}}+\frac{\hbar^2 l(l+1)}{2m_er^2}-\frac{Ze^2}{4\pi \epsilon_0 r}\right){P_{nl}(r)}
= 
EP_{nl}(r).
$$
This is now effectively a one-dimensional problem, expressed as a second order differential eigenvalue problem.
We can apply the boundary conditions, that $P_{nl}$ vanishes at infinity and zero. 
If we measure $r$ in multiples of the Bohr-radius
$$a_0=\frac{4\pi\epsilon_0\hbar^2}{m_e e^2}$$
and the energy in Hartree 
$$
1\mathrm{Ht}=\frac{\hbar^2}{a_0^2 m_e}
$$
we can write the equation dimensionless:
$$
\left(-\frac{\partial_r^2}{2}+\frac{ l(l+1)}{2}-\frac{Z}{r}\right){P_{nl}(r)}
= 
EP_{nl}(r).
$$
The analytical (bound state) eigenvalues to the equations are 
$$
E_n = -\frac{\mathrm{Ht}}{2m^2}.
$$
We see, that they do neither depend on $l$, nor on $m$. 
For a given $n$, $l$ can take the values $l=0, 1, \dots, n-1$ and $m$ can take
the values $m=-l, -l+1,\dots, -1, 0, 1, \dots, l-1, l$.
The (normalized) eigenvectors (in the above units) are 
$$
\frac{P_{n l}(r)}{r}  = -e^{- Z r / {n}} \sqrt {{\left ( \frac{2 Z}{n} \right ) }^3\frac{(n-l-1)!}{2n{(n+l)!}} } \left ( \frac{2 Z r}{n} \right )^{l} L_{n-\ell-1}^{(2\ell+1)} \left ( \frac{2 Z r}{n} \right ).
$$
Where $L$ are the associated Laguerre polynomials.

All bound eigenenergies are negative, but positive, unbound solutions exist too:
they form the continuum of the free eigenstates of the electron.

# The Numerical Method
We will express the (radial) wavefuncions as splines of order $k$, developed in B-Splines
$$
P_{nl}(r) = \sum_i c_i B_{i, k} (x).
$$
By multiplying the eigenvalue problem from the left with $B_{j, k}$ and integrating over 
the non-zero parts,
we can express the eigenvalue equation as 
$$
\sum_i c_i\hspace{-1em} \int\limits_{t_{\mathrm{max}(i,j)}}^{t_{\mathrm{min}(i,j)+k}}\hspace{-1.5em} B_{j,k}(x) \hat H B_{i, k}(x) \mathrm{d}x
=E\sum_i c_i \hspace{-1em}\int\limits_{t_{\mathrm{max}(i,j)}}^{t_{\mathrm{min}(i,j)+k}}\hspace{-1.5em} B_{j,k}(x)  B_{i, k}(x) \mathrm{d}x
$$
where $t_i$ are the knot-sequences defining the splines, and $\hat H$ is the operator infront of the one dimensional eigenvlaue 
problem.
This equation can be expressed as a generalized eigenvalue problem
$$
\mathbf{H}\vec c = E\mathbf{B} \vec c.
$$
Where $(\vec c)_i = c_i$, $\mathbf H_{i,j}$ is the left integral and $\mathbf B_{i,j}$ is the right integral.
By integrating we took out the continousity of the problem, where we re-introduced the lost complexity 
by multipliying with the additonal B-Spline.

The integration can be done numerically, where a good approate is *Gaussian' Quadrature*. The idea of gaussian quadrature is 
to generalize the $h$-appoximation
$$
\int f(x)\mathrm{d} x \approx \sum_i h f(x_i),
$$
by allowing for a grid which does not evaluate $f$ at point which are constant with appart, 
but choosing the $x_i$ in an optimal way and introducing a generalized weight $w_i$ instead of constant $h$.
$$
\int f(x)\mathrm{d} x \approx \sum_i w_i f(x_i).
$$
In *Gaussian' Quadrature* the points are chosen in a way where they optimally 
integrate polynomials. With that it is possible to integrate polynomials of order $2n-1$ and lower exactly, 
by just evaluating $n$ points. 
In our usecase we see that the integrand in the $\mathbf B$ matrix is a polynomial of $2k-2$,
so if we use $n=k$ points, we can evaluate the integral exactly (ignoring numerical effects and machine accuracy limitations).
The integrands in $\mathbf H$ are no simple polynomial, thus it is smart to use more points.
To find the right $w_i$ and $x_i$, we must find the roots of the Legendre polynomials, for further details, there is much 
explanation in standard literature.

# Numerical Implementation
The problem was solved with the `Spectra`-Library in `C++`, which is build on top of 
the `Eigen`-Library. 
It allows to solve generalized eigenvalue problems with sparse matrices, which is very usefull, since our problem only includes
banded matrices with $2k-1$ nonzero elements in each row. 

# Results
The solutions to the numerical problem are not contraints to be bound (negative eigenvalues) solutions, 
so we also find positive eigenvalues. Since the positive eigenvalues are a continoum, they don't
really tell us anything interesting.
In testing I found out, that the amount of negative eigenvalues we find 
scales with the logarithm of $r_{\mathrm{max}} = \max t_{i}$, not so much with the number 
of knot-points, or their distribution. This is illustarted in \autoref{fig:numknots}.
I can not explain this relationship. Naively I would assume that the bound state 
$n_{\mathrm{max}}$ with the lowest energy, needs at least a potential 
barrier of $\mathrm{Ht}/n_{\mathrm{max}}^2$, to be bound. For $l=0$, this means 
that the central potential must at least be $\mathrm{Ht}/r_{\mathrm{max}}=\mathrm{Ht}/n_{\mathrm{max}}^2$,
which would result in square relationship, which is not supported by my data.

![Number of found negative eigenvalues in dependence of the last physical knot-point. The $x$-axis is logarithmic, so the 
apperant relationship $N\propto \log{r_{\mathrm{max}}}$ can be seen. Knot-Points are exponentially distributed. \label{fig:numknots}](build/plots/rmax_v_nbound.pdf)

To get approximately the first 14 negative eigenvalues, I set $r_{\mathrm{max}}=1000a_0$.
In testing I found out that exponentially distributed knot points work very well. This means that the
density of knot points scales with $e^{-cr}$, where $c$ is some prefactor, in our case 10. 
It is illustrated in \autoref{fig:exp}.
Thus, we settle for exponentially distributed knots for the rest of this report.

![Comparison of eigenvalues with linear and exponentially distributed knot points. This uses 30 knots, $r_{\mathrm{max}}=1000a_0$ and splines of order $k=4$. \label{fig:exp}](build/plots/eigvalues.pdf)


I also compared what difference the order of the spline makes, for both $k=4$ and $k=7$ I got very similar results, only the 
very low energy bound states had a small difference, but not very significanly. This can be seen in \autoref{fig:k}.
Thus we remain at $k=4$ for the rest of the report. 

![Comparison of eigenvalues with splines of order $k=4$ and $k=7$. This uses 30 exponentially distributed knots and $r_{\mathrm{max}}=1000a_0$. \label{fig:k}](build/plots/ks.pdf)

Instead, we could also choose fewer knots and compensate by choosing higher order splines. This trade-off is visualized in \autoref{fig:tradeof}.

![Comparison of eigenvalues with splines of order $k=4$ and $k=21$. This uses 5 exponentially distributed knots and $r_{\mathrm{max}}=1000a_0$. \label{fig:tradeof}](build/plots/few_knots.pdf)

\FloatBarrier

We can also visualize the eigenfunctions. For $l=0$, this is done in \autoref{fig:l0} and for some higher eigenfunctions in \autoref{fig:l0_high}. 
The numerical solutions are normalized with a simple numerical integration.

![Lowest three eigenfunctions for $l=0$. The analytical solutions are shown by the dots. \label{fig:l0}](build/plots/eigenfunctions.pdf)

![Three higher $n$ eigenfunctions for $l=0$. The analytical solutions are shown by the dots. \label{fig:l0}](build/plots/eigenfunctions_high.pdf)

<!-- ![Three higher $n$ eigenfunctions for $l=0$. The analytical solutions are shown by the dots. \label{fig:l0}](build/plots/free.pdf) -->

The energy levels for different $l$ are the same, so to test our solutions for different l, 
we can plot the eigenfunctions. This is done in \autoref{fig:l1} to \autoref{fig:l6}.
We see that for high $l$ the accuracy of our simulation decreases a little bit. Since 
we are nowhere near the computational limitations, we could just increase the $k$ order 
the number of knots, to make the calculation more accurate.

![Lowest three eigenfunctions for $l=1$.\label{fig:l1}](build/plots/eigenfunctions_l1.pdf)

![Lowest three eigenfunctions for $l=2$.\label{fig:l2}](build/plots/eigenfunctions_l2.pdf)

![Lowest three eigenfunctions for $l=6$.\label{fig:l6}](build/plots/eigenfunctions_l6.pdf)


# Remarks

I also tested the solutions for different $Z$, and got the right energy-eigenvalues.
The results for different $Z$ have to be tested further. 
Also having a more not point-wise nucleus potential still needs to be investigated.