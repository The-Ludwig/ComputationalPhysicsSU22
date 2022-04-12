---
title: Computational Physics -- 3 -- Quantum Monte Carlo
author:
- Ludwig Neste
header-includes: |
    \usepackage{multirow}
nocite: |
  @*
abstract: |
    'Don't forget to add your abstract!'
...

# The Physical Problem

In quantum mechanics we do not calculate definite experimental outcomes, but probabilities of 
an experimental outcome. 
We first calculate the wavefunction $\psi$, and then the probability density is given by its absolute 
value squared $|\psi|^2$. 
If we use the wavefunction in position-space $x$^[Opposed to e.g. momentum-space], 
the expectation value of a measured quantity, let's say energy $E(x)$ is given by
$$
\langle E \rangle = \langle \psi | \hat H | \psi \rangle = \int_{-\infty}^{\infty} E(x) |\psi(x)|^2 \mathrm{d}x.
$$
This is of course only the one dimensional case. 
If we have for example a 2-particle problem in e.g. 2 dimensions, the wavefunction and the 
energy depend on both of the positions $\vec r_1=(x_1, y_1)^\mathrm{T}$ and $\vec r_2=(x_2, y_2)^\mathrm{T}$
 of the particles and the expression becomes more complicated:
\begin{multline}
\displaystyle{}
\langle E \rangle 
= \int_{\vec r_1, \vec r_2 \in \mathbb{R}^2} E(\vec r_1, \vec r_2) |\psi(\vec r_1, \vec r_2)|^2 \mathrm{d}r_1^2\mathrm{d}r_2^2\\
= \int\limits_{-\infty}^{\infty}\int\limits_{-\infty}^{\infty}\int\limits_{-\infty}^{\infty}\int\limits_{-\infty}^{\infty} 
E(x_1, y_1, x_2, y_2) |\psi(x_1, y_1, x_2, y_2)|^2 \mathrm{d}x_1\mathrm{d}y_1\mathrm{d}x_2\mathrm{d}y_2.
\label{eqn:energy_expectation_value}
\end{multline}
Depending on the wavefunction this integral may not be solvable analytically, even in the one dimensional case.
For that reason we will use the method of *Monte-Carlo Integration* to evaluate this integral (see next chapter).

We can easily see from equation \eqref{eqn:energy_expectation_value}, that if the system is in an energy eigenstate 
$\hat H \psi = E_n \psi$, then the expectation value is just the energy eigenvalue $\langle E \rangle = E_n$, since the wavefunction is normed. 
But for many (especially multidimensional) Hamiltonians it is not possible to exactly solve for the energy eigenstates.
If we still can approximate the energy ground state, but just assuming a form of the wavefunction $\psi_\alpha(x)$, which depends on 
one or multiple parameters $\alpha$. We can then optimize the parameter(s) to minimize the energy expectation value, to 
find an approximation of the ground state. 
To find the right form of the test-wavefunction $\psi_\alpha$, is not trivial and often driven by experience and heuristic arguments. 




# Numerated Section


The analytical expectation value of the energy is 
\begin{multline}
\langle H \rangle =\hbar \omega \sqrt{\frac{2\alpha}{\pi}} \int_{-\infty}^{\infty}
e^{-2\alpha x^2} 
\left( \alpha +x^2\left( \frac{1}{2}-2\alpha^2 \right)\right)
\mathrm{d}x \\
= 
\hbar \omega\sqrt{\frac{2\alpha}{\pi}} \left( 
  \alpha \sqrt{\frac{\pi}{2\alpha}}
  + \left( \frac{1}{2}-2\alpha^2 \right)
  \frac{\sqrt\frac{\pi}{2}}{4 \alpha^{3/2}} 
\right)
=
\hbar\omega \left( 
  \alpha 
  + \left( \frac{1}{2}-2\alpha^2 \right)
  \frac{1}{4 \alpha} 
\right)
\end{multline}

Math works^[This is a footnote]! 
$$
e^{i\pi} = -1
$$

![This is a figure with a label, so it can be referenced! \label{fig}](build/plots/test.pdf){ width=80% }

Reference to \autoref{fig}, and also citations[@realistic]!

# References