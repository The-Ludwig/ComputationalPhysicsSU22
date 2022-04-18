---
title: Computational Physics -- 3 -- Quantum Monte Carlo
author:
- Ludwig Neste
header-includes: |
    \usepackage{multirow}
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
Where we introduced $E(x)=\left(\hat H \psi\right)/\psi$.
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
Thus we can find the ground energy and state with:
$$
  E_0 = \min_{\alpha} \langle E_\alpha \rangle = \min_{\alpha} \langle \psi_\alpha | \hat H | \psi_\alpha \rangle
  = \min_{\alpha} \int_{-\infty}^{\infty} E(x) |\psi_\alpha(x)|^2 \mathrm{d}x.
$$
Where we will evaluate the integral with Monte-Carlo methods.

To find the right form of the test-wavefunction $\psi_\alpha$ for a given $H$, is not trivial and often driven by experience and heuristic arguments. 
In the case of the 1D-Harmonical oscillator $\hat H = \hbar \omega/2 \left(-\partial_x^2+x^2\right)$ (dimensionless), 
$$
\psi_\alpha(x) = N e^{-\alpha x^2}
$$
is a reasonable choice. This is quite obvious, since we know that it is the exact solution for $\alpha=1/2$.
Where the normalization factor $N$ is chosen, such that the wavefunction is normalized to one ($\int |\psi|^2 =1$),
we will ignore the factor for numerical calculations, since for our numerical methods we only need a function which is proportional 
to the probability distribution.
The so-called local energy $E_\alpha(x)$ in this case is
$$
E_\alpha(x) = \hbar \omega \alpha +x^2\left( \frac{1}{2}-2\alpha^2 \right)
$$
which gives the analytical expectation value as
\begin{multline}
\langle E \rangle =\hbar \omega \sqrt{\frac{2\alpha}{\pi}} \int_{-\infty}^{\infty}
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

The two dimensional case we are going to test is that of the two dimensional harmonic oscillator
with two interacting electrons:
$$
\hat H/\hbar\omega = \left(\frac{1}{2}\sum_{z=x,y} \sum_{i=1}^2 -\partial_{z_i}^2 +z_i^2 \right)+\frac{\lambda}{r}.
$$
We have given it in dimensionless form where $\lambda$ is the coupling strength of the electrons and $r=\sqrt{(x_1-x_2)^2+(y_1-y_2)^2}$.
For this case the test function is given as 
$$
\psi_\alpha(z) = \exp\left(-\frac{x_1^2+y_1^2+x_1^2+y_2^2}{2}+\frac{\lambda r}{1+\alpha r}\right).
$$
The second derivative and thus the local energy is analytically calculated from this test function, but is omitted here.

# Numerical Methods
If we have random variable $x$ with probability density $p(x)$, and 
many samples of this variable $x_1, \dots, x_N$, we can 
approximate the probability density with 
$$
p(x) \approx \frac{1}{N} \sum_{i=1}^{N} \delta(x-x_i).
$$
Thus we can approximate the average integral by 
$$
\langle f \rangle = \int_{-\infty}^{\infty} f(x) p(x) \mathrm{d}x
\approx  \int_{-\infty}^{\infty} f(x) \frac{1}{N} \sum_{i=1}^{N} \delta(x-x_i) \mathrm{d}x
= \frac{1}{N} \sum_{i=1}^{N} f(x_i).
$$
If we generate samples $x_i$ out of a given probability density to approximate 
the integral, this is called Monte Carlo integration (and if the probability distribution is 
the square of a wavefunction, it is called Quantum Monte Carlo).

But how do we generate (pseudo) random samples $x_i$ which follow a given distribution $p(x)$?
For special probability distributions there exists highly optimized methods. 
The state-of-the-art method for uniform distributions is the *Mersenne-Prime-Twister*,
for the normal distribution e.g. the Box-Müller-Method exists. 
For arbitrary distributions it is a little more complicated, we will use 
the Metropolis-Algorithm.
Stripped down it works like this:

1. Start with some $x$
2. Generate a random variable $s$ from a symmetrical probability distribution p(s) with mean 0.
3. Set $x':=x+s$
4. If $p(x')/p(x) < 1$
    1. Generate a random variable $e$, which is uniformly distributed in $(0, 1)$
    2. if $p(x')/p(x) < e$, set $x' := x$
5. Set $x := x'$, (save it) and start from 1 with this $x$. 

The algorithm is usually run a few times before saving the values, to make the results more independent of the starting value.
The step size distribution $p(s)$ is one of the easily generated distributions, for example uniformally between $(-c, c)$ or 
just the normal distribution. 
The algorithm can be extended to non- symmetrical step-size distributions, and is then called *Metropolis-Hastings-Algorithm*.

To minimize the average energy depending on $\alpha$, the *golden section* search is used.

# Numerical Implementation
I implemented the metropolis algorithm in a `C++` class called, using uniform step sizes or by choice gaussian step sizes.
The uniform random numbers for the uniform and normal distribution were created using `std::uniform_real_distribution` and 
`std::normal_distribution` from the `<random>` header. 
The seed was chosen with a value from a hardware device (such as the time or cpu temperature).
The algorithm throws the first 100 values away to get roughly independend from the arbitrarily chosen starting position.


I came up wiht the following `C++17` implementation of the golden section search
```C++
std::tuple<double, double> golden_search(std::function<double(double)> f,
                                         double xmin, double xmax, double tol) {
  constexpr double golden_ratio = 1.618033988749895;
  double x1 = xmax - (xmax - xmin) / golden_ratio;
  double x2 = xmin + (xmax - xmin) / golden_ratio;
  double f1 = f(x1);
  double f2 = f(x2);

  if (xmin >= xmax)
    throw std::invalid_argument("xmin must be smaller than xmax");

  while (xmax - xmin > tol) {
    if (f1 < f2) {
      xmax = x2;
      x2 = x1;
      f2 = f1;
      x1 = xmax - (xmax - xmin) / golden_ratio;
      f1 = f(x1);
    } else {
      xmin = x1;
      x1 = x2;
      f1 = f2;
      x2 = xmin + (xmax - xmin) / golden_ratio;
      f2 = f(x2);
    }
  }

  return {xmin, xmax};
}
```

# Results

I first tested if my metropolis algorithm correctly reproduces a Gaussian distribution
and how the step-pdf influences the result. 
In \autoref{fig:test_distr} you can see that for both uniform (between -2 and 2) and Gaussian ($\mu=0$, $\sigma=1$) step-size 
the algorithm correctly reproduces the given normal distribution.
Generating 100000 random numbers with uniform step size took $4.62\pm .09 \mathrm{ms}$ and
with gaussian step size $6.37\pm0.05\mathrm{ms}$. The difference can be explained because 
the generation of gaussian random numbers is more complicated.
The different between the two step sizes is not big. 
In the following I will use always the Gaussian ($\mu=0$, $\sigma=1$) step size.

![Histogram of 100000 generated random variables using the Metropolis algorithm with uniform (in (-2, 2)) and Gaussian ($\mu=0$, $\sigma=1$) step size distribution. \label{fig:test_distr}](build/plots/test_1D_distributions.pdf){ width=80% }

Next thing I tested if the algorithm reproduces the correct 1D-Energy expectation value in 
the harmonic oscilaltor in 
dependence of $\alpha$, since we have an analytical solution for that. 
The results can be seen in \autoref{fig:1D_Energy}. 
In the plot we see that the convidence interval ($=\pm 1.96 \sqrt{\mathrm{variance}}$)
is zero at the minimum, since the minimum here is the analytical solution, so the energy is an 
eigenvalue and thus does not depend on the position. 

![The calculated average energy of the 1D-harmonic ocsillator against test-function parameter $\alpha$. The test function is $\psi_\alpha = \exp(-\alpha z^2)$. Results compare monte carlo methods to analytical solution. The 95% confidence interval is shown. For each of the (300) displayed $\alpha$ the average is calculated using 100000 samples. \label{fig:1D_Energy}](build/plots/1d_energy.pdf){ width=80% }

Now employing the golden search (with tolerance $10^{-6}$, and 1million monte carlo samples for each average energy calculation) for the 1D case gives the interval 
(0.502978, 0.502979) which corresponds to $0.500004 < E <0.500006$.
We thus get approximately the correct solution $E_0=\hbar\omega/2$.

We now move on to the 2D-Harmonic oscillator with interacting electrons. Here generating 100000 
values with the metropolis algorithm took $7.0\pm 0.4 \mathrm{ms}$.
Visualizing the generated distribution is a little harder for this case, since 
to plot the dependence of the probability distribution of each of the four variables, we need
a five dimensional plot. Thus we need to reduce the dimensionality.
In \autoref{fig:dens_p1} to \autoref{fig:dens_r} the probability density for each particle is shown, 
but with this visualization we loose the relationship between the particle positions (they should repell each other).
For that reason in \autoref{fig:dens_distance} to \autoref{fig:hex2} the 
probability distribution of the distance of the particle (and in dependence of other variables) is shown.

![The probability density of the fist particle in the 2D harmonic oscillator with $\lambda = 1$ and $\alpha$-optimized (using golden search). \label{fig:dens_p1}](build/plots/density_first_particle.pdf){ width=50% }

![The probability density of the second particle in the 2D harmonic oscillator with $\lambda = 1$ and $\alpha$-optimized (using golden search). \label{fig:dens_p2}](build/plots/density_second_particle.pdf){ width=50% }

![The probability density of the \text{\color{blue} first } and \text{\color{orange} second} particle in the 2D harmonic oscillator with $\lambda = 1$ and $\alpha$-optimized (using golden search).\label{fig:dens_comb}](build/plots/density_combined.pdf){ width=80% }

![The probability density of the distance to the origin of the particles in the 2D harmonic oscillator with $\lambda = 1$ and $\alpha$-optimized (using golden search).\label{fig:dens_r}](build/plots/r_particles.pdf){ width=80% }

![The probability density of the distance of the particles in the 2D harmonic oscillator with $\lambda = 1$ and $\alpha$-optimized (using golden search).\label{fig:dens_distance}](build/plots/r_distance.pdf){ width=80% }

![The probability density of the distance of the particles in dependence of the first particle location in the 2D harmonic oscillator with $\lambda = 1$ and $\alpha$-optimized (using golden search).\label{fig:hex1}](build/plots/hex_r1.pdf){ width=80% }

![The probability density of the distance of the particles in dependence of the second particle location in the 2D harmonic oscillator with $\lambda = 1$ and $\alpha$-optimized (using golden search).\label{fig:hex2}](build/plots/hex_r2.pdf){ width=80% }


The next thing I tested is if the Monte Carlo method gives the correct energy distribution as a function for alpha in 
the case of $\lambda = 0$. Since the energy does not depend on alpha, \autoref{fig:lam0} is very boring, but gives the correct answer.

In \autoref{fig:lam1} to \autoref{fig:lam8} you can see different versions of the average energy vs $\alpha$ plot, 
with different $\lambda$ and with and without the variance.
The plots all are produced with 100000 Monte Carlo samples for each $\alpha$.

The golden search algorithm (with tolerance $10^{-6}$, and 1million monte carlo samples for each average energy calculation) gave the following intervals:

| $\lambda$ | min$\alpha$ | max$\alpha$  | min$E$   | max$E$| 
| ----        | -----      |----         | ---     | ---- |
| 1           |0.383657    | 0.383658    | 3.00028 | 3.00032|
| 2           |0.473603    | 0.473603    | 3.72358 | 3.72346|
| 8           |0.625232    | 0.625233    | 6.64655 | 6.64521|

![Energy as a function of $\alpha$ for $\lambda = 0$ \label{fig:lam0}](build/plots/2d_energy_0.pdf){ width=80% }

![Energy as a function of $\alpha$ for $\lambda = 1$. The $1\sigma$ interval is shown. \label{fig:lam1}](build/plots/2d_energy_1.pdf){ width=80% }

![Energy as a function of $\alpha$ for $\lambda = 1$. The $1\sigma$ interval is shown. \label{fig:lam1_var}](build/plots/2d_energy_1_variance.pdf){ width=80% }

![Energy as a function of $\alpha$ for $\lambda = 2$. The $1\sigma$ interval is shown. \label{fig:lam2}](build/plots/2d_energy_2_variance.pdf){ width=80% }

![Energy as a function of $\alpha$ for $\lambda = 8$. The $1\sigma$ interval is shown. \label{fig:lam8}](build/plots/2d_energy_8_variance.pdf){ width=80% }


# Final Remarks

A more complicated test function and multiple parameter test-functions are interesting things to explore in the future.
Understanding how exactly the test function has to be chosen also remains as a task.