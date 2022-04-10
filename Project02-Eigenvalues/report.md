---
title: Computational Physics -- 2 -- Eigenvalue Solver
author:
- Ludwig Neste
header-includes: |
    \usepackage{multirow}
    \usepackage{minted}
abstract: |
  Solving eigenvalue problems is a very common task in computational physics. 
  In this report I develop an iterative method to calculate the smallest eigenvalues of an 
  eigenvalue problem and apply it to the quantum harmonic oscillator. 
  I compare the results with standard methods which calcualte every eigenvalue.
...

# Introduction

One of the most famous places where an eigenvalue equation pops up in physics is the 
stationary
*Schrödinger equation*
$$
\hat H \psi = E_n\psi.
$$
Where $\psi$ lives in an infinite dimensional *Hilbert Space* $\psi \in \mathcal{H}$,
for example $\mathcal{H} = L^2 = \left\{f: \mathbb R^N\to \mathbb R \middle | \quad \int_{\mathbb R} f(x)^2 \mathrm{d}x < \infty \right\}$
The form of the Hamiltonian Operator $\hat H$ for a one dimensional problem looks like this:
$$
\hat H = -\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2} +V(x).
$$
One of the simplest potentials is that of the harmonic oscillator:
$$
V(x) = \frac{1}{2} m \omega^2 x^2.
$$
We can turn this harmonic oscillator problem dimensionless by introducting 
the dimensionless position $z = x \sqrt{\frac{m\omega}{\hbar}}$:
$$
\left(-\frac{\partial ^2}{\partial z ^2} +z^2\right)\psi = \tilde E \psi
$$
where $\tilde E = 2E/(\hbar\omega)$, so measured in units of $\hbar \omega/2$.

We can approximate the eigenvalue equation on a finite grid of $N$ points 
$x_1, \dots, x_N$. In our case we will have an equally sized grid, so 
we always have $x_i-x_{i-1} = h$. 
On this equal sized grid we can approximate the second derivative with 
the five point formula:
$$
\frac{\partial^2 \psi(z) }{\partial z^2} = -\frac{1}{12h^2} \sum_{k=-2}^{2} \tilde c_k \psi(z+k\cdot h) + \mathcal{O}(h^4)
$$
with coefficients $\tilde c_{-2}=\tilde c_{2}=1$, $\tilde c_{-1}=\tilde c_{1}=-16$, $\tilde c_0=30$.
We will use $c_i = \frac{-1}{12h^2}\tilde c_i$. 
Where we will from now on use only the coefficients with non negative indices, since they are symmetric anyway.
Let us define the discrete wavefunction vector as $(\vec\psi\ )_i = \psi(x_i)$, or equally
$$
\vec \psi 
= 
\begin{pmatrix}
\psi(x_1)\\
\psi(x_2)\\
\vdots \\
\psi(x_{N-1})\\
\psi(x_{N})
\end{pmatrix}.
$$
Using the five-point approximation the complete eigenvalue problem can be formulated with
$$
\underbrace{\small
\begin{pmatrix}
-c_0 + V(x_1) & c_1 & c_2 & 0 & 0& 0 &  \\
c_1  & -c_0 + V(x_2) & c_1 & c_2 & 0 & 0 & \\
c_2 & c_1  & -c_0 + V(x_3) & c_1 & c_2 & 0 &  \\
0 & c_2 & c_1  & -c_0 + V(x_4) & c_1 & c_2 &  \\
\multicolumn{6}{c}{} & \ddots & 
\end{pmatrix}
}_{\mathbf{M}}
\vec\psi 
= 
\tilde E_N \vec \psi.
$$
We now transformed the eigenvalue problem to an eigenvalue problem on $\mathbb R^N$, so 
we will also only get $N$ of the infinite eigenvalues.
This also contains the approximation, that the wavefunction is zero
outside $(x_1, x_N)$, which is equal to setting the Potential to infinite 
outside this region.

For the harmonic oscillator we now the analytical solution to the problem to be 
$$
\tilde E_N = 2N+1.
$$
These are an infinite amount of solutions, and if we choose our approximation range around the 
minimum of the potential, we expect the first few lowest numerically approximated eigenvalues 
to be close to the first analytical ones.


# Numerical Methods
In physics, we call the eigenfunction with the lowest eigenvalue (Energy), the *ground state*. 
We are often interested in the ground state and a few of the states close to the ground state, 
so we would like a numerical method, which instead of solving the whole eigenvalue-problem, 
only calculates a few of the lower order eigenvalues.

## Power Iteration
The eigenvectors $\vec v_1, \dots, \vec v_N$ of a symmetric, real $N\times N$ matrix $\mathbf M$ with full rank form an orthogonal basis of $\mathbb R^N$.
This means that every vector $\vec u$ can be written as a linear combination of the eigenvectors
$$
\vec u = 
\sum_{i=1}^{N} a_i \vec v_i.
$$
Because of the eigenvalue equation $\mathbf M \vec v_i = \lambda_i \vec v_i$, we get:
$$
\mathbf M^k \vec u = 
\sum_{i=1}^{N} a_i \lambda_i^k \vec v_i.
$$
If we choose the eigenvectors so, that their length is one, the length of this 
new vector is 
$$
\left | \mathbf M^k \vec u\right | = 
\sqrt{
\sum_{i=1}^{N} a_i^2 \lambda_i^{2k}
}
=|\lambda_\mathrm{max}|^k \sqrt{
\sum_{i=1}^{N} a_i^2 \left(\frac{\lambda_i}{\lambda_\mathrm{max}}\right)^{2k}.
}
\overset{k\to\infty}{\rightarrow}
|a_\mathrm{max}|.
$$
Where the index $\mathrm{max}$ is defined as 
$$
\mathrm{max} = \mathrm{argmax}_k |\lambda_k|.
$$
The limit holds, because every eigenvalue fraction is smaller than one, because we chose the 
eigenvalue which has the largest absolute size as the denomitor, except for that eigenvalue itself. 
The power multiplication can be done recursively with
$$ 
\vec u_n = \mathbf{M} \vec u_{n-1} \quad  \text{ with } \vec u_0 = \vec u,
$$
or normalized:
$$ 
\vec u_n = \frac{\mathbf{M} \vec u_{n-1}}{\left|\mathbf{M} \vec u_{n-1}\right|} \quad  \text{ with } \vec u_0 = \vec u.
$$
Then we have
$$
\vec u_n = \frac{
\sum_{i=1}^{N} a_i \left(\frac{\lambda_i}{|\lambda_\mathrm{max}|}\right)^k \vec v_i
}{ \sqrt{
\sum_{i=1}^{N} a_i^2 \left(\frac{\lambda_i}{\lambda_\mathrm{max}}\right)^{2k}
}
}
\overset{k\to\infty}{\rightarrow}
\vec v_{\mathrm{max}}.
$$
This convergence rate depends on the fraction of the second biggest eigenvalue, compared to the biggest eigenvalue,
in absolute terms.


## Inverse Power Iteration
Since we are interested in the ground state, we want to calculate not the highest eigenvalues,
but the lowest ones^[To be correct: highest/lowest in terms of absolute values. But since we have a positive-definite matrix here, all eigenvalues are positive and I will omit this difference.].
The eigenvalues of the inverse matrix are the reciprocals of the eigenvalues of the matrix, and the 'inverse' 
eigenvectors still form an orthonormal basis.
If we do the power iteration with the inverse matrix, we thus get the eigenvectors with the smallest eigenvalues:
$$ 
\vec u_n = \mathbf{M}^{-1} \vec u_{n-1} 
\overset{k\to\infty}{\rightarrow}
\vec v_{\mathrm{min}}
\quad  \text{ with } \vec u_0 = \vec u.
$$
The iteration can also be calculated by solving the linear equation
$$
\mathbf{M}\vec u_n = \vec u_{n-1} 
$$
for $\vec u_n$. This is very efficient, since if we calculate the LU-decomposition for $\mathbf{M}$ once,
we can solve the system very fast. 
We get the eigenvalue by applying the matrix once more:
$$
\lambda_\mathrm{min} = \vec v_{\mathrm{min}}^{\mathrm{T}} \mathbf{M} \vec v_{\mathrm{min}},
$$
asuming we normalized the eigenvectors in every step.
This is equivalent to
$$
\frac{1}{\lambda_\mathrm{min}} = \vec v_{\mathrm{min}}^{\mathrm{T}} \mathbf{M}^{-1} \vec v_{\mathrm{min}}.
$$

If we want to get the second, third, \dots -lowest eigenvalues, we have to shift the matrix by approximately 
the same eigenvalue, so it is now the lowest eigenvalue of the new matrix:
$$
\mathbf{M}_s = \mathbf{M} -s\mathbb{I}.
$$
This shift does not change the eigenvectors, but it changes the eigenvalues, so we 
have to add back the shift, to get the eigenvalues of the original matrix.


# Programming
Since the main system of equation involves the matrix $\mathbf{M}$, which 
has only a fraction of about $5/N$ entries filled, it is very sparse for big $N$.
So in programming I used the sparse matrices of the \texttt{Eigen3 C++} library.
To compare the power iteration method with the full solution to the eigenvectors 
I used the \texttt{SelfAdjointEigenSolver} of \texttt{Eigen3}, which 
relies on the fact, that the matrix is self-adjoint. Sadly it is not able 
to make use of the sparse feature of the matrix, so its runtime does not scale good with higher 
$N$.

The main algorithm of the inverse power iteration is shown in the following code-snippet:
```C++
std::tuple<VectorXd, double, unsigned int> 
inverse_power_iteration(SparseMatrix<double>& matrix, 
                        VectorXd& start_vector, 
                        double shift = 0,
                        double tol = 1e-4, 
                        unsigned int max_iterations = 10000) 
{
  // first do the (sparse) LU decomposition
  Eigen::SparseLU<Eigen::SparseMatrix<type>, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(mat);
  solver.factorize(mat);

  // now the iteration
  Eigen::VectorXd y_now(start_vector);
  Eigen::VectorXd y_last;

  // Set the init value to maximal allowed  double,
  // so it does not converge on accident
  double lambda_now = std::numeric_limits<double>::max();
  double lambda_last;

  // keep track of iterations
  unsigned int i = 0;

  do {
    // normalize vector
    y_now.normalize();
    y_last = y_now;

    // iteration
    y_now = solver.solve(y_last);

    // calculate eigenvector
    lambda_last = lambda_now;
    lambda_now = 1 / (y_last.dot(y_now)) + shift;

  } while ((std::abs(lambda_last - lambda_now) > tol) && i++ < max_iterations);

  return {y_now, lambda_now, i};
}
```

In each iteration it compares the eigenvalue to that of the previous iteration,
and stops the iteration if it does not change more than some tolerance `tol`. 
It also stops the iteration after some predefined amount of iterations, should it 
not converge. It returns the eigenvalue, eigenvector and number of needed iterations.

Apart from that general layout, I also measured the time each function needs.
Since I used `gcc`s `g++` compiler, I also activated the 
best optimization settings with the command line options `-O3 -ffast-math`. 


# Results

I started by discretising the (dimensionless) space from -10 to 10, with different amount of points $N$.
I got the five lowest eigenvalues with shifts $1, 3, 5, 7, 9$ and set the tolerance 
to $10^{-6}$. 
In the following table you can see how much faster the iterative method for the 
small eigenvalues is, compared to the full diagonalization. It also shows the number 
of iterations took to converge (given the tolerance) for each shift.
I aborted the execution of the 10000 point full diagonalization, since 
it took to long. 

| N | $t_{\mathrm{Eigen}}/\mathrm{ms}$ | $t_{\mathrm{Iterative}}/\mathrm{ms}$ | Iterations |
| ----------- | ----------- | ----------- | ----------- | 
|100|2.9|8.9|2,6,2,7,2|
|1000|1320|3.1|1,3,1,4,2|
|10000|XXX|41.1|1,2,1,2,1|

Both algorithms gave the correct eigenvalues of $1, 3, 5, 7, 9$ to at least 6 decimal 
places for $N\geq1000$. The $N=100$ solution is correct to about 3 decimal places.

The comparison is a bit unfair, since we are only calculating 5 eigenvalues, instead of every 
eigenvalue, but this is usually what we are interested in.
The comparison is also unfair, since we fine-tuned the shifts to the eigenvalues, 
in reality we might not now them.
The bigger 
Admitedly, I chose these shifts because they are exactly the five lowest eigenvalues, which is 
cheating a little bit. 
But if I use not-fine tuned shifts of $1.5, 3.5, 5.5, 7.5, 9.5$, 
the overall picture stays the same:

| N | $t_{\mathrm{Eigen}}/\mathrm{ms}$ | $t_{\mathrm{Iterative}}/\mathrm{ms}$ | Iterations |
| ----------- | ----------- | ----------- | ----------- | 
|100|2.9|1.2|5,16,5,16,5|
|1000|1320|3.4|5,16,5,16,5|
|10000|XXX|50.5|5,16,5,16,5|

The wavefunctions belonging to the five lowest eigenvalues and the potential can be seen in 
\autoref{fig:wavefunctions}. They are normalized, so that
$$
\int_{-\infty}^{\infty} |\psi|^2 \mathrm{d}x \approx h{\vec \psi} ^2 = 1.
$$
Since we see, that the wavefunctions are practically 0 outside the range of $(-5,5)$, we 
will restrict the range in future testing to $(-6,6)$, which had no 
noticeable influence on the quality of the eigenvalues in testing.

![Normalized wavefunction (labeled, right scale) \label{fig:wavefunctions}. The blue parabola is the potential (left scale).](build/plots/plot.pdf){ width=80% }

Now we want to test the solver for a problem, we can't solve analytically, so we 
introduce a disturbance to our original harmonic oscillator potential, and see how it affects the 
solution:
$$
V_{\mathrm{disturbed}} = z^2+C_1e^{-C_2z^2}.
$$
The different values I used are listed in the following table:

| $C_1$ | $C_2 |
| ----- | -----|
|    1 | 10 |
|    5 | 10 |
|    50 | 100 |
|    5 | 1 |

The plots for the resulting wavefunctions, can be seen in \autoref{fig:disturbed_small}
to \autoref{fig:disturbed_wide}. 
If the disturbance is too big, we don't reach different wavefunction with the used shifts,
so some of the wavefunctions are actually the same and overlap. With more time 
differents shifts should be used.
The results from the complete solver of \texttt{Eigen} and the iterative method agree, 
they can't be compared to an analytical solution, though. Looking at the picture, the 
solutions seem to make sense: the small bumb in the middle results 
in that the wavefunctions have smaller values in that region, since it needs to tunnel through 
the small bumb. 

![Normalized wavefunction (labeled, right scale) \label{fig:disturbed_small}. The potential (blue) is the disturbed potential with $C_1=1$ $C_2=10$ (left scale).](build/plots/disturbed_small.pdf){ width=80% }

![Normalized wavefunction (labeled, right scale) \label{fig:disturbed_big}. The potential (blue) is the disturbed potential with $C_1=5$ $C_2=10$ (left scale).](build/plots/disturbed_big.pdf){ width=80% }

![Normalized wavefunction (labeled, right scale) \label{fig:disturbed_pointy}. The potential (blue) is the disturbed potential with $C_1=50$ $C_2=100$ (left scale).](build/plots/disturbed_pointy.pdf){ width=80% }

![Normalized wavefunction (labeled, right scale) \label{fig:disturbed_wide}. The potential (blue) is the disturbed potential with $C_1=5$ $C_2=1$ (left scale).](build/plots/disturbed_wide.pdf){ width=80% }


I tested every solution vector, if they are orthogonal to each other, and indeed they are, 
ignoring numerical fluctiations.


# Conclusion
As expected the iterative method is way faster than the 
complete eigenvalue solver, if $N$ is big enough.
Since 
we use sparse matrices, the scaling with $N$ is also better.

There are still some things to test and implement, if I had more time:
For example one could test other analytically solvable potentials and compare 
the solutions.
The biggest thing left to do is to figure out reasonable shifts to use, without having 
an analytical solution. A starting point would be to calculate the smallest eigenvalue with 
the inverse iterative method, and the biggest eigenvalue with the iterative method. 
Then we know in which range the eigenvalues lay, and we can, for example, divide this range 
by $N$ and thus we have the shifts. 
With this we will probably skip some eigenvalues, so there needs to be more thorough thought, of 
how to precisely choose the shifts.

For now I also just use a 1-vector as the starting value, there is surely a smarter way of choosing a starting vector.




<!-- # References -->