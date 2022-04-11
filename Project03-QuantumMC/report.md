---
title: Computational Physics -- N -- The Project
author:
- Ludwig Neste
header-includes: |
    \usepackage{multirow}
nocite: |
  @*
abstract: |
    'Don't forget to add your abstract!'
...

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