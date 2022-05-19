---
title: Computational Physics -- 08 -- Fourier Transform of Images
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

# Discrete Fourier Transform

The *Discrete Fourier Transform (DFT)* of a series of complex numbers $x_0, \dots, x_{N-1}$
are the $N$ numbers defined as
$$
X_k = \frac{1}{\sqrt{N}}\sum_{n=0}^{N-1} x_n e^{-in\frac{2\pi}{N} k}.
$$
For intuitive understanding one can think of $X_k$ as the intensity of the frequency, 
$\omega = k 2\pi / N$ in the signal. 

For a deeper understanding let's take a step back and remember that Fourier
found that a periodic function can be represented as a weighted sum of cosine and 
sine functions. 
This does not only apply for continuous, but also discontinuous functions, 
such as a square wave. 
To put it in a more mathematical way: The 
$L^2([0, T])$-space^[This means every function defined on $[0, T]$, which's square is integrable over $[0, T]$.]
has the basis functions 
$$
\exp\left({-2ikx\frac{2\pi}{T}}\right) \quad k\in \mathbb{Z}
$$