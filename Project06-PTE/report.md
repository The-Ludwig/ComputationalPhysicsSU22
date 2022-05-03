---
title: Computational Physics -- N -- The Project
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
nocite: |
  @*
abstract: |
    'Don't forget to add your abstract!'
...

# Numerated Section
Math works^[This is a footnote]! 
$$
e^{i\pi} = -1
$$

![This is a figure with a label, so it can be referenced! \label{fig}](build/plots/test.pdf){ width=40% }

Reference to \autoref{fig}, and also citations[@realistic]!

\Blindtext

# References