---
title: Computational Physics -- 7 -- The Temperature of Earth
author:
- Ludwig Neste
header-includes: |
    \usepackage{blindtext}
    \usepackage[section, below]{placeins}
    \usepackage{siunitx}
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

# The Physical Problem
The earth's climate is a complex system. 
By climate, we do not mean the small scale weather effects like wind, local temperature, air humidity and so on. 
Rather, we mean the weather and other thermodynamical variables averaged over large times and volumes. 
This makes the problem of predicting the earth's climate easier than predicting tomorrow's weather, since the small variations in e.g. the temperature 
are averaged out. 

In this report we will try to make a first reasonable guess for the earth's temperature, by just modelling the radiation balance in the atmosphere.
The earth's temperature is mainly coming from the sun's radiation. Other effects like cosmic radiation and the temperature of earth's core also have a small contribution, which is negligible for us.
A body of constant temperature, which is only heated up due to outer radiation, has to emit as much radiation $E$ as it 
absorbs $A=E$. 
According to the Stefan-Boltzmann law of radiation, a black body with temperature $T$ emits a radiative energy flux of 
$$
F = \sigma T^4\quad \text{ with } \sigma = \SI{5.67e-8}{\watt\per\meter\tothe{2}\per\kelvin\tothe{4}}.
$$
With this law we can calculate the temperature of the surface of the earth, if no atmosphere is present.
If we assume an average incoming energy flux from the sun of $F=\SI{344}{\watt\per\meter\tothe{2}}$ ^[This is already averaged over the whole surface (day and night).] and a 
loss due to reflection of $\eta = 0.3$ (called 'Albedo effect'), we get a temperature of 
$$
T = \sqrt[4]{\frac{(1-\eta)F}{\sigma}}  = \sqrt[4]{\frac{0.7 \SI{344}{\watt\per\meter\tothe{2}}}{\sigma}} \approx \SI{-18}{\celsius}.
$$

This is obviously too cold, since we didn't take the atmosphere into account and the resulting greenhouse-effect. 
We can think of one atmospheric layer as a semitransparent medium, which absorbs some light and thus heats up and re-emitting the light 
into the infrared spectrum. The infrared light of the layer thus also heat up other layers. 
In our model we will ignore reflected light, only assume two (averaged) light wavelengths (visual and infrared) and ignore any effects 
like air movement, convection. We will also assume that a layer radiates half of its radiation to the top and half to the bottom.
Let's summarize the whole balance for one middle layer with index $i$ (not at the ground nor at space):
Incoming radiation from above: Visible light from the sun $T_{i+1}^{\mathrm{in, V}}$, 
infrared radiation from the upper adjacent cell $E_{i+1}/2$, infrared radiation from all the upper cells not 
already absorbed by other cells $T_{i+1}^{\mathrm{in, IR}}$. 
Incoming radiation from below: Infrared radiation from the lower adjacent cell $E_{i-1}/2$, infrared radiation from all the lower cells not 
already absorbed by other cells $T_{i-1}^{\mathrm{out, IR}}$. 
If we summarize this the cell $i$ absorbs (and thus emits, since we are in equilibrium) in total 
$$
A_i = E_i 
= \left(\frac{E_{i+1}+E_{i-1}}{2}+T_{i-1}^{\mathrm{out, IR}}+T_{i+1}^{\mathrm{in, IR}}\right)
 \left(1-e^{-\sigma_{\mathrm{IR}} \rho_i h}\right)
+  T_{i+1}^{\mathrm{in, V}}\left(1-e^{-\sigma_{\mathrm{V}} \rho_i h}\right).
$$

Where $\rho_i$ is the density of the cell, $h$ it's width, $\sigma_{\mathrm{IR}}$ the cross section for 
the infrared light and $\sigma_{\mathrm{V}}$ the cross section for visual light.
For the sake of simplicity we will calculate the density with the barometric formula
$$
\rho_i = \rho_i \exp\left[\frac{-g_0 M \left(h_i-h_0\right)}{R^* T_0}\right].
$$

The self-consistency equations for the 3 types of transmitted rations are
\begin{align*}
T_{i}^{\mathrm{in, V}} &= T_{i+1}^{\mathrm{in, V}}e^{-\sigma_{\mathrm{V}}\rho_i h}\\
T_{i}^{\mathrm{in, IR}} &= \left(T_{i+1}^{\mathrm{in, IR}}+\frac{E_{i+1}}{2}\right)e^{-\sigma_{\mathrm{IR}}\rho_i h}\\
T_{i}^{\mathrm{out, IR}} &= \left(T_{i-1}^{\mathrm{out, IR}}+\frac{E_{i-1}}{2}\right)e^{-\sigma_{\mathrm{IR}}\rho_i h}.
\end{align*}

The top layer is a special layer, since we have to keep in mind that we assume no infrared radiation is coming from space.
The two bottom-most layers are also special layers, since no radiation is transmitted below the lowest layer and 
the lowest layer radiates *all* (not half) of its radiation to the layer above.

# Numerical Methods
We will solve the system of equations by starting with an arbitrary starting vector and repeatedly
evaluating the equations from previous section until we get to a stable point in the temperatures. 
We can then evaluate Stefan-Boltzman's law in every layer to get the temperature in each layer.


# Results