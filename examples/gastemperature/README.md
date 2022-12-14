This example builds a standalone APECSS code for an ultrasound-driven bubble, with the addition of the gas-energy model of [Stricker et al., _J. Acoust. Soc. Am._ 130, 3243-3251](https://doi.org/10.1121/1.3626132). To this end, the additional ODE describing the gas temperature is added to the set of ODEs that APECSS solves. The provided [run.apecss](./run.apecss) file reproduces the case shown in Fig. 7 of the paper of Stricker et al. using the following execution command: ````./build/gastemperature_apecss -options run.apecss -freq 20e3 -amp 70e3 -tend 0.001````. Please note, the bubble needs a few periods to attain a steady oscillation.

Applying the single-ODE energy model of Stricker and co-workers, the gas temperature $T_\mathrm{G}$ is computed as
$$\dot T_\mathrm{G} = \dfrac{Q_\mathrm{cond} - p \dot{V}}{c_{v,\mathrm{g}} m_\mathrm{G}},$$
where $m_\mathrm{G}$ is the gas mass contained in the bubble. This ODE is solved together with the ODEs of the Rayleigh-Plesset models and, if applicable, any other ODEs. The net heat absorbed by the gas due to conduction, assuming the temperature of the liquid is constant and corresponds to the ambient temperature $T_\infty$, is
$$Q_\mathrm{cond} = 4 \pi R^2 k_\mathrm{g}  \frac{T_\infty-T_\mathrm{G}}{l_\mathrm{th}},$$
where the thickness of the thermal boundary layer is estimated as
$$l_\mathrm{th} = \mathrm{min} \left(\sqrt{\frac{k_\mathrm{g} R}{\rho_\mathrm{G} c_{v,\mathrm{g}} \Gamma_\mathrm{g} |\dot{R}|}}, \frac{R}{\pi}\right).$$
