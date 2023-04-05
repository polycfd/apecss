---
title: 'APECSS: A software library for cavitation bubble dynamics and its acoustic emissions'
tags:
  - Python
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Fabian Denner
    orcid: 0000-0001-5812-061X
    # equal-contrib: true
    affiliation: 1
  - name: Sören Schenke
    orcid: 0000-0001-8765-3722
    # equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: Department of Mechanical Engineering, Polytechnique Montr\'eal, Montr\'eal, H3T 1J4, QC, Canada
   index: 1
 - name: Chair of Mechanical Process Engineering, Otto-von-Guericke-Universit\"{a}t Magdeburg, 39106 Magdeburg, Germany
   index: 2
date: 04 April 2023
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The dynamics of cavitation bubbles and the acoustic emissions they produce are important in a broad range of engineering applications and natural phenomena, either because the strong energy focusing of the bubble collapse is to be avoided, as it may cause damage to surfaces, or to be exploited, such as in emerging medical applications. APECSS (Acoustic Pulse Emitted by Cavitation in Spherical Symmetry) is a software library to simulate the dynamic behavior and acoustic emissions of cavitation bubbles using an efficient state-of-the-art numerical framework. APECSS supports different Rayleigh-Plesset models for bubble dynamics in incompressible and compressible media with Newtonian or viscoelastic rheology, considering clean or coated bubbles. Acoustic emissions may be modeled under different modeling assumptions using a tailored Lagrangian wave tracking method, including the formation and attenuation of shock waves. APECSS can be extended easily to include custom functionality and may be incorporated into other software frameworks.

# Statement of need

The pressure-driven dynamics of bubbles, a process commonly referred to as _cavitation_, and the acoustic emissions these bubbles produce play a central role in a large variety of engineering applications and natural phenomena [@Plesset1977;@Lauterborn2010]. Cavitation drives material erosion of propellers and hydro turbines [@Blake1987;@Reuter2022], can be used to improve the fatigue strength of materials by peening [@Gu2021;@Soyama2022] or to clean surfaces and membranes [@Reuter2017a], and is utilized as a building block of smart materials [@Athanassiadis2022] as well as in emerging diagnostic and therapeutic medical applications [@Wan2015;@Kooiman2020]. The extreme conditions produced during a violent bubble collapse, with pressures of $\mathcal{O}(10^{10}) \, \mathrm{Pa}$ and gas temperatures in excess of $10^4 \, \mathrm{K}$, in reproducible benchtop experiments promotes cavitation as a microlaboratory for the study of high-pressure and high-temperature fluid dynamics [@Brenner2002;@Liang2022], and attracts interest as a microfluidic facility for sonochemistry [@Meroni2022;@Tandiono2011] and material synthesis [@Barcikowski2019]. Cavitation is also key to understanding many natural phenomena, such as the hunting behavior of pistol shrimps, which generate collapsing cavitation bubbles to stun their prey [@Koukouvinis2017], or the proliferation of ferns, where cavitation bubbles catapult the spores into the air [@Llorens2016]. Furthermore, the dynamic behavior of bubbles subject to pressure changes is highly nonlinear and exhibits routes to chaos [@Behnia2009]. The acoustic emissions produced by such oscillating or collapsing bubbles are similarly complex [@Liang2022;@Denner2023], being highly nonlinear and admitting the formation of shock fronts. These emissions are, for instance, used in medical applications to enhance the contrast of ultrasound imaging of the vasculature [@Tang2011] or to monitor and control high-intensity focused ultrasound cancer treatments [@OReilly2012], and are studied in the context of material damage from cavitation bubbles [@Gonzalez-Avila2021], geological exploration using seismic airguns [@MacGillivray2019] and the controlled perforation or lysis of biological cells [@Helfield2016;@Tandiono2012]. Although by no means exhaustive, this versatile list of applications illustrates the importance of a comprehensive understanding of cavitation bubble dynamics and its acoustic emissions.

APECSS fills a gap in the open-source software domain by combining established models for bubble dynamics with a state-of-the-art numerical framework for the acoustic emissions that supports different modeling assumptions in a portable software library and with minimal computational overhead. APECSS is written in `C` and, aside from the standard `math` library, has no external dependencies. Compared to the Python and Matlab implementations we used previously for research and teaching, we observe a speed-up of more than $100 \times$ with APECSS for representative examples, enabling, for instance, large-scale parameter studies. The flexible and modular design of APECSS allows to easily extend and change its functionality, and integrate it into other software frameworks. Because of its straightforward installation and ease of use, APECSS is also attractive for undergraduate and postgraduate education, as a virtual laboratory for nonlinear mechanics, chaos and acoustics. 

# Features

At the heart of APECSS lies a numerical solver for the pressure-driven radial dynamics of a bubble described by a Rayleigh-Plesset-type equation, such as the standard Rayleigh-Plesset equation for bubbles in incompressible media [@Rayleigh1917;@Plesset1949], the Keller-Miksis equation for bubbles in weakly-compressible media [@Keller1980], or the Gilmore equation for bubbles in compressible media [@Gilmore1952;@Denner2021]. In addition, models for viscoelastic media (e.g.~Oldroyd-B [@Jimenez-Fernandez2005], Kelvin-Voigt [@Yang2005b] and Zener [@Hua2013] models) and for lipid monolayer coatings of the bubble [@Marmottant2005;@Guemmer2021] are available in APECSS. The radial bubble dynamics are solved using a custom implementation of the embedded Runge-Kutta RK5(4) scheme of [@Dormand1980]. The numerical framework underpinning APECSS is agnostic to the applied equations of state (EoS) for the gas and, if applicable, the compressible liquid [@Denner2023], whereby the ideal-gas EoS (gas), with or without van-der-Waals hard-core radius, the Tait EoS (liquid) [@Cole1948;@Gilmore1952], the Noble-Abel EoS (gas) [@Toro1999;@Denner2021], and the Noble-Abel-stiffened-gas EoS (gas and liquid) [@LeMetayer2016;@Denner2021] are readily supported.

The acoustic emissions of an oscillating or collapsing bubble can be simulated assuming the surrounding medium is incompressible, with $c \rightarrow \infty$, or compressible, with $c=\mathrm{constant}$ or $c=f(p)$, where $c$ is the speed of sound and $p$ is the pressure of the medium surrounding the bubble. For cases in which acoustic waves are emitted into a compressible medium with finite speed of sound, the information associated with the emissions is tracked in the radial direction by a Lagrangian wave tracking method [@Denner2023], which we developed specifically for the acoustic emissions of cavitation bubbles and which is, at present, a unique feature of APECSS. Explicit expressions for the radial position, local velocity and pressure are available for small and moderate Mach numbers, assuming constant fluid properties (quasi-acoustic assumption) [@Trilling1952;@Gilmore1952] or pressure-dependent fluid properties [@Denner2023]. For large Mach numbers, including trans- and supersonic flows, the radial coordinate is integrated in time and the local velocity is determined by either integrating its spatial [@Gilmore1952] or temporal [@Hickling1963] derivative along the outgoing characteristic, with the enthalpy and pressure readily obtained by explicit expressions [@Denner2023]. The formation and attenuation of shock fronts, should they occur, is accounted for by treating multivalued solutions at runtime.

APECSS provides a compact and efficient simulation tool for the prediction of bubble dynamics and its acoustic emissions for research, as demonstrated by the studies that have driven and accompanied the development of APECSS [@Denner2021;@Guemmer2021;@Denner2023], and higher education. A variety of results of the bubble dynamics and the acoustic emissions can be readily written into text files at the end of a simulation or at predefined intervals, as illustrated by the examples shown in Figure \autoref{fig:1}. The APECSS repository includes representative examples that introduce the capabilities of the software library and serve to test the correct functionality of APECSS. The examples also demonstrate how APECSS may be extended with additional functionality or be used to build models for more complex scenarios, such as incorporating an energy model for the gas temperature or modeling the acoustic interaction of multiple bubbles. Designed as a software library, the capabilities of APECSS can be adopted in other software frameworks, and its modular design enables a straightforward extension with custom functionality. A software tool for the simulation of bubble dynamics, cavitation processes and the associated acoustic emissions with features and attributes comparable to APECSS is currently not available open source.


![Results of an argon bubble with initial radius $R_0 = 5 \, \mu \mathrm{m}$ in water, driven by ultrasound with a frequency of $23.5 \, \mathrm{kHz}$ and a pressure amplitude of $145 \, \mathrm{kPa}$, as previously considered by [@Holzfuss2010] in the context of sonoluminescence. (a)-(c) The bubble radius $R(t)$, bubble-wall Mach number $M(t)=\dot{R}(t)/c_\mathrm{L}(t)$, where $c_\mathrm{L}(t)$ is the speed of sound of the liquid at the bubble wall, and gas pressure $p_\mathrm{G}(t)$ as a function of time $t$. (d)-(e) The pressure amplitude $\Delta p(r,t)$ and velocity $u(r,t)$  of the acoustic wave generated by the primary collapse of the bubble as a function of the radial coordinate $r$. (f) The pressure amplitude $\Delta p(r,t)$ emitted by the bubble at a fixed radial distance $r=100 \, \mu \mathrm{m}$ from the bubble center as a function of time $t$. (g)-(h) Spatial profiles of the pressure amplitude $\Delta p(r,t)$ and the velocity $u(r,t)$ at selected time instances, where $t_0$ is the time at which the bubble assumes its minimum radius. The radial bubble dynamics are simulated using the Gilmore-NASG model [@Denner2021] and the acoustic emissions are simulated using the in-built Lagrangian wave tracking method [@Denner2023], integrating ordinary differential equations for $r(t)$ and $u(r,t)$ along the outgoing characteristic and treating the multivalued solutions associated with the formed shock front at runtime. Argon is modeled by the ideal-gas EoS with a polytropic exponent of $1.666$ and a van-der-Waals hard-core radius of $R_0/8.86$ [@Holzfuss2010], and water is modeled using the Noble-Abel-stiffened-gas EoS with a polytropic exponent of $1.11$, a Tait pressure constant of $6.48 \times 10^8 \, \mathrm{Pa}$, a co-volume of $6.8 \times 10^{-4} \, \mathrm{m}^3\mathrm{/kg}$ and a reference density of $997 \, \mathrm{kg/m}^3$ [@Denner2023]. The ambient and reference pressure is $10^5 \, \mathrm{Pa}$.\label{fig:1}](Fig1.png){ width=100% }


# Acknowledgements

We thank Jonas Gümmer for his help with the Python scripts during the initial development stages and Rishav Saha for his efforts on testing early implementations of APECSS. This research was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation), grant number 441063377.