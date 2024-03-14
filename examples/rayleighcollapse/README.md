This example builds a standalone APECSS code for the Rayleigh collapse of a bubble. 

#### Simple Rayleigh collapse of a spherical bubble with Rayleigh-Plesset
Using [simple.apecss](simple.apecss) with ````./build/rayleigh_apecss -options simple.apecss -tend 1```` runs a simple Rayleigh collapse of a bubble, whereby surface tension and viscosity are neglected. Since there is no damping mechanism, the bubble oscillates with constant amplitude indefinitely.

#### Rayleigh collapse of a cylindrical bubble
Using [cylindrical.apecss](cylindrical.apecss) with ````./build/rayleigh_apecss -options cylindrical.apecss -tend 5e-4```` runs a Rayleigh collapse of a cylindrical bubble using the Gilmore model with the Tait EoS. This example reproduces the case shown in Figure 3(a) of [Zhang et al., _Physics of Fluids_ 35 (2023), 103302](https://doi.org/10.1063/5.0167537).

#### Rayleigh collapse including acoustic emissions
Using [emissions.apecss](emissions.apecss) with ````./build/rayleigh_apecss -options emissions.apecss -tend 0.2```` runs a Rayleigh collapse of a bubble and records its acoustic emissions at different locations in space. The Noble-Abel (NA) EoS is used to describe the gas and the Tait EoS is to describe the liquid. The emissions are modelled based on the Kirkwood-Bethe hypothesis. This example reproduces the case studied in Section IV.A of [Denner & Schenke, _Physics of Fluids_ 35 (2023), 012114](https://doi.org/10.1063/5.0131930).

#### Rayleigh collapse that produces a shock front (including acoustic emissions)
Using [shock.apecss](shock.apecss) with ````./build/rayleigh_apecss -options shock.apecss -tend 0.11```` runs a Rayleigh collapse of a bubble and records its acoustic emissions at different locations in space. Contrary to the previous case, the emitted pressure wave forms a shock front. The Noble-Abel (NA) EoS is used to describe the gas and the Tait EoS is to describe the liquid. The emissions are modelled based on the Kirkwood-Bethe hypothesis.

