This example builds a standalone APECSS code for the Rayleigh collapse of a bubble. 

#### Simple collapse with Rayleigh-Plesset
Using [simple.apecss](simple.apecss) with ````./build/rayleigh_apecss -options simple.apecss -tend 1```` runs a simple Rayleigh collapse of a bubble, whereby surface tension and viscosity are neglected. Since there is no damping mechanism, the bubble oscillates with constant amplitude indefinitely.

#### Collapse including acoustic emissions
Using [emissions.apecss](emissions.apecss) with ````./build/rayleigh_apecss -options emissions.apecss -tend 0.2```` runs a Rayleigh collapse of a bubble and records its acoustic emissions at different locations in space. The Noble-Abel (NA) EoS is used to describe the gas and the Tait EoS is to describe the liquid. The emissions are modelled based on the Kirkwood-Bethe hypothesis.