This example builds a standalone APECSS code for an ultrasound-driven bubble.

#### Lipid-coated microbubble (simple)
This examples simulates a lipid-coated microbubble using the Rayleigh-Plesset model including acoustic radiation damping. The lipid-monolayer coating is modelled using the Marmottant model. Using [lipidcoated_simple.apecss](lipidcoated_simple.apecss) with ````./build/ultrasound_apecss -options lipidcoated_simple.apecss -freq 2.9e6 -amp 130e3 -tend 2e-6````, this example reproduces Figure 5b of [Marmottant et al., _Journal of the Acoustical Society of America_ 118 (2005), 3499](https://doi.org/10.1121/1.2109427). 

#### Lipid-coated microbubble (with emissions)
This examples simulated a lipid-coated microbubble using the Gilmore model and records its acoustic emissions based on the Kirkwood-Bethe hypothesis. The lipid-monolayer coating is modelled using the Marmottant-Gompertz model and the bubble is assumed to contain sulfur hexafluoride (SF6). Using [lipidcoated_emissions.apecss](lipidcoated_emissions.apecss) with ````./build/ultrasound_apecss -options lipidcoated_emissions.apecss -tend 11e-6 -freq 1e6 -amp 600e3````, a bubble excited by ultrasound with a frequency of 1 MHz and a pressure amplitude of 600 kPa is simulated.

#### Sonoluminescence (simple)
Using [sonolumin_simple.apess](sonolum_simple.apecss) with ````./build/ultrasound_apecss -options sonolum_simple.apecss -tend 40e-6 -freq 25e3 -amp 135e3````, this example reproduces Figure 2 of [Denner, _Ultrasonics Sonochemistry_ 79 (2021), 105307](https://doi.org/10.1016/j.ultsonch.2020.105307).

#### Sonoluminescence (with emissions)
Using [sonolumin_emissions.apess](sonolum_emissions.apecss) with ````./build/ultrasound_apecss -options sonolum_emissions.apecss -tend 40e-6 -freq 23.5e3 -amp 145e3````, this example

#### Kelvin-Voigt
Using [kelvinvoigt.apess](kelvinvoigt.apecss) with ````./build/ultrasound_apecss -options kelvinvoigt.apecss -tend 6e-6 -freq 1e6 -amp 3e6````, this example reproduces Figure 5b of [Yang & Church, _Journal of the Acoustical Society of America_ 118 (2005), 3595](https://doi.org/10.1121/1.2118307).

#### Zener
Using [zener.apess](zener.apecss) with ````./build/ultrasound_apecss -options zener.apecss -tend 5e-6 -freq 1e6 -amp 1e6````, this example reproduces Figure 5b of 
[Zilonova et al., _Ultrasonics Sonochemistry_ 40 (2018), 900](https://doi.org/10.1016/j.ultsonch.2017.08.017).

#### Oldroyd-B
Using [oldroydb.apess](oldroydb.apecss) with ````./build/ultrasound_apecss -options oldroydb.apecss -tend 3e-6 -freq 3e6 -amp 400e3````, this example reproduces Figure 1 (De = 1, lower solid line) of [Jimenez-Fernandez & Crespo, _Ultrasonics_ 43 (2005), 643](https://doi.org/10.1016/j.ultras.2005.03.010).
