/*
various ways to analyse TRISST
IW 16sep2025
*/

cd "C:\ian\git\pllf"
myadopath pllf
pda

use TRISST, clear
cs outcome modality [fw=n]
binreg outcome modality [fw=n], rd
rdci outcome modality [fw=n]

* binreg succeeds:
pllf, list profile(modality): binreg outcome modality [fw=n], rd

* glm fails for parm<=-0.03 ...
* [use iter(50) to make it fail faster]
pllf, list profile(modality): glm outcome modality [fw=n], family(binomial) link(identity) iter(50)
pllf, list profile([outcome]modality): glm outcome modality [fw=n], family(binomial) link(identity) iter(50)

* ... unless you use irls ...
pllf, list profile(modality): glm outcome modality [fw=n], family(binomial) link(identity) irls

* ... or set starting values
pllf, profileoptions(from(0.04, copy)) list profile(modality): glm outcome modality [fw=n], family(binomial) link(identity)
