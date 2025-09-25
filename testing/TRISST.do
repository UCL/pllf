/*
various ways to analyse TRISST
IW 16sep2025
*/

cd "C:\ian\git\pllf"
myadopath pllf
pda

use TRISST, clear
binreg outcome modality [fw=n], rd
rdci outcome modality [fw=n]

* binreg succeeds:
pllf, list profile(modality): binreg outcome modality [fw=n], rd

* glm fails for parm<=-0.03 ...
pllf, list profile(modality): glm outcome modality [fw=n], family(binomial) link(identity)
pllf, list profile(modality): glm outcome modality [fw=n], family(binomial) link(identity) difficult
pllf, list profile([outcome]modality): glm outcome modality [fw=n], family(binomial) link(identity)

* ... unless you use irls ...
pllf, list profile(modality): glm outcome modality [fw=n], family(binomial) link(identity) irls

* ... or set starting values
pllf, profileoptions(from(0.04, copy)) list profile(modality): glm outcome modality [fw=n], family(binomial) link(identity)
