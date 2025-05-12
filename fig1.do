set scheme s1mono
/*
	Example 1 for pllf paper.
*/
use brcancer, replace
stset rectime censrec
fracgen x1 -2 -0.5
fracgen x6 0.5
pllf stcox x1_1 x1_2 x4a x5e x6_1 hormon, profile(x4a) gropt(saving(fig1, replace))

gph2emf fig1
