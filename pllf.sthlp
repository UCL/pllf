{smcl}
{* 12may2025}{...}
{hline}
help for {hi:pllf}{right:Patrick Royston}
{hline}


{title:Profile log likelihood function}


{phang2}
{opt pllf}
{it:regression_cmd}
{it:regression_cmd_stuff}
{ifin}
{weight}
[{cmd:,} {it:options}]


{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :{it:Syntax 1}}
{synopt :{opt pro:file(xvarname)}}PLL is required for variable {it:xvarname}, or{p_end}
{synopt :{opt pro:file([eqname]paramname)}}PLL is required for parameter
{it:paramname} or {opt [}{it:eqname}{opt ]}{it:paramname}{p_end}
{syntab :{it:Syntax 2}}
{synopt :{opt form:ula(formula)}}defines a transformation involving at least
one variable in the dataset{p_end}
{synopt :{opt range(#1 #2)}}evaluates PLL over {it:#1} <= {opt X} <= {it:#2}{p_end}
{syntab :{it:Options in common}}
{synopt :{opt nocons:tant}}suppresses the regression constant{p_end}
{synopt :{opt dev:iance}}requests minus 2 times PLL function{p_end}
{synopt :{opt diff:erence}}computes PLL minus maximised log likelihood{p_end}
{synopt :{opt gen(beta_var pll_var)}}creates {it:beta_var} and {it:pll_var}{p_end}
{synopt :{opt cilin:es(cline_options)}}specifies rendition of confidence interval{p_end}
{synopt :{opt gropt(cline_opts twoway_opts)}}supplies graph options to enhance PLL plot{p_end}
{synopt :{opt lev:el(#)}}sets the confidence level to {it:#}{p_end}
{synopt :{opt levlin:e(cline_options)}}specifies rendition of horizontal line{p_end}
{synopt :{opt maxc:ost(#)}}sets an upper limit of 2 * {it:#} on the additional evaluations of the PLL{p_end}
{synopt :{opt n:_eval(#)}}evaluates PLL at {it:#} equally spaced {it:X} values{p_end}
{synopt :{opt noci}}suppresses calculation of PLL-based CI{p_end}
{synopt :{opt nodot:s}}suppresses (supposedly entertaining) dots{p_end}
{synopt :{opt nograph}}suppresses the line plot of the results{p_end}
{synopt :{opt pl:aceholder(string)}}sets the placeholder in Syntax 2 to {it:string}{p_end}
{synopt :{it:regression_cmd_options}}options appropriate to {it:regression_cmd}{p_end}


{pstd}
where, in essence, any {it:regression_cmd} for which the parameters are estimated
by maximum likelihood may be used. This includes
{help clogit},
{help cnreg},
{help glm},
{help heckman},
{help logistic},
{help logit},
{help mlogit},
{help nbreg},
{help gnbreg},
{help ologit},
{help oprobit},
{help poisson},
{help probit},
{help regress},
{help reg3},
{help stcox},
{help streg}
{help stpm}
{help stpm2}, and probably others.

{pstd}
All weight types supported by {it:regression_cmd} are allowed; see help
{help weight}.


{title:Description}

{pstd}
{opt pllf} has two basic syntaxes, depending on which of the options
{opt profile()} or {opt formula()} is used. Let us call these
syntaxes 1 and 2.

{pstd}
With syntax 1, {opt profile()} must be specified and {opt formula()} is
not allowed. {it:regression_cmd_stuff} typically takes the simple form
[{it:depvar}] {it:varlist}, although more complex syntax is supported
according to the needs of {it:regression_cmd}. For example, for
{it:regression_cmd} {cmd:ivreg}, {it:regression_cmd_stuff} takes the form
{it:depvar} [{it:varlist1}] {cmd:(}{it:varlist2} {cmd:=}
{it:varlist_iv}{cmd:)}.

{pstd}
{opt pllf} with the {opt profile()} option (syntax 1)
computes the profile log likelihood (PLL) function
for the regression coefficient of a covariate defined by
{cmd:profile(}{it:xvarname}{cmd:)} or a parameter or variable
defined by {cmd:profile(}[{cmd:[}{it:eqname}{cmd:]}]{it:paramname}{cmd:)}
within a model specified by {it:regression_cmd} {it:yvar} {it:xvarlist},
{it:regression_cmd_options}.
{opt pllf} also reports
PLL-based confidence limits, computed by a simple grid search.
{it:xvarname} does not have to be a member of {it:xvarlist}, although
it is harmless to include it. The results
are saved to new variables assigned by the {cmd:gen()} option.
The dataset length is increased if {cmd:n()} exceeds the number
of observations ({cmd:_N}).

{pstd}
With syntax 2, {opt formula()} and {cmd:range()}
must be specified and {opt profile()} is 
not allowed. {it:regression_cmd_stuff} is similar to that for
syntax 1, except that {it:regression_cmd_stuff} must include the
placeholder {cmd:X}. This is substituted by a variable calculated 
according to the formula defined by {opt formula()}, which must 
also include {cmd:X} at least once.

{pstd}
{opt pllf} with the {opt formula()} option (syntax 2)
computes the PLL function of a non-linear parameter denoted
by {cmd:X}. {cmd:X} is symbolically included where necessary
in {it:regression_cmd_stuff}. In effect, {cmd:X} is replaced on the fly
by the variable created by substituting the current value of {cmd:X} in
{it:formula}. {opt pllf} reports the MLE and PLL-based confidence limits,
computed by a simple grid search. Normal-based confidence limits
are not computed. Other features are similar to
those with syntax 1.

{pstd}
Finally, before using {cmd:pllf} please read an important comment
in {it:Remarks} under the heading {it:Handling 'equations' correctly}.


{title:Options}

{phang}
{cmd:profile(}{it:xvarname} | [{cmd:[}{it:eqname}{cmd:]}]{it:paramname}{cmd:)}
(syntax 1) is not optional. In the first format,
the PLL function for the regression
coefficient for {it:xvarname} is calculated. {it:xvarname}
is a covariate in the main response model.
In the second format, the PLL function for the parameter defined
by {cmd:[}{it:eqname}{cmd:]}{it:paramname} is calculated. Typically,
{it:paramname} will be an auxiliary parameter of some kind, such as
a scale or shape parameter, with its own 'equation'. For example,
for the Weibull model, {cmd:profile([ln_p]_cons)} would give
the PLL function for the log of the shape parameter, p.

{phang}
{opt formula(formula)} (syntax 2) is not optional.
{it:formula} defines
a transformation involving at least one variable in the dataset and the
parameter denoted by placeholder {cmd:X}.
{it:formula} may be any valid Stata expression.
Example: {cmd:formula(exp(-X*x5))}.

{phang}
{opt cilines(cline_options)} specifies the rendition of the vertical
lines representing the bounds of the profile-likelihood-based confidence
interval (CI).  See {help cline_options:{it:cline_options}}.

{phang}
{opt deviance} requests minus 2 times the PLL function, i.e.
the profile deviance function. If {cmd:difference}
is also specified, {cmd:deviance} produces the profile deviance
difference, i.e. minus 2 times the PLL difference.

{phang}
{opt difference} computes the PLL function minus the maximised
log likelhood for the model. See also {cmd:deviance}. Except
in pathological cases, the resulting values are negative or zero.
Pathological cases denote likelihood functions with multiple maxima
or no maximum.

{phang}
{opt gropt(cline_options twoway_options)} supplies graph 
options to enhance the plot of the PLL (or a transformation of it)
against beta. The default graph is a line plot ({help twoway line})
showing the PLL-based confidence interval for beta as vertical lines
parallel to the Y-axis and the corresponding PLL value (or a
transformation of it) as a horizontal line parallel to the X-axis.
Appropriate linear transformation of the PLL is applied when the
{cmd:deviance} and/or {cmd:difference} options are specified.  For a more on
these options, see {help cline_options:{it:cline_options}} and 
{help twoway_options:{it:twoway_options}}.

{phang}
{opt gen(beta_var pll_var)} creates two new variables:
{it:beta_var} to contain the values of the regression coefficient
over which the PLL is evaluated, and {it:pll_var}, to contain the PLL values.
If {opt gen()} is not specified, the variables are created
with default names of {opt _beta} and {opt _pll}, respectively.

{phang}
{opt level(#)} sets the confidence level to {it:#}; default is
{cmd:level(95)}.

{phang}
{opt levline(cline_options)} specifies the rendition of the horizontal
line showing the profile-likelihood at the confidence level for the the
profile-likelihood-based CI.  See {help cline_options:{it:cline_options}}.

{phang}
{opt maxcost(#)} sets an upper limit of 2 * {it:#}
on the number of additional evaluations
of the PLL, when searching for the PLL-based confidence limits. You
should rarely if ever need this option. The program is prevented
by {cmd:maxcost()} from cycling "for ever" when trying to find
confidence limits in pathological cases (see
{opt difference}. Default {it:#} is {cmd:n()}/2.

{phang}
{opt n_eval(#)} evaluates the PLL function at {it:#} equally spaced X values;
default is 100.

{phang}
{cmd:noci} suppresses calculation of the PLL-based confidence limits.

{phang}
{opt nodots} suppresses dots. By default, a dot is
displayed at each evaluation of the PLL.

{phang}
{opt nograph} suppresses the line plot of the results.

{phang}
{opt placeholder(string)} defines the placeholder
character(s) used in Syntax 2. Spaces or other punctuation
characters are not allowed. Default {it:string} is {cmd:X}
(capital-x).

{phang}
{opt range(#1 #2)} with syntax 1 evaluates the PLL function
over {it:#1} <= beta <= {it:#2}, where beta is the regression coefficient
for {it:xvarname}. Default is for {it:#1} and {it:#2} to be the
confidence limits for beta defined by the option {cmd:level()} and
the usual assumption of a normal distribution for beta_hat, the
maximum likelihood estimate of beta.

{phang}
{opt range(#1 #2)} with syntax 2 is not optional.
It evaluates the PLL function over {it:#1} <= {cmd:X} <= {it:#2},
where {cmd:X} is the non-linear parameter of interest. {opt pllf}
also seeks the MLE of {cmd:X}, but if the values of {it:#1 #2} are
ill-chosen or the PLL function behaves 'badly',
it may fail to find the MLE, or give an inaccurate estimate.
The most satisfactory situation is when the MLE lies between {it:#1}
and {it:#2}, and this may be judged from the plot of the PLL function.
In many cases, particularly with large sample sizes,
the PLL function is approximately quadratic with a single maximum.

{phang}
{it:regression_cmd_options} may be any of the options appropriate to
{it:regression_cmd}. If supported by {it:regression_cmd}, this could
include {opt offset(varname)}.


{title:Remarks}

{pstd}
The PLL function is used for two purposes: (1) to estimate likelihood-based
confidence intervals for parameters, and (2) to study the behaviour of the
likelihood function in pathological situations, i.e. when there is no
unique maximum or no maximum at all. Departure of the shape of the
PLL function from that of a quadratic indicates non-normality in the 
distribution of the parameter estimate of interest.

{pstd}
Note that sometimes, the MLE cannot be found at some values of the parameter
being profiled. In those cases, the maximization continues by default for
16000 iterations, and that can take a very long time. If the program appears
to "freeze" in this way, it is best to halt it by pressing ctrl/break. Then
restart it with a more suitable range of values (see the {cmd:range()}
option).

{pstd}
The pseudo standard error of beta or {cmd:X} is computed as
upper PLL confidence limit minus lower PLL confidence limit,
divided by twice t, where t is the appropriate quantile
of the t or normal distribution used in calculating normal
based confidence limits. When the sampling distribution of
the parameter of interest is close to normal, the usual standard
error and the PLL-based standard error will be approximately equal.

{pstd}
To get a profile likelihood for the regression constant ({hi:_b[_cons]}),
you should create a new variable equal to 1 for all observations
and include that in the model together with the {opt noconstant} option
(see example below).

{pstd}
Some estimation commands do not return a log likelihood, that is, after
running the command, {cmd:e(ll)} is missing. Such programs may
return a log pseudo-likelihood in the form of a deviance, which is
defined as minus twice the log likelihood. When {cmd:e(ll)} is
missing and {cmd:e(deviance)} is non-missing, {cmd:pllf} utilizes the
log pseudo-likelihood, defined as -0.5*{cmd:e(deviance)}, for
estimating confidence intervals. An example is {cmd:binreg}; if
the {opt ml} option is not specified, it returns {cmd:e(deviance)}
but not {cmd:e(ll)}.

    {bf:{ul:Handling 'equations' correctly}}

{pstd}
Sometimes, when one seems to have a single unnamed 'equation'
whose parameters are to be estimated, Stata has a hidden name
for the equation. It is essential in such cases to specify the
equation name in the {opt profile()} option, otherwise a syntax
error will probably occur. The most obvious example
is the {cmd:poisson} command (also the {cmd:glm}
command). The main equation in {cmd:poisson} and {cmd:glm} is called
{cmd:[}{it:varname}{cmd:]} where {it:varname} is the name of
the dependent variable. For example, the command

{phang}{cmd:. pllf count group, exposure(y) range(0.2 1.6) profile(count)}

{pstd}
fails, whereas 

{phang}{cmd:. pllf count group, exposure(y) range(0.2 1.6) profile([count]group)}

{pstd}
works as intended.

{pstd}
You can discover hidden equation names by running the regression command
in question and then typing {cmd:matrix list e(b)}. Any equation names will
appear in the top row of the description, followed by a colon. For example:

{phang}{cmd:. pllf count group, exposure(y)}{p_end}
{phang}{cmd:. matrix list e(b)}

    {bf:e(b)[1,2]}
    {bf:         count:      count:}
    {bf:         group       _cons}
    {bf:y1   .89797677  -5.0930282}

{pstd}
Unfortunately, there is no simple way around this naming requirement.


{title:Examples}

{phang}Syntax 1

{phang}{cmd:. webuse brcancer, clear}{p_end}
{phang}{cmd:. stset rectime, failure(censrec) scale(365.24)}

{phang}{cmd:. pllf stcox x1 x4a x5e x6 hormon, profile(x5e) range(-3 -1)}

{phang}{cmd:. pllf stpm2 x1 x4a x5e x6 hormon, df(2) scale(hazard) gen(X Y) profile(x1)}

{phang}{cmd:. pllf streg x1 x4a x5e x6 hormon, distribution(weibull) profile([ln_p]_cons) n(50)}

{phang}{cmd:. pllf streg x1 x4a x5e x6 hormon, distribution(weibull) ancillary(x4b) profile([ln_p]x4b) deviance difference n(20)}

{phang}{cmd:. pllf poisson d group, exposure(y) profile([d]group) range(0.27 1.55)}

{pstd}
{cmd:. sysuse auto}{break}
{cmd:. gen one = 1}{break}
{cmd:. pllf logit foreign mpg one, noconstant profile(one)}

{pstd}
{cmd:. clear}{break}
{cmd:. input group pyears events}{break}
{cmd:{space 2}0 200 38}{break}
{cmd:{space 2}1 100 19}{break}
{cmd:{space 2}end}{break}
{cmd:. poisson events group, exposure(pyears)}


{phang}Syntax 2

{phang}{cmd:. pllf logit y x1 X, formula(exp(-X*x2)) range(.05 .25)}


{title:Stored}

{pstd}
{opt pllf} stores in {cmd:r()}. Entries in which only beta is mentioned
apply only to syntax 1, otherwise to both syntaxes:

	scalar {cmd:r(n)}         Number of obsevations in the estimation sample
	scalar {cmd:r(b)}         MLE of beta or {cmd:X}
	scalar {cmd:r(se)}        Usual standard error of beta
	scalar {cmd:r(pse)}       Pseudo standard error of beta or {cmd:X}
	scalar {cmd:r(n_llci)}    Lower normal-based confidence limit for beta
	scalar {cmd:r(n_ulci)}    Upper normal-based confidence limit for beta
	scalar {cmd:r(l_llci)}    Lower PLL-based confidence limit for beta or {cmd:X}
	scalar {cmd:r(l_ulci)}    Upper PLL-based confidence limit for beta or {cmd:X}
	scalar {cmd:r(ll)}        Maximised log likelihood
	scalar {cmd:r(ll_limit)}  Value of log likelihood for likelihood based CI calc
	scalar {cmd:r(cost)}      "Cost" of evaluating PLL-based confidence limits


{title:References}

{phang}
D. J. Venzon and S. H. Moolgavkar. 1988. A method for computing
profile-likelihood-based confidence intervals. Applied Statistics
37: 87-94.

{phang}
M. S. Pearce. 2000. Profile likelihood confidence intervals
for explanatory variables in logistic regression. 
Stata Technical Bulletin STB-56, 45-47; STB reprints
vol 10, 211-214.


{title:Author}

{pstd}
Patrick Royston{break}
MRC Clinical Trials Unit, London NW1 2DA, UK.

{pstd}Email: {browse "mailto:pr@ctu.mrc.ac.uk":patrick.royston@ctu.mrc.ac.uk}


{title:Also see}
