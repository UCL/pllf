{smcl}
{* 03nov2025}{...}
{viewerjumpto "Description" "pllf##description"}{...}
{viewerjumpto "Options" "pllf##options"}{...}
{viewerjumpto "Remarks" "pllf##remarks"}{...}
{viewerjumpto "Examples" "pllf##examples"}{...}
{viewerjumpto "Stored results" "pllf##stored_results"}{...}
{viewerjumpto "References" "pllf##references"}{...}
{viewerjumpto "Authors" "pllf##authors"}{...}

help for {hi:pllf}{right:Patrick Royston, Ian White}
{hline}


{title:Profile log likelihood function}


{phang}
Syntax 1: fit a regression command and compute the profile log likelihood (PLL) function for one of its parameters

{phang2}
{opt pllf, profile(parameter)} [{it:options}]: {it:regression_cmd}

{phang}
Syntax 2: fit a regression command and compute the profile log likelihood (PLL) function for a parameter in an extra non-linear term,
which is given by a formula involving observed variable(s) and one unknown parameter

{phang2}
{opt pllf, formula(formula)} [{it:options}]: {it:regression_cmd}



{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :{it:Syntax 1}}
{synopt :{opt pro:file(xvarname)}}PLL is required for the coefficient of variable {it:xvarname}, or{p_end}
{synopt :{opt pro:file([eqname]paramname)}}PLL is required for parameter
{it:paramname} or {opt [}{it:eqname}{opt ]}{it:paramname}{p_end}

{syntab :{it:Syntax 2}}
{synopt :{opt form:ula(formula)}}PLL is required for the unknown parameter in {it:formula}, which is a model term involving at least
one variable in the dataset and one unknown parameter{p_end}
{synopt :{opt pl:aceholder(string)}}sets the placeholder for the unknown parameter in {it:formula} to {it:string}{p_end}

{syntab :{it:Evaluation options: both syntaxes}}
{synopt :{opt dropc:ollinear}}drops any collinear variables{p_end}
{synopt :{opt lev:el(#)}}sets the confidence level to {it:#}{p_end}
{synopt :{opt maxc:ost(#)}}sets an upper limit of 2 * {it:#} on the additional evaluations of the PLL{p_end}
{synopt :{opt n:_eval(#)}}evaluates PLL at {it:#} equally spaced parameter values{p_end}
{synopt :{opt noci}}suppresses calculation of PLL-based CI{p_end}
{synopt :{opt range(#1 #2)}}evaluates PLL over {it:#1} <= parameter <= {it:#2}{p_end}

{syntab :{it:Output options: both syntaxes}}
{synopt :{opt dev:iance}}requests minus 2 times PLL function{p_end}
{synopt :{opt diff:erence}}computes PLL minus maximised log likelihood{p_end}
{synopt :{opt eform}[{opt (string)}]}reports results on the exponentiated scale (with name {it:string}){p_end}
{synopt :{opt gen(newvarlist)}}names the new variables to hold the parameter and the PLL; 
if {cmd:normal()} is specified, also the Normal approximation to the PLL{p_end}
{synopt :{opt nodot:s}}suppresses (supposedly entertaining) dots{p_end}
{synopt :{opt tr:ace}}displays the result of each log-likelihood evaluation{p_end}
{synopt :{opt ver:bose}}displays extended output including results of initial maximum likelihood fitting{p_end}

{syntab :{it:Graph options: both syntaxes}}
{synopt :{opt cilin:es(cline_options)}}specifies rendition of confidence interval{p_end}
{synopt :{opt gropt(cline_opts twoway_opts)}}supplies graph options to enhance PLL plot{p_end}
{synopt :{opt levlin:e(cline_options)}}specifies rendition of horizontal line{p_end}
{synopt :{opt mlel:ine}}adds a horizontal line at the MLE{p_end}
{synopt :{opt nograph}}suppresses the line plot of the results{p_end}
{synopt :{opt norm:al}[{opt (line_options)}]}adds a Normal approximation to the PLL plot.{p_end}
{synoptline}


{pstd}
Any {it:regression_cmd} for which the parameters are estimated
by maximum likelihood may be used. This includes
{help clogit},
{help glm},
{help heckman},
{help intreg},
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
{help stcox},
{help streg},
{help stpm},
{help stpm2}, and probably others.

{pstd}
All weight types supported by {it:regression_cmd} are allowed; see help
{help weight}. All options supported by {it:regression_cmd} should be allowed.


{title:Description}{marker description}


{pstd}
{opt pllf} computes and plots the profile log likelihood function for a
single predictor in a regression model. The term "profile" implies
maximisation over other parameters in the model, if any. It is assumed that
the likelihood function exists.

{pstd}
Plentiful general material on the PLLF is available on the Internet.
This article was inspired by {help pllf##VenzonMoolgavkar1988:Venzon and Moolgavkar (1988)}.

{pstd}
{opt pllf} has two basic syntaxes, depending on which of the options
{opt profile()} or {opt formula()} is used. 

{pstd}
With syntax 1, {opt profile()} must be specified and {opt formula()} is
not allowed. {opt pllf} computes the
profile log likelihood (PLL) function for the regression coefficient
of a covariate defined by {cmd:profile(}{it:xvarname}{cmd:)} or a
parameter defined by {cmd:profile(}[{cmd:[}{it:eqname}{cmd:]}]{it:paramname}{cmd:)}
within a model specified by {it:regression_cmd}.
{opt pllf} also reports
PLL-based confidence limits, computed by a simple grid search.
{it:xvarname} must be a covariate in {it:regression_cmd}. 

{pstd}
With syntax 2, {opt formula()} and {cmd:range()}
must be specified and {opt profile()} is 
not allowed. {opt pllf} computes the PLL function for a non-linear parameter denoted
(by default) by {cmd:@}. {cmd:formula()} is symbolically added to the regression model 
{it:regression_cmd}. {opt pllf} reports the MLE and PLL-based confidence limits,
computed by a simple grid search. Normal-based confidence limits
are not computed. 

{pstd}
With both syntaxes, the results
are saved to new variables assigned by the {cmd:gen()} option.
The dataset length is increased if {cmd:n_eval()} exceeds the number
of observations ({cmd:_N}).


{title:Options}{marker options}

{p 2}{bf:Syntax 1}

{phang}
{cmd:profile(}{it:xvarname} | [{cmd:[}{it:eqname}{cmd:]}]{it:paramname}{cmd:)}
with syntax 1 is not optional. In the first format,
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
{opt range(#1 #2)} with syntax 1 evaluates the PLL function
over {it:#1} <= beta <= {it:#2}, where beta is the 
parameter of interest. Default is for {it:#1} and {it:#2} to be the
confidence limits for beta defined by the option {cmd:level()} and
the usual assumption of a normal distribution for beta_hat, the
maximum likelihood estimate of beta.

{p 2}{bf:Syntax 2}

{phang}
{opt formula(formula)} with syntax 2 is not optional.
{it:formula} defines
a transformation involving at least one variable in the dataset and the
parameter denoted by placeholder {cmd:@}.
{it:formula} may be any valid Stata expression.
Example: {cmd:formula(exp(-@*x5))}.

{phang}
{opt range(#1 #2)} with syntax 2 is not optional.
It evaluates the PLL function over {it:#1} <= {cmd:@} <= {it:#2},
where {cmd:@} is the non-linear parameter of interest. {opt pllf}
also seeks the MLE of {cmd:@}, but if the values of {it:#1 #2} are
ill-chosen or the PLL function behaves 'badly',
it may fail to find the MLE, or give an inaccurate estimate.
The most satisfactory situation is when the MLE lies between {it:#1}
and {it:#2}, and this may be judged from the plot of the PLL function.
In many cases, particularly with large sample sizes,
the PLL function is approximately quadratic with a single maximum.

{phang}
{opt placeholder(string)} defines the placeholder
character(s) used in syntax 2. Spaces or other punctuation
characters are not allowed. Default {it:string} is {cmd:@}.

{p 2}{bf:Evaluation options: both syntaxes}

{phang}
{opt dropcollinear} drops any collinear variables. The default is to 
stop with an error if the x variables are collinear.

{phang}
{opt level(#)} sets the confidence level to {it:#}; default is the current system level,
see {help level}.

{phang}
{opt maxcost(#)} sets an upper limit of 2 * {it:#}
on the number of additional evaluations
of the PLL, when searching for the PLL-based confidence limits. You
should rarely if ever need this option. The program is prevented
by {cmd:maxcost()} from cycling "for ever" when trying to find
confidence limits in pathological cases (see
{opt difference}). Default {it:#} is {cmd:n()}/2.

{phang}
{opt n_eval(#)} evaluates the PLL function at {it:#} equally spaced parameter 
values. Default is 100.

{phang}
{cmd:noci} suppresses calculation of the PLL-based confidence limits.

{p 2}{bf:Output options: both syntaxes}

{phang}
{opt deviance} requests minus 2 times the PLL function, 
i.e. the profile deviance function. If {cmd:difference}
is also specified, {cmd:deviance} produces the profile deviance
difference, i.e. minus 2 times the PLL difference.

{phang}
{opt difference} computes the PLL function minus the maximised
log likelhood for the model. See also {cmd:deviance}. Except
in pathological cases, the resulting values are negative or zero.
Pathological cases denote likelihood functions with multiple maxima
or no maximum.

{phang}{opt eform}[{opt (string)}] reports results for the exponential of the 
coefficient. If {it:string} is specified then this is used as the parameter name.

{phang}
{opt gen(beta_var pll_var [pllnorm_var])} names the new variables created as
{it:beta_var} to contain the values of the regression coefficient
over which the PLL is evaluated, and {it:pll_var} to contain the PLL values, and 
(if {cmd:normal()} is specified) also {it:pllnorm_var} to contain the Normal 
approximation to the PLL values.
If {opt gen()} is not specified, the variables are created
with default names of {opt _beta}, {opt _pll} (and {opt _pllnorm}), respectively.

{phang}
{opt nodots} suppresses dots. By default, a dot is
displayed at each evaluation of the PLL.

{phang}{opt trace} displays the result of each log-likelihood evaluation.

{phang}{opt verbose} displays extended output including results of initial maximum likelihood fitting.

{p 2}{bf:Graph options: both syntaxes}

{phang}
{opt cilines(cline_options)} specifies the rendition of the vertical
lines representing the bounds of the profile-likelihood-based confidence
interval (CI).  See {help cline_options:{it:cline_options}}.

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
{opt levline(cline_options)} specifies the rendition of the horizontal
line showing the profile-likelihood at the confidence level for the the
profile-likelihood-based CI.  See {help cline_options:{it:cline_options}}.

{phang}
{opt mleline} adds a horizontal line at the MLE.

{phang}
{opt nograph} suppresses the line plot of the results.

{phang}
{opt normal}[{opt (line_options)}] adds a Normal approximation to the 
plot of the PLL. {it:line_options} are options valid for {help line}.


{title:Remarks}{marker remarks}

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
many iterations (see {help set iter}), 
and that can take a very long time. If the program appears
to "freeze" in this way, it is best to halt it by pressing ctrl/break. Then
restart it with a more suitable range of values (see the {cmd:range()}
option).

{pstd}
The pseudo standard error of the parameter beta is computed as
upper PLL confidence limit minus lower PLL confidence limit,
divided by twice t, where t is the appropriate quantile
of the t or normal distribution used in calculating normal
based confidence limits. When the sampling distribution of
the parameter of interest is close to normal, the usual standard
error and the PLL-based standard error will be approximately equal.

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

{title:Limitations}{marker limitations}

{pstd}
{cmd:pllf} only handles models with factor variables and interactions 
if the prefix command {help xi} is used,
and not if the {help fvarlist} syntax is used.

{pstd}
Profile likelihood calculation is a likelihood-based method, so it 
does not respect non-likelihood-based variance calculations such as
the robust (sandwich) or cluster-robust variance. 

{pstd}
{cmd:pllf} is unable to profile over a variance component.

{pstd}
Some commands may not be supported. These include {help reg3} and (at the time of writing)
{help stpm3}.


{title:Examples}{marker examples}

{p 2}{bf:Syntax 1}

{phang}Input two-group data with event outcome

{phang}. {stata "clear"}{p_end}
{phang}. {stata "input group pyears events"}{p_end}
{phang}{space 4}{stata "      0 200 8"}{p_end}
{phang}{space 4}{stata "      1 100 3"}{p_end}
{phang}{space 4}{stata "      end"}{p_end}

{phang}Fit Poisson model

{phang}. {stata "poisson events group, exposure(pyears)"}{p_end}

{phang}Explore profile likelihood for the coefficient of group

{phang}. {stata "pllf, profile(group): poisson events group, exposure(pyears)"}

{phang}Load TRISST trial data ({help pllf##Joffe2022:Joffe et al, 2022}). The data 
are from a non-inferiority trial of MRI vs CT for surveillance after testicular 
cancer. The PLL CI does not cross zero while the Normal CI does. However this 
is a non-inferiority trial with margin a risk difference of +0.057, so both CIs 
clearly establish non-inferiority.

{phang}. {stata "use TRISST, clear"}{p_end}
{phang}. {stata "pllf, norm verbose profile(modality): binreg outcome modality [fw=n], rd"}

{phang}Load breast cancer data

{phang}. {stata "webuse brcancer, clear"}{p_end}
{phang}. {stata "stset rectime, failure(censrec) scale(365.24)"}

{phang}Explore profile likelihood for coefficient of x5e

{phang}. {stata "pllf, profile(x5e) range(-3 -1): stcox x1 x4a x5e x6 hormon"}

{phang}Explore profile likelihood for coefficient of x1

{phang}. {stata "pllf, profile(x1) gen(X Y): stcox x1 x4a x5e x6 hormon"}

{phang}Explore profile likelihood for the constant of Weibull distribution

{phang}. {stata "pllf, profile(_cons) n(50): streg x1 x4a x5e x6 hormon, distribution(weibull)"}

{phang}Explore profile likelihood for the constant of flexible parametric model (stpm2)

{phang}. {stata "pllf, profile(_cons) n(50): stpm2 x1 x4a x5e x6 hormon, df(2) scale(hazard) "}

{phang}Explore profile likelihood for Weibull shape parameter

{phang}. {stata "pllf, profile([ln_p]_cons) n(50): streg x1 x4a x5e x6 hormon, distribution(weibull)"}

{phang}Explore profile deviance for predictor of Weibull shape parameter

{phang}. {stata "pllf, profile([ln_p]x4b) deviance difference n(20): streg x1 x4a x5e x6 hormon, distribution(weibull) ancillary(x4b) "}

{p 2}{bf:Syntax 2}

{pstd}The following two commands with syntax 2 are equivalent and give
identical results. The {cmd:stcox} models are nonlinear in the parameter @ but conditionally linear
in exp(-@*x5).

{phang}. {stata "pllf, formula(exp(-@*x5)) range(.05 .25): stcox x1 x4a x6 hormon"}{p_end}
{phang}. {stata "pllf, placeholder(X) formula(exp(-X*x5)) range(.05 .25): stcox x1 x4a x6 hormon"}

{pstd}
In both commands, values of {it:beta} in {cmd:range(.05 .25)} are substituted in turn into
the extra model term exp(-{it:beta}*x5), and the corresponding log likelihood values are determined.


{title:Stored results}{marker stored_results}

{pstd}
{opt pllf} stores the following scalars in {cmd:r()}:

{synoptset 23 tabbed}{...}
{synopt:{cmd:r(nobs)}}      number of observations in the estimation sample{p_end}
{synopt:{cmd:r(b)}}         MLE of parameter{p_end}
{synopt:{cmd:r(se)}}        usual standard error of parameter (syntax 1 only){p_end}
{synopt:{cmd:r(pse)}}       pseudo standard error of parameter{p_end}
{synopt:{cmd:r(n_llci)}}    lower normal-based confidence limit for parameter (syntax 1 only){p_end}
{synopt:{cmd:r(n_ulci)}}    upper normal-based confidence limit for parameter (syntax 1 only){p_end}
{synopt:{cmd:r(l_llci)}}    lower PLL-based confidence limit for parameter{p_end}
{synopt:{cmd:r(l_ulci)}}    upper PLL-based confidence limit for parameter{p_end}
{synopt:{cmd:r(ll)}}        maximised log likelihood{p_end}
{synopt:{cmd:r(ll_limit)}}  value of log likelihood for likelihood based CI calculation{p_end}
{synopt:{cmd:r(cost)}}      "cost" of evaluating PLL-based confidence limits{p_end}
{synopt:{cmd:r(asym)}}      estimated asymmetry of PLL function{p_end}


{title:References}{marker references}

{phang}Please use this reference to cite this software:

{phang}
P. Royston. 2007. Profile Likelihood for Estimation and Confidence Intervals. The Stata 
Journal 7, 376–387. {browse "https://doi.org/10.1177/1536867X0700700305"}.

{phang}Other references:

{phang}{marker Joffe2022}
J. K. Joffe et al. 2022. Imaging Modality and Frequency in Surveillance of Stage I 
Seminoma Testicular Cancer: Results From a Randomized, Phase III, Noninferiority 
Trial (TRISST). Journal of Clinical Oncology 40, 2468–2478. 
{browse "https://doi.org/10.1200/JCO.21.01199"}.

{phang}{marker VenzonMoolgavkar1988}
D. J. Venzon and S. H. Moolgavkar. 1988. A method for computing
profile-likelihood-based confidence intervals. Applied Statistics
37: 87-94.
{browse "https://www.jstor.org/stable/10.2307/2347496"}


{title:Authors}{marker authors}

{pstd}Patrick Royston, MRC Clinical Trials Unit at UCL, London, UK.{break}

{pstd}Ian White, MRC Clinical Trials Unit at UCL, London, UK.{break}
Email: {browse "mailto:ian.white@ucl.ac.uk":Ian White}


{title:Also see}
