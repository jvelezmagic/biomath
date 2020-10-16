# Biomath in times of Coronavirus

Since the outbreak of the coronavirus disease 19 (COVID-19), many mathematical models have been proposed, various of which several governments had used to implement their policy to lockdown the population. In particular, Mexico City has posed for an ordinary differential equations model to keep the spread of the virus under control. But, no optimal span of the lockdown was proposed, and therefore, the city is in an indefinite quarantine, affecting the economy of the whole country. Here we assess different aspects of the coronavirus epidemic in Mexico City, studying the dynamics in response to different spans and intensities of the lockdown, performing simulations using the ordinary differenltial equations model proposed by the government of Mexico City, and a new individual agent-based model that takes into account hygiene measures of individuals.

## Examples of models
### ABM - Periodic lockdown with vaccination

![Periodic lockdown with vaccination](./plots/social_distancing/periodic_lockdown_with_vaccination_slow_β=50-50.gif)
![Example of SIRV results](plots/social_distancing/fraction_isolated=0.5_interval_between_locks=336_nsteps=12000_periodic_lockdown=true_start_lockdown=720_start_vaccination=7200_stop_lockdown_after=720_βmax=4_βmin=0.png)

### ODE - Shield Immunity

#### Scheme
![](./plots/shield_immunity/scheme.png)

#### Epidemic differences depending on the severity of sanitary measures.
![](./plots/shield_immunity/covid_epidemic_US-UK.jpg.png)
