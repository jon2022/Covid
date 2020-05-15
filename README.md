# Covid simulation model
This is a model of infection spread in Covid-19 based upon a simulation using a SEIR structure.  My concern is that many models use an average population R0, which can be altered by various measures to suppress the spread of the virus.  However, the evidence is that suppression measures don't have an immediate effect, but the R0 appears to reduce gradually.

My suggestion is that there are two particular effects that are worth modelling.  The first is that there is a subgroup of population (for example, key workers and those requiring personal care), who are unable to suppress to the same extent as the rest of the population, and have a partially reduced R0.

The second issue that the model addresses is that 'lockdown', at least in most countries, occurs in households rather than individually, so that the R0 is unlikely to reduce within the household group (in fact it may increase within a locked down household).

The model considers four different R0 values that apply to different portions of the population, the general value for an unsuppressed population, the value within a locked down household and the values for  partially and fully suppressed populations.  The household size, proportion of population partially and fully suppressed and values such as the number of seeding cases and duration of incubation and infectivity can be set.

These are potentially important issues because different suppression measures may affect different aspects of the spread.  For example, allowing wider family groups to combine may increase the effective household size, allowing greater freedom of movement or intensive contact tracing may increase the general R0 of the supressed population.  Re-opening certain industries may increase the size of the population that are partially suppressed and changes to PPE availability may alter the R0 of the partially suppressed population.

This is a model in development and comments, corrections and suggestions would be welcome.

There are two files - the first is the model itself "SEIR_InitialModel.R" and "example_simulations.R", which has an example of testing a set of parameters and providing graphical comparisons of the output.
