# AdaptiveAntivenomDesign

## Overview

This repository has a set of R scripts which simulate antivenom dose-finding trials under adaptive designs.
The main design of interest is a model-based adaptive design. For comparison, we also have coded up a modified `3+' design.


## Background

The clinical context of many antivenoms is as follows: 

* A binary toxicity outcome of interest: e.g. occurence of anaphylactic shock shortly after administration of the antivenom (usually defined as within 180 minutes of administration)

* A binary efficacy outcome of interest: e.g. restoration of blood coagulability (for Russell's viper envenoming, this is usually defined by the 20 minute whole blood clotting test 6 hours after antivenom administration)

Therefore a desirable antivenom dose should have a low probability of causing toxicity, and a high probability of being efficacious. If we state a maximum tolerated toxicity (MTT) and a target level efficacy (TEL), then there is a dose which has as average toxicity the MTT: this is the maximum tolerated dose (MTD); and there is a dose which has as average efficacy the TEL: this is the target efficacious dose (TED). The optimal dose can be defined as the lower between the MTD and TED.


## Details of the simulation

The workhorse scripts can be found in the file *Simulation_functions.R*


## Output of main scripts

The main RMarkdown file runs 7 simulation scenarios as detailed in the paper - see preprint at XXX.

