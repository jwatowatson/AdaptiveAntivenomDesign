# AdaptiveAntivenomDesign

## Overview

This repository has a set of R scripts which simulate antivenom dose-finding trials under adaptive designs.
The main design of interest is a model-based adaptive design. For comparison, we also have coded up a modified `3+' design.


## Background

The clinical context of many antivenoms is as follows: 

* A binary toxicity outcome of interest: occurence of anaphylactic shock shortly after administration of the antivenom (usually defined as within 180 minutes of administration)

* A binary efficacy outcome of interest: restoration of blood coagulability (for Russell's viper envenoming, this is usually defined by the 20 minute whole blood clotting test 6 hours after antivenom administration)

Therefore a desirable antivenom dose should have a low probability of causing the toxicity outcome, and a high probability of causing the efficacy outcome. If we state a maximum tolerated toxicity (MTT) and a target level efficacy (TEL), then there is a dose which has as average toxicity the MTT: this is the maximum tolerated dose (MTD); and there is a dose which has as average efficacy the TEL: this is the target efficacious dose (TED). The optimal dose can be defined as that which is the lower of the MTD and TED:
$V_{\text{star}} = \min(MTD, TED)$


## Details of the adaptive designs

TODO

## Output of main scripts

The main RMarkdown file runs 4 simulation scenarios.
