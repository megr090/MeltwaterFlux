# MeltwaterFlux
This repository contains code that corresponds to the recently submitted paper "Meltwater generation in ice stream shear margins: case study in Antarctic ice streams".

findIceTemperature.m computes ice temperature from inputted Brinkman number, whereas findIceTemperaturefromSR.m computes ice temperature from inputted strain rate. Both are based on the model presented in Meyer and Minchew 2018.

findOuterPressureProfile.m, findInnerPressureProfile.m, and findCompositePressureProfile.m find effective pressure from inputted ice temperature and select variables described in the paper.

findOuterPorosityProfile.m, findInnerPorosity0Profile.m, findInnerPorosity1Profile.m, findInnerPorosityProfile.m, and findCompositePorosityProfile.m find ice porosity from ice temperature and effective pressure, along with variables outlined in the paper.

findMeltwaterFlux.m computes meltwater flux with depth from ice temperature, porosity, and effective pressure.

Run main_meltwaterflux_comparenumerics.m to generate Figure 1 in the submitted paper, which compares the model described above to results from a numerical model presented in Meyer et al. 2018. 

Run main_meltwaterflux_margins_downstream.m to generate the plots in Figure 4 of the submitted paper, which computes porosity, effectiv epressure, and meltwater flux with depth for 3 ice streams: Bindschadler Ice Stream, Byrd Glacier, and Pine Island Glacier.
