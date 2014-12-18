Action-Potential-Refractory-Period-Model
========================================
This is an action potential model written in MATLAB as a project for one of my classes, BME 4641, Bioelectricity.  Included in the code is a refractory period model based upon parameters and differential equations outlined in the Hodgkin and Huxley model of the neuron (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/), modeling the elicited amplitude response of the AP at certain restimulus times.  The simultaneous differential equations are solved with a numerical Runge-Kutta method.

actionPotential.m contains the source code for generating the action potential with conductance constants and ionic currents.

refractory.m contains source code for modeling the AP response in the relative and absolute refractory period of the action potential.
