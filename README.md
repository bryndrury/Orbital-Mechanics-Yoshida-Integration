# Orbital-Mechanics-Yoshida-Integration

**USE AT YOUR OWN RISK**
(idk what risk there is but if it goes t*ts-up don't blame me)

Simulating orbits using Newtonian mechanics and the Yoshida integrations.

It doesnt use optimazation algorithums like the Barnes-Hut approximation algorithum (I'm woring on that) so it's just an N-body problem 'solution'.

The Computation part of the code is done in C++ but the analysis, (i.e. the graphs and plots) are done in matplotlib in python.

If you run the sim.py, it 'should' compile and run the code with some included intial conditions. 

Some initial conditions are include as csv files (sun_earth.csv, inner_planets.csv and solar_system.csv) but as long as the csv file is in the right format you should be good to put what starting conditions you like. Just run the code from command line as shown in the sim.py file.

This was an adaptation of a project that I had for an uni project, but I decided to make a few changes myself. 
(I've been working on a render engine, so I might add live graphics if I get around to it.)
