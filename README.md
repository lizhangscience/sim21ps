# sim21ps
==========
In order to precisely detect and analysis the 21cm HI signal, which performs a good probe to reveal the epoch of reionization (EoR), a simulation is conducted. Among the foreground whose temperature brightness are orders higher than the target 21cm signal, the extragalactic point sources are my field of interest. The point sources are usually glaxies, quasars, blazers and other glaxy-like sources, which diffuse bright light among different radio bands. 

Usually, there is a large massive black hole in the center of a point source. Since the black hole fullfills itself by absorb high energy dust or gas around it, an accreation disk is generated around the black hole. Also, these absorbed energy will be ejected by the black hole as two anti-direction jets (e.g. M87). I am going to simulate the point sources, as well as the jets, so as to mimid them precisely.

![M87](https://upload.wikimedia.org/wikipedia/commons/thumb/0/07/Messier_87_Hubble_WikiSky.jpg/250px-Messier_87_Hubble_WikiSky.jpg)

## Point sources simulation
The script namely SimAGNReal.py is to simualte AGN like point sources in the radio sky map. In order to successfully take advantage of this code, some python modules should be installed ahead. They are listed as follows,

1. pyfits
2. numpy,scipy
3. matplotlib

### How to use it?
    import SimAGNReal
    SimAGNReal.GenMultiFRs(Rows,Cols,Freq,NumFR)
    

