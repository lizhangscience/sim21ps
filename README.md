# sim21ps
***
In order to precisely detect and analyze the 21cm HI signal, which performs as a good probe to reveal the epoch of reionization (EoR), a simulation is conducted. Among the foregrounds whose  brightness are orders higher than the target 21cm signal, the extragalactic point sources (PS)  are my field of interest. Those point or point-like sources are usually glaxies, quasars, blazers and other glaxy-like sources, which diffuse bright light among different radio bands.

Usually, there is a large massive black hole in the center of a point source. Since the black hole fullfills itself by absorbing high energy dust or gas around it, an accreation disk is generated around the black hole. However, the absorbed energy will be ejected by the black hole as a form of two anti-direction jets (e.g. M87), which are usually vertical to the disk. This kind of strategy is called active galaxy nuclei (AGN).

![M87](https://upload.wikimedia.org/wikipedia/commons/thumb/0/07/Messier_87_Hubble_WikiSky.jpg/250px-Messier_87_Hubble_WikiSky.jpg)

Works to classify the point sources according to their spectral and spatial features were conducted. The PSs are usually categorized into star forming (SF), star bursting (SB), radio-quiet (RQ) and the Fairnoff-Railey I and II AGN.

|ClassType| SF | SB |RQ AGN|FRI |FRII|
|---------|----|----|------|----|----|
|  Code   | 1  |  2 |   3  | 4  | 5  |

In this work, I am going to simulate the point sources as well as the jets, so as to mimid them precisely.Codes of the five types PSs are listed above.

## Point sources simulation
The folder namely sim21ps holds python3 scripts is to simualte the AGN like point sources in the radio sky map. In order to successfully take advantage of this code, some python modules should be installed ahead. They are listed as follows, and the version of `Python` should be `3.4.x` or above.

- astropy or pyfits : processing on the fits files
- numpy,scipy: : processing with ndarray formation data under numpy
- matplotlib,PIL : displaying results and write or save the images
- pandans : processing on the point sources catalog in `csv` formation.
- fg21sim: a three party package designed by [Weitian Li](https://github.com/liweitianux/fg21sim)

### How to use it?
There are three ways to use my package, (1) install it in your computer (TODO), (2) import the modules when simulating PS, and (3) automatically generate PS catelogue and fits files in the terminal.However, at present,the simulation scripts can only worked under the `linux` system.

1. Install the pakcage
   <TODO>

2. Import as a module
	- `basid_params`
a module to generate cosmological constants and calute objects scales at the provided redshift (z).
	- `psCatelogue`
a module to generate PS lists or catelogue, the PS type and amount can be mannually settled.
	- `psDraw`
a module to fill the simulated PS with surface britness at provided frequency, and ourput the simulated all sky map.

3. Automaticaaly generation
	- To generate PS catelogues
	````sh
	python3 ./psCatelogue.py -c <ClassType> -n <NumPS> -o <Output file>
	```
	- To generate the fits file
	````sh
	python3 ./psDraw.py -i <PS_catelogue(csv)> -o <Output fits> -n <nside> -f <frequency> 
	````

## Author
- Zhixian MA <`zxma_sjtu(at)qq.com`>

## License
Unless otherwise declared:

- Codes developed are distributed under the [MIT license](https://opensource.org/licenses/mit-license.php);
- Documentations and products generated are distributed under the [Creative Commons Attribution 3.0 license](https://creativecommons.org/licenses/by/3.0/us/deed.en_US);
- Third-party codes and products used are distributed under their own licenses.

