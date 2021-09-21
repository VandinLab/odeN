# Code for review of the paper: odeN Simultaneous Approximation of Multiple Motif Counts in Large Networks (CIKM 2021)

In this repository you will find the code used as proof of concepts for the paper odeN: Simultaneous Approximation of Multiple Motif Counts in Large Temporal Networks. 
Additionally we report a file with the parameters used in our experimental evaluation in the folder parameters/parameters.pdf.
If you use this code for your purpose, please cite our work!

## Building the code
NOTE: only Ubuntu >=18.04 with g++9 is supported!
Building the code of odeN requires various steps, which are now described:

- First you will need to compile and install [Lemon](https://lemon.cs.elte.hu/trac/lemon) if you want to execute the algorthm based on the VF2++ algorithm. We provide the version we used of lemon in the folder lemon-main, to build it and install it correctly follow the instruction on [this](https://lemon.cs.elte.hu/trac/lemon/wiki/InstallLinux) page, that we will also summarize here
    - ```cd lemon-main```
	- ```mkdir build```
    - ```cd build```
    - ```cmake ..```
    - ```make```
    - Now you have to install it in order to use the API in our progaram, to do so we suggest to execute the following line ```cmake -DCMAKE_INSTALL_PREFIX=<path> ..``` where ```<path>``` it the path to the folder ```lemon``` in our main folder (where this README is placed). Then execute ```make install```
- Now you should build SNAP, to do so just run ```make``` in the current folder.
- Once SNAP is built you can go to the folder ```examples/staticsampler``` and edit the file ```../../Makefile.config``` by setting the version of the compiler to ```g++2a``` then run ```make```. This should build correctly the code

Now you can run the executable ```./staticsamplermain``` to print a list of the parameters. We report also in ```example-usage.txt``` some example of how to run it properly.

## Known Issues and notes
There are some known compilation issues due to the SNAP library that is required in the project, depending on the OS version where you compile the whole package, in particular:
* If you run the `make` command on a Ubuntu OS with a version lower than 18.04 (i.e., 14.04, 16.04) then you need to comment the struct `__exception` in lines 17-23 of the file `glib-core/bd.h`

We will soon update the repository by rendering the dependecy on Lemon optional and by providing better documentation on odeN functionalities. Stay tuned!
