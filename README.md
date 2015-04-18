[![DOI](https://zenodo.org/badge/8411/kallisons/MicrobeReminModel_v1.0.svg)](http://dx.doi.org/10.5281/zenodo.16145)

Microbial Remineralization Model v1.0
======================
The Microbial Remineralization Model v1.0 simulates the interactions between sinking particles and heterotrophic bacteria in the ocean water column in a 1-dimensional Eulerian framework. The model has 9 state variables including particulate organic carbon, particle-attached bacteria, free-living bacteria, active exoenzyme in the particle, inactive exoenzyme in the particle, hydrolysate in the particle, hydrolysate in the dissolved environment, active exoenzyme in the dissolved environment, and inactive exoenzyme in the dissolved environment.

Please cite the following paper if you use this code:

Mislan KAS, CA Stock, JP Dunne, and JL Sarmiento. 2014. Group behavior among model bacteria influences particulate carbon remineralization depths.  Journal of Marine Research. 72:183-218

----------------------
Software dependencies
----------------------
All the required software is open source:

gfortran version 5.0.0:   [https://gcc.gnu.org/wiki/GFortran](https://gcc.gnu.org/wiki/GFortran)

R version 3.1.2: [http://www.r-project.org/](http://www.r-project.org/)

Mac OS X and Unix-like operating systems should be able to install gfortran and R without any additional dependencies.  Information on running the model on computers with Microsoft Windows operating systems is located towards the end of this README file.

---------------------
Configurations
---------------------
**0-dimensional configuration:**  The 0D configuration of the Microbial Remineralization Model was developed to determine the effect of exoenzyme production on the growth of bacteria on particles.  This configuration was derived from the 1D configuration and used to determine optimal exoenzyme production.

**1-dimensional configuration:**  The 1D configuration of the Microbial Remineralization Model was developed to determine the effect of particle-attached bacteria on particle decomposition and remineralization in the deep ocean water column.  The depth range of the model is from 150 to 4000 m.  

---------
Folders
---------
**Analysis** - Results from the optimal exoenzyme production analysis are saved to this folder.

**Graphs** -  Graphs of the model output plotted using R are saved to this folder.

**Input** -  The 1-dimensional (MRM1D) configuration of the model requires two input files with environmental data. One input file has particulate organic carbon flux at the base of the euphotic zone for each time step. Another input file has a long-term average vertical profiles of temperature, salinity, density, and semi-labile dissolved organic carbon.  

**Model** -  This folder contains 0-dimensional (MRM0D) and 1-dimensional (MRM1D) configurations of the Microbial Remineralization Model v1.0.

**Output** -  The shell script that runs the model will put output files from the model in this folder.  The R code will access the model output from this folder.  

**RCode** - The R code files are contained in this folder.

**TestFiles** - Contains test files to verify the example output was generated correctly by the model.  

---------------------
Compiling the model
---------------------

**0-dimensional configuration:**  
Open a shell window and change the directory to the MRM0D folder inside the MicrobeReminModel_v1.0/Model folder.  

Change directory:  

    cd MicrobeReminModel_v1.0/Model/MRM0D

Command to compile model:

    gfortran MRM0D_v1.0.F90 -o MRM0D -ffpe-summary=invalid,zero,overflow

**1-dimensional configuration:**  
Open a shell window and change the directory to the MRM1D folder inside the MicrobeReminModel_v1.0/Model folder.  

Change directory:  

    cd MicrobeReminModel_v1.0/Model/MRM1D

Command to compile model:

    gfortran MRM1D_v1.0.F90 -o MRM1D -ffpe-summary=invalid,zero,overflow

------------------
Running the model
------------------
The model is run using shell scripts.

**0-dimensional configuration:**

Change directory:  

    cd MicrobeReminModel_v1.0/Model/MRM0D

Command to run model:

    sh RunModel_BATS_MRM0D.sh

**Important Note:**  The resolution of max_epsilon values in the RunModel_BATS_MRM0D.sh file was reduced decrease the time required to perform functionality tests of the model output. The resolution can be changed in the RunModel_BATS_MRM0D.sh file.  

**1-dimensional configuration:**

Change directory:  

    cd MicrobeReminModel_v1.0/Model/MRM1D

Command to run model:

    sh RunModel_BATS_MRM1D.sh

**Important Note:**  MRM1D is set to run for 15 days for the tests of the model output to decrease the time required to perform functionality tests.  The model needs to be run for at least 180 days to reach steady state.  The number of time intervals can be changed in the PRINT FINAL OUTPUT TO FILES section towards the bottom of the MRM1D_v1.0.F90 file.

----------------------------
Test model output
----------------------------
Compare example model output created using the commands above to a set of test files to make sure the results are the same.

**0-dimensional configuration:**

Command to run comparison tests:

    RScript RunTest_BATS_MRM0D_Output.R

**1-dimensional configuration:**

Command to run comparison tests:

    RScript RunTest_BATS_MRM1D_Output.R

----------------------------------------------
Optimal exoenzyme production analysis MRM0D
----------------------------------------------
Change directory:  

    cd MicrobeReminModel_v1.0/RCode

Calculate and plot the optimal enzyme production rate (h<sup>-1</sup>) for different initial bacterial abundances per particle using R:

    Rscript MRM0D_OptimalEpsilon.R

The optimal epsilon values are written to text file in the MicrobeReminModel_v1.0/Analysis folder.  A plot of optimal epsilon is saved to the MicrobeReminModel_v1.0/Graphs folder.  

Command to run comparison tests:

    RScript RunTest_BATS_MRM0D_OptimalEpsilon.R

**Important Note:**  The resolution of max_epsilon values in the RunModel_BATS_MRM0D.sh file was reduced to decrease the time required to perform functionality tests of the model output.  The resolution can be changed in the RunModel_BATS_MRM0D.sh file.

------------------------
Plot output from MRM1D
------------------------
Change directory:  

    cd MicrobeReminModel_v1.0/RCode

Plot State Variables using R:

    Rscript MRM1D_StateVariables.R

Plot Bacterial Rates using R:

    Rscript MRM1D_BacteriaRates.R

Plot Mass Transfer Rates using R:

    Rscript MRM1D_MassTransferRates.R

Plots are saved as postscript files in the MicrobeReminModel_v1.0/Graphs folder.

-------------------
Changing scenarios
-------------------
Three scenarios were explored:

(1) **Interior** - bacterial uptake is separate from diffusive flux from particles  
(2) **Interception** - bacterial uptake is from the diffusive flux from particles  
(3) **Retention** - exoenzyme and hydrolysate flux is stopped by particle-attached bacteria

The default setting of the Microbial Remineralization Model (0D and 1D configurations) for the tests of the model output is Interception (2).  The scenario setting can be changed in the SCENARIO sections in the MRM0D_v1.0.F90 and MRM1D_v1.0.F90 files.

-------------------
Input files
-------------------
**Forcing File**  
The model forcing input file for MRM1D has no header and contains the following comma delimited columns:

* day
* hour and minute in the hhmm format
* seconds
* carbon flux at 150 m (mg m<sup>-2</sup> h<sup>-1</sup>)

The number of time steps in the forcing input file needs to be less than or equal to the number of time steps in the entire time period being simulated.

**Depth Profile**
MRM1D requires a depth profile file with the following space delimited columns at 10 m depth intervals:

* depth (m)
* temperature (degrees C)
* salinity (psu)
* density (kg m<sup>-3</sup>)
* semi-labile doc (mg m<sup>-3</sup>)

--------------------------------
Microsoft Windows
--------------------------------
Windows operating systems require a unix environment to be installed in order for gfortran to be installed.  Options include:

Cygwin: [www.cygwin.com](www.cygwin.com)  
MinGW: [www.mingw.org](www.mingw.org)

In order to use the Rscript command in Microsoft Windows, the file path to the Rscript.exe needs to be specified.

Command to run R code from a Windows CMD window:

c:\progra~1\R\R-3.1.2\bin\Rscript.exe MRM1D_StateVariables.R

Command to run R code from a Windows cygwin shell:

/cygdrive/c/progra~1/R/R-3.1.2/bin/Rscript.exe MRM1D_StateVariables.R

------------------
Acknowledgements
------------------
Ben Marwick, University of Washington, and Rahul Biswas, University of Washington, vetted this code release.

Code Release:  
KAS was supported by the Gordon and Betty Moore Foundation, the Alfred P. Sloan Foundation, and the Washington Research Foundation through the eScience Institute at the University of Washington.

Scientific Research and Code Development:  
The project was supported by the Carbon Mitigation Initiative at Princeton University which is sponsored by BP and the NOAA Cooperative Institute for Climate Science (NA08OAR4320752).
