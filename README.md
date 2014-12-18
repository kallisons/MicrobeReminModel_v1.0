MicrobeReminModel_v1.0
======================
The Microbial Remineralization Model v1.0 simulates the interactions between sinking particles and heterotrophic bacteria in the ocean water column in a 1-dimensional Eulerian framework. The model has 9 state variables including particulate organic carbon, particle-attached bacteria, free-living bacteria, active exoenzyme in the particle, inactive exoenzyme in the particle, hydrolysate in the particle, hydrolysate in the dissolved environment, active exoenzyme in the dissolved environment, and inactive exoenzyme in the dissolved environment.

Please cite the following paper if you use this code:

Mislan KAS, CA Stock, JP Dunne, and JL Sarmiento. 2015. Group behavior among model bacteria influences particulate carbon remineralization depths.  Journal of Marine Research.

----------------------
Software dependencies
----------------------
gfortran version 5.0.0:   [https://gcc.gnu.org/wiki/GFortran](https://gcc.gnu.org/wiki/GFortran)

R version 3.1.2: [http://www.r-project.org/](http://www.r-project.org/)

---------
Folders
---------
**Data** -  The model uses data from the Bermuda Atlantic Time Series (BATS).  Data downloaded from the [BATS](http://bats.bios.edu/) website is in this folder.  

**Graphs** -  Graphs of the model output plotted using R are saved to this folder.  

**Input** -  The 1-dimensional (MRM1D) configuration of the model requires two input files with environmental data. One input file has particulate organic carbon flux at the base of the euphotic zone for each time step. Another input file has a long-term average vertical profiles of temperature, salinity, density, and semi-labile dissolved organic carbon.  

**Model** -  This folder contains 0-dimensional (MRM0D) and 1-dimensional (MRM1D) configurations of the Microbial Remineralization Model v1.0.

**Output** -  The shell script that runs the model will put output files from the model in this folder.  The R code will access the model output from this folder.  

**RCode** - The R code files to format the input data (DatatoInput folder) and analyze output (Analysis folder) and generate figures (Figures folder) are contained in this folder.

**TestFiles** - Contains test files to verify the example output was generated correctly by the model and R scripts.  

---------------------
Compiling the model
---------------------

**0-dimensional configuration:**  
Open a shell window and change the directory to the MRM0D folder inside the MicrobeReminModel_v1.0/Model folder.  

Change directory:  

    cd MicrobeReminModel_v1.0/Model/MRM0D

Command to compile model:

    gfortran MRM0D_v1.0.F90 -o MRM0D1 -ffpe-summary=invalid,zero,overflow

**1-dimensional configuration:**  
Open a shell window and change the directory to the MRM1D folder inside the MicrobeReminModel_v1.0/Model folder.  

Change directory:  

    cd MicrobeReminModel_v1.0/Model/MRM1D

Command to compile model:

    gfortran MRM1D_v1.0.F90 -o MRM1D1 -ffpe-summary=invalid,zero,overflow


------------------
Running the model
------------------
**0-dimensional configuration:**

Change directory:  

    cd MicrobeReminModel_v1.0/Model/MRM0D

Command to run model:

    sh RunModel_BATS_MRM0D.sh

**1-dimensional configuration:**

Change directory:  

    cd MicrobeReminModel_v1.0/Model/MRM1D

Command to run model:

    sh RunModel_BATS_MRM1D.sh


----------------------------
Test model output
----------------------------
Compare example model output created using the commands above to a set of test files to make sure the results are the same.  

**0-dimensional configuration:**

Command to run comparison tests:

    sh RunTest_BATS_MRM0D_Output.sh

**1-dimensional configuration:**

Command to run comparison tests:

    sh RunTest_BATS_MRM1D_Output.sh

--------------------------------
Exoenzyme optimization analysis
--------------------------------


-------------------
Changing scenarios
-------------------
The interior, interception, and retention scenarios.


------------------
Acknowledgements
------------------

Code Release:  
KAS was supported by the Gordon and Betty Moore Foundation, the Alfred P. Sloan Foundation, and the Washington Research Foundation through the eScience Institute at the University of Washington.

Scientific Research and Code Development:  
KAS was supported by the Carbon Mitigation Initiative at Princeton University which is sponsored by BP and the NOAA Cooperative Institute for Climate Science (NA08OAR4320752).
