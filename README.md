# Introduction:
5GPencil Simulator is able to synthesize the pencil beams for each gNB and to evaluate the overall impact in terms of EMF exposure and UE throughput levels over the territory under consideration.


## Requirements:
- Matlab: 
Minium version tried MATLAB Version: 9.5.0.1298439 (R2018b) Update 7 Also tried on MATLAB Version: 9.9.0.1495850 (R2020b) Update 1
- Matlab Toolbox List:
	* Simulink
	* DSP System Toolbox
	* Mapping Toolbox
	* Optimization Toolbox
	* Phased Array System Toolbox
	* Signal Processing Toolbox
	* Statistics and Machine Learning Toolbox
	* Symbolic Math Toolbox


## Project Structure:
The project is composed by the following files
- grid_generator.m
The goal of this module is to define the cells layout and to generate 
	* the 'Deployment spots' (points where the beam will be directed) and
	* the 'Measure spots' (points that measure the power density and the electric field)
- simulator_5GPencil.m
This module is the core of the project that is the simulator. 
The simulator has 4 goals: 
	1.  Synthesizing the beams, aimed at the deployment spots defined in the previous module; 
	2.  Generating each user within a circular area of diameter epsilon (localization uncertainty) inside each deployment spot;
	3. Computing SINR and throughput over users and deployment spots
	4.  Computing the contribution in terms of electromagnetic field over the serving cell, which is the central one surrounded by 6 cells.


## Execution:
1. Run grid_generator.m with `measure=1;`
**NOTE:**  default radius is set to 100 meters, and it represents the radius of the cell, you can change it by changing `r `variable.
In this step you are generating the measurement spots in the whole central cell. In order to limit the computational complexity, the plots are deactivated by default, but you can activate them by imposing `plot_on=1;`
**NOTE:**   whenever you change cell radius you have to repeat step 1. 
2. Run grid_generator.m with `measure=0;`
Now you generate the deployment spot. 
At the end of these two steps, 4 files are generated (needed for next step):
	-  layout.mat: contains the basic configuration: gNBs, cells and sectors
	-  useful_spot.mat:  contains all points
	-  measure.mat: containing all measurement spots
	-  deploy.mat: containing all deploy spots

3. Run simulator_5GPencil.m
Default setting:
    ``` matlab
    measure_on=1; %1 to enable computation of emf, 0 to turn off
    fixed_emf = 0; %1 to enable computation of emf when no beamforming is used
    plot_on=0; %1 to enable plot, 0 to disable
    NLOS_on=1; %1 for NLOS, 0 for LOS
    int_intrasec=1; %1 to consider intrasector interference, 0 to disable
    fig_style=1; %to enable plot of user tolerance areas
    accuracy=2; %diameter (m) of area around deploy spot in which the UE is deployed
    ```
    **NOTE:**  In lines 178 and 179 you can choose to set az3dB equal to 65° or 90° respectively and in lines 209 and 210 you can choose to set el3dB equal to 65° or 7° respectively.

At the end of these steps, the following files are generated:
- deploy_comp.mat: all the deployment spots with their infos
- user_comp.mat: all the users with their infos
- fixed_emf_data.mat 
- angle_3dB.mat: the azimuth and elevation 3dB angles for each beams of each sector
- emf_data.mat: all measure spots with their infos

Moreover there is a directory "Data" in which the following files will be saved: `emf<accuracy>.mat` (for example: emf2.mat or emf16.mat) this file contains the avg emf, the emf confidence interval, and the emf.

**NOTE for macOS users:**  all backslashes should be changed with forward slashes, because of Linux-friendly notation. For example line 649 should be changed as follows:
`fname=join(['Data/emf',num2str(accuracy),'.mat']);`


### 5GPENCIL Team

Main contributors:
- Simone Rossetti <simone.rossetti@cnit.it>
- Sara Saida <sara.saida@cnit.it>

Other contributors:
- Matteo Arciuli <marciuli96@gmail.com>
- Luca Chiaraviglio <luca.chiaraviglio@uniroma2.it>

