For running these script either a docker environment needs to setup or a workstation with openMDAO and ADFlow with other utilities (pygeo, PyHyp etc is need).
Each folder contains the script for a different objective.

### Project Overview
In this project we implement the Aerodynamics Shape Optimization (ASO) of a
727 200 Adv wing in cruise conditions. The wing is designed in OpenVSP with the

sketches and dimensions available in literature. In this optimization problem we in-
troduce another discipline/component in the form of performance. At first solely the

aerodynamics shape optimization is done with the objective to minimize the coefficient

of drag for a given coefficient of lift constraint and structural sizing variables consider-
ing aeroelastic deflections using the framework for multidisciplinary design optimiza-
tion (MDO) of aircraft configurations with high fidelity (MACH).Then we introduce

changes in objective of the previous problem and make it Range (with volume of fuel as
constraint) of the aircraft and mass of fuel burnt (using range constraint) respectively,
this helps arrive at a design that considers the performance.

#### Surface and Vol mesh
![image](https://github.com/vishwaaaaa/High-Fidelity-Aerodynamics-Shape-Optimization/assets/122648757/56ff2090-3976-439b-b716-91f5c5ccc56f)
Upon creating the geometry in OpenVSP the surface mesh was created in Pointewise.
Further the volume mesh is generated using `pyHyp`

#### Optimization Architecture 
![image](https://github.com/vishwaaaaa/High-Fidelity-Aerodynamics-Shape-Optimization/assets/122648757/2dd57194-470f-4a84-a9bf-6ad8cccd97f6)

#### Optimization Formulation
![image](https://github.com/vishwaaaaa/High-Fidelity-Aerodynamics-Shape-Optimization/assets/122648757/9801a579-9ba7-4f23-a9e5-6ec1fd25681b)
