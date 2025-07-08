# HTHP
This directory contains the dynamic models developed and used to produce the results in Vecchi et al. [ref]. These packages require prior installation of **OpenModelica** [1] with default libraries, as well as the installation and import of the **ExternalMedia** library [2].

## Dynamic Models

The dynamic models are organised into two main packages: **Fluids** and **NewHTHP**.

### Fluids Package

The *Fluids* package contains some predefined fluid packages which allow evaluation of the fluid thermophysical properties as a function of the thermodynamic state. The following fluid packages are predefined:

- Water  
- R1233zd(E)  
- R1234ze(E)  
- R601 (n-pentane)  
- R601a (iso-pentane)

Additional fluids among those in the *CoolProp* library [3] can be defined by extending the `ExternalMedia.Media.CoolPropMedium` package in an analogous fashion.

### NewHTHP Package

The `NewHTHP` package contains the component models necessary for bottom-up modelling of the High Temperature Heat Pump (HTHP).  
The following models are included:

- **Machines**: a compressor  
- **Heat exchangers**: a co-flow evaporator; a counterflow condenser; and a counterflow internal heat exchanger (IHX)  
  - *The full mathematical formulation and validation of the counterflow condenser model are presented in Vecchi et al. [ref]*  
- **Valves**: a Joule-Thomson valve  
- **Volumes**: a steam accumulator and a mixed tank  
  - *The full mathematical formulation and validation of the steam accumulator model are presented in Vecchi et al. [ref]*  
- **Boundaries**: different types of sources and sinks  
- **Networks**: different types of Tee-junctions  
- **Utilities**: a pressure reduction valve (PRV)  
- **Connections**: fluid inlet and outlet connectors used in the above models


## Examples

The `Examples` directory includes 6 integrated models used to generate the results in Vecchi et al. [ref].

- **`Condenser`**: HTHP only; can simulate *Turn-up* or *Turn-down* cases depending on the simulation setup  
- **`Condenser_and_SA`**: HTHP + Steam Accumulator (SA) operated between 2 and 4 bar  
  - *Turn-up* and *Turn-down* are separate model setups: `TU` and `TD`  
- **`Condenser_and_8barSA`**: HTHP + SA operated between 4 and 8 bar  
  - *Turn-up* and *Turn-down* are separate model setups: `TU` and `TD`  

These integrated models are fully defined and ready to run. They use standard Modelica blocks and the custom models developed in this work.  
Further integrated models can be defined by extending or combining the components provided.

## References

[1] Fritzson P, Aronsson P, Pop A, Lundvall N, Nystrom K, Saldamli L, et al.  
*OpenModelica - A free open-source environment for system modeling, simulation, and teaching.*  
IEEE Conference on Computer Aided Control System Design, IEEE International Conference on Control Applications, IEEE International Symposium on Intelligent Control, Munich, Germany: IEEE; 2006.

[2] Casella F, Richter C.  
*ExternalMedia: A Library for Easy Re-use of External Fluid Property Code in Modelica.*  
Proceedings of the 6th International Modelica Conference, Bielefeld: 2008, p. 157–161.

[3] Bell IH, Wronski J, Quoilin S, Lemort V.  
*Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp.*  
Ind Eng Chem Res 2014;53:2498–508. https://doi.org/10.1021/ie4033999
