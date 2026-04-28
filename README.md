# Automotive Radar Antenna Simulation Using CST and MATLAB

## Project Introduction

This repository contains the MATLAB scripts and CST-exported antenna gain data used in the third-year individual project **Automotive Antenna Design and Performance Analysis Using CST and MATLAB**.

This project investigates the performance of square, circular, and hexagonal patch antennas for automotive radar applications at **25 GHz**. CST Studio Suite is used to model and simulate individual patch antennas, while MATLAB is used to build array models, import CST-derived gain data, generate coverage diagrams and gain maps, and analyse representative automotive radar scenarios.

The MATLAB part first uses isotropic antenna elements to establish baseline radar array models. These models are used to analyse the spatial power density distribution under single-element, five-element, twenty-element, and phased array conditions. The CST-exported realistic antenna gain data are then introduced into the MATLAB array model to compare the coverage and gain distribution characteristics of different patch antenna geometries under array conditions.

## Repository Structure

```text
baseline_radar_models/
    Baseline MATLAB array models used in Chapter 3.
    These include the single isotropic antenna model, five-element array model,
    twenty-element array model, and phased array examples.

cst_gain_data/
    CST-exported antenna gain data files.
    These include selected cut-plane gain data for the square, circular,
    and hexagonal patch antennas.

antenna_array_analysis/
    MATLAB array analysis scripts used in Chapter 4.
    These scripts import CST-derived gain data and generate coverage diagrams
    and gain maps for the three patch antennas.

scenario_simulation/
    Representative automotive radar scenario simulation scripts used in Chapter 5.
    These include front/rear arrays, side arrays, target vehicles, bicycles,
    shadow effects, and blind region analysis.

## Required Software

The code was developed using MATLAB.

Required software:

- MATLAB

## Input Data

The CST-derived antenna gain files are used as input data for the MATLAB array model.

The main input files are:

```text
SQU_Gain_Theta25.xlsx
CIR_Gain_Phi90.xlsx
HEX_Gain_Phi0.xlsx
```

Each CST-exported gain file should follow the basic format below:

```text
Column 1: Theta [deg]
Column 2: Phi [deg]
Column 3: Gain or directivity [dBi]
```

In MATLAB, these data are read and filtered to obtain the corresponding selected cut-plane data. The gain values are then converted from dBi to linear scale, and an interpolation method is used to construct an element gain function. This allows the element gain to be estimated at different observation angles.

## How to Run the Code

1. Open MATLAB.
2. Open the repository folder.
3. First run the scripts in `baseline_radar_models/` to generate the Chapter 3 baseline isotropic antenna and array model results.
4. Then run the scripts in `antenna_array_analysis/` to generate coverage diagrams and gain maps using CST-derived antenna gain data.
5. Finally, run the scripts in `scenario_simulation/` to generate received power distribution plots for the representative automotive radar scenarios.

Recommended running order:

```text
1. baseline_radar_models/
2. antenna_array_analysis/
3. scenario_simulation/
```

The scripts generate the following types of result figures:

- single isotropic antenna power density distribution
- five-element isotropic array power density distribution
- twenty-element isotropic array power density distribution
- phased array power density distribution
- coverage diagrams
- gain maps
- scenario-based received power distribution plots

## Technical Details

### Baseline radar array models

The baseline models in Chapter 3 are built using isotropic antenna elements. The single-element model is used to show the basic free-space power density attenuation with distance. The five-element and twenty-element models further consider the coherent superposition of multiple elements in a two-dimensional observation plane.

In the multi-element models, MATLAB calculates the propagation distance, single-element power density, and propagation phase from each element to each observation point. The model does not directly add the power from each element. Instead, coherent field summation is performed using both amplitude and phase, and the total power density is then obtained from the resultant field.

Compared with the five-element model, the twenty-element model shows a more distinct spatial structure and a more concentrated high-power region in the forward direction. Therefore, it is used as the main baseline array model for the later introduction of CST-derived antenna radiation patterns and scenario analysis.

### Phased array examples

Based on the twenty-element isotropic array model, the code further introduces an adjacent-element feed phase difference to observe the effect of phase control on the direction of the spatial power density distribution.

Fixed phase difference examples, such as 60° and 120°, are used in the project. The results show that, when the array geometry remains unchanged, changing the feed phase difference between adjacent elements can significantly change the direction of the high-power-density region. This part provides a basis for the later beam scanning analysis.

### CST-derived gain data integration

In Chapter 4, the square, circular, and hexagonal patch antennas are first modelled and simulated in CST. Representative cut angles are then selected from the CST far-field results, and the corresponding gain data are exported.

In MATLAB, the CST-exported gain data are converted into linear gain and used to construct an interpolated element gain function. In this way, when the array model calculates the contribution of each element to an observation point, it can query the corresponding antenna gain according to the observation direction instead of continuing to use the isotropic element assumption.

The representative cut-plane data used for the three antennas are:

```text
Square patch antenna: theta = 25°
Circular patch antenna: phi = 90°
Hexagonal patch antenna: phi = 0°
```

### Coverage diagram and gain map

The coverage diagram represents the maximum received power distribution that can be achieved at each spatial sample point under beam-scanning conditions. In the MATLAB implementation, multiple steering angles are scanned within a predefined steering angle range, and the maximum received power is retained at each sample point.

The gain map records the equivalent array gain corresponding to the maximum received power at each sample point. Therefore, the coverage diagram mainly reflects the maximum detection capability distribution, while the gain map mainly reflects the array gain distribution characteristics. These two results are related, but they represent different physical quantities.

### Scenario simulations

The scenario simulations in Chapter 5 are used to analyse the received power distribution of different array configurations in representative automotive scenarios. The ego vehicle and target vehicles are simplified as rectangular models, while bicycle targets are simplified as line targets.

Scenario 1 uses front and rear arrays with circular patch antennas. It is used to analyse the response of front and rear targets and the shadow effect caused by target occlusion.

Scenario 2a still uses only front and rear arrays. It is used to show that a lateral target is difficult to detect effectively when side arrays are not included.

Scenario 2b adds side arrays to the left and right sides of the ego vehicle. It is used to verify the improvement provided by lateral arrays for side target detection. In this scenario, the front/rear arrays use circular patch antennas, while the side arrays use square patch antennas.

Scenario 3 uses front, rear, left, and right array configurations. It is used to show that even with four-direction arrays, the current array configuration and scanning range still cannot completely eliminate local blind regions.

## Assumptions and Limitations

The current MATLAB model includes the following simplifying assumptions:

- The observation space is simplified as a two-dimensional plane.
- A free-space propagation assumption is mainly used.
- The effect of the vehicle body on antenna radiation and electromagnetic propagation is simplified.
- Radar cross section is simplified for relative comparison.
- Target vehicles are represented using sampled rectangular edge points.
- Bicycle targets are represented using simplified line targets.
- The shadow effect is mainly based on geometric line-of-sight blockage judgement.
- CST-derived gain data are taken from selected cut planes rather than full 3D antenna pattern interpolation.

These assumptions make the model more suitable for comparing the relative performance of different antenna geometries and array configurations, rather than predicting the absolute performance of a real production automotive radar system.

## Future Improvements

Future work could be extended in the following directions:

- investigating more antenna geometries, such as triangular, elliptical, or other patch structures;
- further analysing different array configurations, including different numbers of elements, non-uniform element spacing, MIMO systems, and planar arrays;
- optimising the array arrangement and scanning range to further reduce blind regions around the vehicle;
- introducing full-vehicle models to further consider vehicle body structures, mounting positions, and interactions between multiple antennas;
- extending the analysis to dynamic scenarios involving moving targets and multiple targets;
- introducing machine learning or optimisation algorithms to automatically optimise antenna dimensions, array structures, or element excitations.

## Author

Wenhao Zhu
BEng
The University of Manchester  
2025–2026
