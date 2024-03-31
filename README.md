## Response to von Schmalansee et al. 

**Authors**: Guillermo Garcia-Costoya<sup>1</sup>, Claire E. Williamns<sup>1</sup>, Trevor M. Faske<sup>1</sup>, Jacob D. Moorman<sup>2</sup>, Michael L. Logan<sup>1</sup>

**Affiliations**: <sup>1</sup>University of Nevada, Reno, NV, USA. <sup>2</sup>University of California, Los Angeles, CA, USA. 

This repository contains the data and code necessary to generate Figure 1 of the manuscript titled "**Response to von Schmalensee et al.**" by Garcia-Costoya et al. 2023. This manuscript was submitted to *Ecology Letters* as a reply to the technical note written by von Schmalensee et al. titled "Methodological artefacts [mispelling in the original title] cause counter-intuitive evolutionary conclusions in a simulation study". This note was in reference to Garcia-Costoya et al. 's 2023 *Ecology Letters* article titled "Evolutionary constraints medaite extinction risk under climate change". 

The structure of this repository is as follows: 

- `plot_fig1.R` : `R` script to replicate Figure 1.
- `re_run_sims_&_plot_fig2.R`: `R` script to re-run simulations and plot Figure 2. To plot figure 2 the `sum_sim_results.RData` data set from the [supplementary materials of our original paper](https://datadryad.org/stash/dataset/doi:10.5061/dryad.2fqz612t3) is needed,
- `data` : folder containing all data used to replicate Figure 1. This folder contains the following files:
    - `apletophallus_sprint_speed.csv` : Sprint print speed data of Panamanian slender anoles ( *Anolis apletophallus* ) used by Logan et al. 2021 *J. Exp. Biol* and Neel et al. 2021 *Biotropica*. 
    - `cage_temperatures.txt` : Cage temperature data from von Schmalansee et al. 2021 *Ecology Letters*. 
    - `eggs.txt` : Egg developmental rate data from von Schmalansee et al. 2021 *Ecology Letters*. 
    - `larvae.txt` : Larvae developmental rate data from von Schmalansee et al. 2021 *Ecology Letters*. 
    - `monthly_weather_station_data_bahamas.csv` : Weather station data from Great Exhuma, Bahamas used by Neel et al. 2021 *Biotropica*. 
    - `montly_weather_station_data_panama.csv` : Weather station data from Barro Colorado Island, Panama used by Neel et al. 2021 *Biotropica*. 
    - `sagrei_sprint_speed.csv` : Sprint speed data of Bahamian brown anoles ( *Anolis sagrei* ) used by Logan et al. 2021 *J. Exp. Biol* and Neel et al. 2021 *Biotropica*. 
    - `weather_station_temperatures.txt` : Weather station data from Floda, Västra Götaland County, Sweden from von Schmalansee et al. 2021 *Ecology Letters*.
