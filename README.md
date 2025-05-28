# SEIR-model

This repository presents an extended SEIR (Susceptible-Exposed-Infectious-Removed) epidemic model, developed to analyze and predict the evolution of infectious diseases, with a particular focus on the COVID-19 pandemic in Russia and Finland.

## Overview

Compartmental transmission models are powerful tools to study the dynamics of epidemics. This project uses a seasonally forced SEIR model, which incorporates chaotic behavior and noise sensitivity, making it more biologically realistic. The model's key objectives are to study the evolution of COVID-19 using extended SEIR dynamics and estimate key parameters that govern epidemic progression.

## Application Context

The model is applied to COVID-19 data:
- **Countries Analyzed**: Russia and Finland
- **Data Period**: Until June 12, 2021
- **Source**: WHO statistics and other public datasets

The model uses curve-fitting techniques to estimate how the disease evolves and assess healthcare demands.

## Repository Structure

```
SEIR-model/
│
├── report/        # Project report including detailed methodology and results
├── codes/         # Source code for model implementation and simulations
├── data/          # Raw and processed data used for analysis
├── figures/       # Visualizations such as epidemic curves and model diagrams
└── README.md      # Project description and documentation
```

## Usage

To run the model and reproduce the results:

1. Clone the repository:
   ```bash
   git clone https://github.com/nehabinish/SEIR-model.git
   cd SEIR-model
   ```
2. Navigate to the `codes/` directory and run the simulation scripts.
3. Review the `report/` for an in-depth explanation and findings.

> This project is a part of the M2 CompuPhys masters program at [UBFC](https://www.ubfc.fr/en/).

 
