# Top Quark Pair Production Analysis Using CMS Open Data

## Overview

This project focuses on the analysis of top quark pair production (\(t\bar{t}\)) using proton-proton collision data from the **CMS experiment at CERN**.

The analysis uses CMS Run 2016 DoubleMuon data and Monte Carlo simulations to investigate top quark pair production in the dilepton decay channel, particularly events containing two muons.

The main objectives are to identify candidate \(t\bar{t}\) events, suppress background contributions, study relevant kinematic distributions, and compare experimental observations with simulated predictions.

The analysis was conducted as an undergraduate particle physics course project using **ROOT and C++**.

## Physics Background

Top quarks are predominantly produced in pairs through strong interactions in proton-proton collisions at the LHC.

Each top quark decays predominantly into a W boson and a bottom quark:

\[
t\rightarrow W^+b,\qquad
\bar{t}\rightarrow W^-\bar{b}
\]

In the dimuon channel considered in this analysis, both W bosons decay leptonically:

\[
t\bar{t}\rightarrow
b\bar{b}\mu^+\mu^-\nu_\mu\bar{\nu}_\mu
\]

The resulting final state contains:

- Two oppositely charged muons
- Jets originating from bottom quarks
- Missing transverse momentum associated with neutrinos

These characteristics provide the basis for identifying top quark pair production events.

## Dataset and Monte Carlo Samples

### Experimental Data

- **Experiment:** CMS
- **Collision system:** Proton-proton
- **Center-of-mass energy:** 13 TeV
- **Data-taking period:** Run 2016
- **Primary dataset:** DoubleMuon, Run2016H

### Monte Carlo Simulation

Simulated samples are used to model the expected signal and background contributions.

The investigated processes include:

- \(t\bar{t}\) production
- Drell–Yan production
- Single top production
- Diboson production (WW, WZ, ZZ)

Monte Carlo event yields are normalized using the relevant production cross sections, integrated luminosity, and number of generated events.

## Analysis Methodology

### 1. Event Selection

A sequence of selection criteria is applied to identify candidate top quark pair events.

The analysis considers:

- Muon multiplicity and charge requirements
- Muon transverse momentum and pseudorapidity
- Jet multiplicity and kinematic requirements
- Dilepton invariant mass
- Missing transverse energy (MET)
- Event-level kinematic variables

A Z boson mass veto is applied to suppress Drell–Yan background events.

### 2. Cutflow Analysis

A cutflow analysis is performed to monitor the number of events remaining after each selection requirement.

This allows the effects of individual selection criteria to be investigated and provides a framework for evaluating signal and background reduction.

### 3. Kinematic Distributions

Selected events are analyzed through distributions of relevant observables, including:

- Dimuon invariant mass
- Jet multiplicity
- Missing transverse energy (MET)
- Scalar sum of transverse momenta (HT)
- Other reconstructed event-level quantities

These distributions are used to study the characteristics of candidate \(t\bar{t}\) events.

### 4. Data and Monte Carlo Comparison

Experimental data are compared with simulated signal and background samples.

The analysis includes:

- Monte Carlo event normalization
- Background contributions
- Comparison of reconstructed kinematic distributions
- Data/MC ratios
- Investigation of agreement between measured and simulated event yields

This comparison is used to evaluate the modeling of the selected events and identify differences between experimental observations and simulation predictions.

### 5. Top Quark Reconstruction Studies

Additional ROOT macros investigate reconstructed top quark mass distributions and related fitting procedures.

These studies complement the event selection and kinematic analyses.

## Software and Tools

- **ROOT:** Data processing, histograms, fitting, and visualization
- **C++:** Event selection and analysis algorithms
- **CMS Open Data:** Experimental collision data and simulated samples

## Repository Structure

The repository contains ROOT/C++ analysis macros and selected analysis outputs.

Key analysis files include:

- `ttbar_dilepton_analysis.C` — Main dilepton analysis
- `cutflow.C` — Event selection and cutflow analysis
- `reco_top_mass.C` — Top quark mass reconstruction studies
- `fit_top_mass.C` — Mass distribution fitting
- `cross_section.C` — Cross-section-related analysis
- `cutflow.csv` — Event selection summary
- `plots/` — Analysis plots and distributions

Large ROOT datasets are not included in the repository.

## Results and Discussion

The analysis investigates the effects of event selection and background suppression on the identification of top quark pair candidates.

Normalized Monte Carlo samples are compared with CMS collision data through kinematic distributions and event yields.

The project provides practical experience in collider data analysis, including event reconstruction, background modeling, statistical comparisons, and the interpretation of experimental particle physics data.

## Data Availability

The analysis is based on publicly available CMS collision data and Monte Carlo samples provided through the CERN Open Data Portal.

**Data source:** https://opendata.cern.ch/

## Research Context

This work was completed as an undergraduate course project in particle physics at **Middle East Technical University (METU)**.
