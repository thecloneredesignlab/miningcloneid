# Gemcitabine PK/PD Kinetics

## Overview
This directory contains experimental data and mathematical modeling parameters for Gemcitabine pharmacokinetics and pharmacodynamics (PK/PD). The modeling effort combines empirical kinetic rates with clinical trial scaling parameters to predict drug behavior, specifically using the clinical data to generalize the experimental kinetics to arbitrary administered doses.

---

## 1. In Vitro Experiments (Drug Build-up and Decay)
**Purpose:** To establish the baseline kinetic rates—specifically, how fast the drug builds up in the system and subsequently decays over time.

Two experiments were conducted to capture these kinetics:
* **24-hour experiment:** Administered dose was 0.1 mM
* **48-hour experiment:** Administered dose was 1.0 mM

---

## 2. Clinical Trial Data (Mapping Intake to Peak Concentration)
**Purpose:** To map the administered drug intake to the peak concentration (`C_peak`). This allows us to generalize the build-up and decay kinetics observed *in vitro* to arbitrary doses (such as the 30, 60, and 120 mg/kg doses used in mouse models).

Within a moderate dose range, clinical Phase I PK/PD studies indicate that the relationship between drug intake and peak concentration is a linear, dose-proportional function:

`C_peak = a + b * drug_intake`

Using approximate calculations derived from these clinical datasets and established pharmacological trends, the specific equations for different biological compartments are:

* **Plasma Concentration:** 
  `C_peak = 0.1943 * drug_intake`

* **Tumor Interstitial Fluid (TIF) Concentration:** 
  `C_peak = 0.0971 * drug_intake`
  *(Context: TIF drug concentrations are typically 2-5% of serum levels).*

* **Intracellular Concentration:** 
  `C_peak = 0.00971 * drug_intake`
  *(Context: Intracellular concentration is lower than TIF, utilizing a 10:1 TIF-to-intracellular ratio here).*

Source: https://www.tandfonline.com/doi/abs/10.1179/joc.2007.19.2.212
