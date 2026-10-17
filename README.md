# Satellite Thermal Downscaling via MHUTS (Landsat-8 & Sentinel-2)

[![Paper DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.rsase.2022.100702-blue.svg)](https://doi.org/10.1016/j.rsase.2022.100702)
[![Journal](https://img.shields.io/badge/Journal-Elsevier%20RSASE-orange.svg)](https://www.sciencedirect.com/journal/remote-sensing-applications-society-and-environment)
[![Language](https://img.shields.io/badge/Language-MATLAB%20(Vectorized)-green.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Official Implementation of the Research Publication:**  
> **"Use of the downscalling method in satellite images for the analysis of heat islands in landfills"**  
> *Remote Sensing Applications: Society and Environment (Elsevier)*, Vol. 26, 2022, 100702.  
> **Authors:** Atílio Grondona, Luciana P. Gomes, Luis M. Schiavo, Marcelo Caetano, Bruno J.B.L. Barbosa.  
> **Direct DOI:** [10.1016/j.rsase.2022.100702](https://doi.org/10.1016/j.rsase.2022.100702)

---

## 📌 Executive Summary

Thermal Infrared (TIR) Earth Observation is indispensable for environmental monitoring, industrial landfill management, and urban heat island mitigation. In solid waste landfills, anaerobic biological decomposition generates exothermic biogas emissions, resulting in pronounced surface thermal anomalies. 

However, satellite thermal infrared imagery is severely constrained by spatial resolution tradeoffs (Landsat-8 TIRS: native $100\text{ m}$, resampled to $30\text{ m}$). This spatial smearing obscures localized hot spots across individual landfill cells and leachate evaporation lagoons.

This codebase provides an end-to-end multi-sensor super-resolution framework coupling **Landsat-8 Thermal Infrared Radiance** with **Sentinel-2 Multi-Spectral Imagery (10 m)** via the **Modified High-resolution Urban Thermal Sharpener (MHUTS)**. The algorithm achieves a **300% spatial resolution enhancement ($30\text{ m} \to 10\text{ m}$)** while rigorously preserving radiometric physics, reducing localized temperature retrieval error to an **RMSE of $1.03\text{ K}$**.

---

## 🔬 Mathematical & Physical Formulation

```
                                  MULTI-SENSOR DATA FUSION PIPELINE
                                  
  Sentinel-2 MSI (10 m)                                               Landsat-8 TIRS (30 m)
  ┌──────────────────────────────────────────────┐                   ┌──────────────────────────────┐
  │ 9 Optical / SWIR Surface Reflectance Bands   │                   │ Band 10 TOA Radiance (L_10)  │
  └──────────────────────┬───────────────────────┘                   └──────────────┬───────────────┘
                         │                                                          │
          ┌──────────────┴──────────────┐                                           │
          ▼                             ▼                                           │
  Broadband Surface Albedo       Biophysical Predictors                             │
  (D'Urso & Calera, 2006)       (Liu & Zhang; Stathopoulou)                         │
     [ α at 10 m ]                  [ PV, ε at 10 m ]                               │
          │                             │                                           │
          │                             │ Resample ε to 30 m                        │
          │                             └───────────────────► GSCM Atmospheric LST  │
          │                                                   (Jiménez-Muñoz, 2014) │
          │                                                       [ LST at 30 m ]   │
          │                                                                 │
          │                                                                 ▼
          │                                                      Factor-3 Resampling (Kronecker)
          │                                                                 │
          └─────────────────────────────┬───────────────────────────────────┘
                                        ▼
                  Least-Squares Multivariate Cubic Regression (MHUTS)
                  fitlm([PV, ε, α], LST_coarse, 'poly333') -> 20 parameters
                                        │
                                        ▼
                  High-Resolution Sharpened Thermal Map (10 m)
                           [ LST_10m (RMSE ≤ 1.03 K) ]
```

### 1. Land Surface Temperature (LST) Retrieval — GSCM
Top-of-Atmosphere (TOA) spectral radiance is atmospherically corrected via the Generalized Single-Channel Method (Jiménez-Muñoz et al., 2014):

$$LST = \gamma \left[ \frac{\psi_1 L_\lambda + \psi_2}{\varepsilon} + \psi_3 \right] + \delta$$

where:
* $L_\lambda$: Landsat-8 Band 10 spectral radiance ($\text{W}/(\text{m}^2 \cdot \text{sr} \cdot \mu\text{m})$).
* $\psi_1, \psi_2, \psi_3$: Atmospheric functions modeled as 2nd-degree polynomials of total atmospheric water vapor ($\omega = 2.8861\text{ g/cm}^2$, derived via ENVI FLAASH).
* $\gamma, \delta$: Planck linearization parameters dependent on effective brightness temperature $T_b$:
  $$\gamma = \frac{T_b^2}{b_\gamma L_\lambda}, \quad \delta = T_b - \frac{T_b^2}{b_\gamma} \quad (b_\gamma = 1324\text{ K})$$

### 2. Multi-Spectral Biophysical Predictors (Sentinel-2 at 10 m)
* **Broadband Surface Albedo ($\alpha$):** Derived from 9 Sentinel-2 optical/SWIR bands via D'Urso & Calera (2006):
  $$\alpha = 0.1324 \rho_{B2} + 0.1269 \rho_{B3} + 0.1051 \rho_{B4} + 0.0971 \rho_{B5} + 0.0890 \rho_{B6} + 0.0818 \rho_{B7} + 0.0722 \rho_{B8} + 0.0167 \rho_{B11} + 0.0002 \rho_{B12}$$

* **Fractional Vegetation Cover ($PV$):** Scaled quadratic vegetation proportion (Stathopoulou & Cartalis, 2009):
  $$PV = \left( \frac{\text{NDVI} - \text{NDVI}_{\min}}{\text{NDVI}_{\max} - \text{NDVI}_{\min}} \right)^2$$

* **Land Surface Emissivity ($\varepsilon$):** Piecewise semi-empirical formulation (Liu & Zhang, 2011):
  $$\varepsilon = \begin{cases} 
  0.990 & \text{if } \text{NDVI} \ge 0.727 \\ 
  1.0094 + 0.047 \ln(\text{NDVI}) & \text{if } 0.157 \le \text{NDVI} < 0.727 \\ 
  0.970 & \text{if } -0.185 \le \text{NDVI} < 0.157 \\ 
  0.995 & \text{if } \text{NDVI} < -0.185 
  \end{cases}$$

### 3. MHUTS Multivariate Polynomial Sharpener
A 3rd-degree trivariate polynomial in biophysical space ($PV, \varepsilon, \alpha$) is fitted using the Least Squares Method (LSM) over the $3\times 3$ replicated coarse grid, yielding 20 estimated parameters ($P_1$ to $P_{20}$):

$$LST_{10\text{m}} = \sum_{i+j+k \le 3} P_{ijk} \cdot (PV)^i \cdot (\varepsilon)^j \cdot (\alpha)^k$$

---

## 📊 Scientific Validation Results (São Leopoldo Landfill AOI)

The algorithm was validated over the CRVR Solid Waste Landfill (São Leopoldo, Rio Grande do Sul, Brazil). The table below replicates the exact empirical results published in **Table 3** of the Elsevier RSASE paper:

| Scope / Experimental Condition | Resolution | Mean LST (K) | Std Dev (K) | RMSE (K) | Max LST (K) | Min LST (K) | Accuracy Gain |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Full Scene (Landsat-8 GSCM)** | 30 m | 298.20 | 2.52 | — | 307.87 | 279.80 | Baseline |
| **Full Scene (MHUTS Downscaled)**| 10 m | 298.20 | 2.03 | **1.51** | 334.18 | 290.60 | — |
| **Landfill AOI (Landsat-8 GSCM)** | 30 m | 296.03 | 1.62 | — | 304.44 | 293.42 | Baseline |
| **Landfill AOI (MHUTS Downscaled)**| **10 m** | **296.32** | **1.27** | **1.03** | **302.70** | **293.87** | **31.79% Error Reduction** |

### Key Physical Insights:
* **Spatial Super-Resolution:** Localized optimization across active waste disposal phases (Phases 1 to 6) lowered the RMSE from $1.51\text{ K}$ to $1.03\text{ K}$ (a **31.8% accuracy gain**).
* **Thermal Anomaly Isolation:** Resolved distinct thermal boundaries between recently deposited municipal waste (exothermic biological digestion) and inactive covered phases, which were completely blended in standard $30\text{ m}$ Landsat pixels.

---

## 📁 Repository Structure

```
Satellite_Thermal_Downscaling_MHUTS/
├── matlab/                         # MATLAB Engine
│   ├── src/                        # High-Performance Vectorized Library
│   │   ├── albedo_sentinel.m       # Broadband Albedo (D'Urso / Liang / Knap)
│   │   ├── vegetation_emissivity.m # Vectorized PV and Liu & Zhang Emissivity
│   │   ├── lst_gscm.m              # Generalized Single-Channel Method (GSCM LST)
│   │   ├── mhuts_downscale.m       # Trivariate Cubic Polynomial Sharpener (poly333)
│   │   ├── evaluate_metrics.m      # Table 3 Statistical Metric Suite
│   │   └── pipeline_downscale.m    # Integrated Multi-Sensor Execution Wrapper
│   │
│   ├── legacy/                     # Original 2022 Publication Code (Historical Provenance)
│   │   ├── MMQ3.m                  # Original Full-Scene LSM Script
│   │   ├── MMQ4.m                  # Original Landfill Sub-AOI Script (Table 3)
│   │   ├── calc_temp_updated.m     # Original GSCM LST Loop
│   │   ├── sentinel_pve.m          # Original Emissivity Loop
│   │   ├── TES.m                   # Reference Channel Method
│   │   └── UHS.m                   # Urban Thermal Inversion
│   │
│   └── run_demo.m                  # One-Click Reproducible Validation Script
│
├── data/                           # Sourcing guides and synthetic benchmark
│   └── README.md                   # Instructions for Landsat-8 and Sentinel-2 L2A scenes
│
├── docs/                           # Published Paper and Documentation
│   └── Use of the downscalling method in satellite images for the analysis.pdf
│
├── figures/                        # High-resolution benchmark figures
├── CITATION.cff                    # Academic citation metadata
├── LICENSE                         # MIT License
└── README.md                       # Main Scientific Documentation
```

---

## 🚀 Quickstart: Reproducing Table 3 in MATLAB

To reproduce the downscaling pipeline and Table 3 validation without downloading external satellite scenes, execute the self-contained demonstration:

```matlab
% In MATLAB command window:
cd matlab
run_demo
```

This script:
1. Simulates a calibrated Landsat-8 ($30\text{ m}$) and Sentinel-2 ($10\text{ m}$) scene matching the physical characteristics of the São Leopoldo landfill.
2. Runs the modernized, vectorized MHUTS downscaling engine.
3. Prints the formatted validation table directly to the console.
4. Generates a multi-panel visual comparison saved to `figures/demo_thermal_downscaling.png`.

---

## 📖 Citation

If you use this algorithm, code, or methodology in your research, please cite our paper:

```bibtex
@article{grondona2022use,
  title={Use of the downscalling method in satellite images for the analysis of heat islands in landfills},
  author={Grondona, At{\'i}lio and Gomes, Luciana P and Schiavo, Luis M and Caetano, Marcelo and Barbosa, Bruno JBL},
  journal={Remote Sensing Applications: Society and Environment},
  volume={26},
  pages={100702},
  year={2022},
  publisher={Elsevier},
  doi={10.1016/j.rsase.2022.100702}
}
```

---

## 📄 License
This project is licensed under the [MIT License](LICENSE).
