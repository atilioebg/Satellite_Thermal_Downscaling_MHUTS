# Data Management & Satellite Earth Observation Sources

This directory contains instructions and utilities for sourcing the multi-sensor satellite imagery required by the **MHUTS** (Modified High-resolution Urban Thermal Sharpener) pipeline.

---

## 🛰️ Required Satellite Products

The MHUTS downscaling model couples two satellite systems:

| Sensor / Satellite | Spectral Product | Native Resolution | Target Resolution | Source Platform |
| :--- | :--- | :---: | :---: | :--- |
| **Landsat-8 TIRS** | Band 10 Top-of-Atmosphere (TOA) Radiance | $100\text{ m}$ (resampled to $30\text{ m}$) | $30\text{ m}$ input grid | [USGS EarthExplorer](https://earthexplorer.usgs.gov/) |
| **Sentinel-2 MSI** | Level-2A Bottom-of-Atmosphere Reflectance (B2–B12) | $10\text{ m}$ / $20\text{ m}$ | $10\text{ m}$ master grid | [Copernicus Open Access Hub](https://dataspace.copernicus.eu/) |

---

## 📐 Data Preprocessing Workflow

1. **Sentinel-2 Bands Alignment:**
   - Resample $20\text{ m}$ bands (B5, B6, B7, B11, B12) to the $10\text{ m}$ reference grid (B2, B3, B4, B8) using bilinear or bicubic interpolation.
   - Scale reflectance values to $[0, 1]$ (divide raw ESA L2A integer values by $10{,}000$).

2. **Landsat-8 Thermal Calibration:**
   - Convert raw Digital Numbers (DN) of Band 10 to TOA Spectral Radiance ($L_\lambda$ in $\text{W}/(\text{m}^2 \cdot \text{sr} \cdot \mu\text{m})$) using the MTL metadata gains and offsets:
     $$L_\lambda = M_L \cdot \text{DN} + A_L$$
   - Estimate total atmospheric water vapor content ($\omega$) from ENVI FLAASH, MODIS atmospheric profiles, or atmospheric radiosondes ($\omega = 2.8861\text{ g/cm}^2$ for the publication scene).

3. **Spatial Co-Registration:**
   - Reproject both datasets to the local UTM coordinate system (e.g., WGS 84 / UTM Zone 22S for São Leopoldo, RS, Brazil).
   - Co-register pixels to guarantee sub-pixel spatial congruence.

---

## 🧪 Synthetic Benchmark Demo

For quick testing and validation without downloading gigabytes of satellite tiles, run:
```matlab
% In MATLAB:
cd matlab
run_demo
```
This automatically simulates a calibrated multi-sensor Earth Observation scene with identical spatial and statistical properties to the published paper.
