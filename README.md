# Instrumented Indentation Data Analysis

This repository contains a Python code developed for the processing and analysis of **instrumented indentation test data**, aiming at the extraction of mechanical properties from force–depth curves.

The methodology implemented in this project is based on the following reference:

> **Instrumented indentation for the evaluation of mechanical properties of materials**  
> https://www.sciencedirect.com/science/article/pii/S0167844222004931

This work was developed as part of an internship / undergraduate research project, following technical guidelines provided by the academic advisor.

---

## 🎯 Project Objectives

- Import and preprocess experimental indentation data;
- Automatically identify loading and unloading cycles;
- Calibrate penetration depth to correct instrumental displacement;
- Compute unloading stiffness and contact geometry parameters;
- Extract representative stress and strain values;
- Fit a power-law constitutive material model using least squares;
- Estimate mechanical properties such as:
  - Elastic modulus
  - Constitutive model parameters
  - Yield stress and yield strain
- Generate graphical visualizations of the results.

---

## 🧪 Implemented Methodology

The computational workflow follows these main steps:

1. **Experimental data loading**
   - Import force (`load`) and depth (`depth`) columns from a `.txt` file.

2. **Pre-processing**
   - System compliance correction;
   - Unit conversion;
   - Automatic identification of loading/unloading cycles.

3. **Depth calibration**
   - Linear regression applied to 70–95% of the first loading cycle;
   - Horizontal shifting of the curve so that the extrapolated line passes through the origin.

4. **Unloading stiffness calculation**
   - Linear regression on the upper unloading segment (≥ 75% of maximum force).

5. **Representative stress and strain calculation**
   - Based on geometric and mechanical relations from the literature.

6. **Constitutive model fitting**
   - Power-law model:
     \[
     \sigma = K \cdot \varepsilon^n
     \]
   - Parameter estimation via nonlinear least squares (`scipy.optimize.least_squares`).

7. **Additional property estimation**
   - Elastic modulus;
   - Yield parameters;
   - Final stress–strain curve.

---

## 📊 Generated Outputs

The code automatically produces:

- **Force × Depth plot**
  - Raw experimental data
  - Calibrated curve
  - Maximum and minimum points of each cycle

- **Stress × Strain plot**
  - Fitted experimental points
  - Initial elastic response

---

## 🧰 Technologies Used

- **Python 3**
- **NumPy** – numerical computation
- **SciPy** – optimization and curve fitting
- **Matplotlib** – data visualization

---

