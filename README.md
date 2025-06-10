# SUnWNS_TV: Papanicolaou stain Unmixing Using Nonnegativity, Weighted Nucleus Sparsity, and Total Variation Regularization

This repository contains the MATLAB implementation of the **SUnWNS_TV** algorithm, a stain unmixing method tailored for Papanicolaou-stained RGB images. The algorithm incorporates **nonnegativity**, **weighted nucleus norm sparsity**, and **total variation (TV)** to robustly decompose stain abundances of 4 dyes from RGB images.

---

## 🔍 Overview

**Key Features:**
- Converts RGB images to optical density (OD) space using the Beer–Lambert law.
- Estimates stain abundance matrix under:
  - Nonnegativity constraint
  - Weighted nucleus sparsity promoting low-rank structure in each stain
  - Total variation for spatial smoothness
- Implements an **ADMM**-based optimization strategy for efficient convergence.

---

## 📁 Directory Structure

```
SUnWNS_TV/
├── test.m                   # Main script for running the algorithm
├── SUnWNS_TV.m              # Core unmixing algorithm
├── Dependency/
│   ├── RGB_absor.mat        # stain absorption coefficient matrix
│   └── Calibration.mat      # calibration parameters
├── data/
│   └── Case05_E-L_04.mat    # RGB intensities and ground truth of stain abundance
└── README.md
```

---

### Run the Code

In test.m

You can modify the parameters.

| Parameter | Description | Default |
|----------|-------------|---------|
| `Lam1`     | for weighted nucleus sparsity | `2e-6` |
| `LamTV`    | for TV regularization         | `1e-3` |
| `gamma`    | augmented Lagrangian weight   | `0.01` |
| `maxIter`  | Max iterations for ADMM       | `2000` |

---

## 📈 Output

- Estimated stain abundance
- Quantitative plot

---

## 📄 Reference

xxxxx

---

## 📬 Contact

For questions, please contact:

**Nanxin Gong**  
Department of Information and Communications Engineering  
Institute of Science Tokyo  
📧 gong.n.aa@m.titech.ac.jp

---