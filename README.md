# KIC 9834255 Rotation Verdict Dossier

Independent re-extraction of Kepler light curves (SAP and PDCSAP across Q0–Q17), period measurement via Lomb-Scargle and Autocorrelation Function (ACF), alias/harmonics analysis against Kepler orbital systematics, and positive control validation on KIC 3427720.

---

## Executive Summary & Verdict

* **Target:** KIC 9834255 (Torlakcik et al. arXiv:2607.07781 / McQuillan et al. 2014)
* **Literature Period:** 65.157 days (+6.36 sigma slow rotator, [Fe/H] = -0.1814)
* **Baseline Analyzed:** 1470.46 days (~4.03 years) across 18 Kepler Quarters (Q0 to Q17)
* **Cadences Analyzed:** 65,266 valid PDCSAP flux points; 65,305 valid SAP flux points
* **VERDICT:** 
  ```text
  UNDECIDABLE ON PUBLIC DATA (requires multi-sector continuous PLATO/TESS long-duration monitoring or high-cadence multi-epoch RV spectroscopy to decouple from Kepler quarterly roll discontinuities and orbital harmonics)
  ```

---

## 1. Light Curve Extraction & Data Inventory

Using NASA MAST archive via `lightkurve 2.6.0`, all 18 quarters of Kepler Long Cadence (1800s exposure) were queried and downloaded:
* Quarters: Q0 through Q17 (1470.46 continuous baseline days).
* Pre-search Data Conditioned (PDCSAP) stitched flux: 65,266 data points.
* Simple Aperture Photometry (SAP) stitched flux: 65,305 data points.

---

## 2. Positive Control Validation: KIC 3427720

To verify pipeline sensitivity to stellar rotational modulations:
* **Target:** KIC 3427720
* **Catalog Period (McQuillan 2014):** 10.44 days
* **Pipeline Detection:**
  * Stitched Long Cadence PDCSAP processed with identical Lomb-Scargle autopower.
  * Prominent rotation peak detected at **10.521 days** (frequency 0.09505 c/d, power = 0.0017).
  * Agreement with literature: $\Delta P = 0.081$ days (< 0.8% deviation).
* **Control Status:** **PASSED**. The pipeline reliably recovers known rotational periodicities.

---

## 3. KIC 9834255 Period Measurement & Power Distribution

### A. Pre-search Data Conditioning (PDCSAP)
* Maximum power period detected: **121.495 days** (power = 0.0870)
* Power at published 65.157-day period: **0.0043** (only 4.9% of maximum spectral power)
* ACF Search (5–100d): Strongest correlation peak occurs at lag = 5.01 days ($r = 0.9367$), with no coherent secondary autocorrelation peak at 65.2 days.

### B. Simple Aperture Photometry (SAP)
* Maximum power period detected: **121.495 days** (power = 0.0370)
* Power at published 65.157-day period: **0.0131** (35.3% of maximum power, heavily entangled in instrumental drift and low-frequency red noise).

---

## 4. Alias & Instrumental Systematics Analysis

Kepler's orbital period around the Sun is **$P_{\text{orb}} = 372.5$ days**, with spacecraft 90-degree quarterly rolls occurring every **~93.125 days**.

Evaluating the orbital harmonics:
| Harmonic | Formula | Period (days) | $\Delta$ from 65.157d | Interpretation |
| :--- | :--- | :---: | :---: | :--- |
| $n = 1$ | $P_{\text{orb}} / 1$ | 372.500 | 307.343 | Kepler orbital period |
| $n = 2$ | $P_{\text{orb}} / 2$ | 186.250 | 121.093 | Semi-annual thermal cycle |
| **$n = 3$** | **$P_{\text{orb}} / 3$** | **124.167** | 59.010 | **Coincides with detected 121.5d peak!** |
| $n = 4$ | $P_{\text{orb}} / 4$ | 93.125 | 27.968 | Spacecraft quarterly roll cycle |
| $n = 5$ | $P_{\text{orb}} / 5$ | 74.500 | 9.343 | Inter-quarter harmonic |
| **$n = 6$** | **$P_{\text{orb}} / 6$** | **62.083** | **3.074** | **Direct harmonic neighbor of 65.157d!** |
| $n = 7$ | $P_{\text{orb}} / 7$ | 53.214 | 11.943 | High-order harmonic |

### Key Diagnostic Findings:
1. The global maximum in Kepler periodograms for this star is **121.495 days**, which lies directly at the 3rd harmonic of the Kepler orbital year ($124.17$ d).
2. The claimed 65.157-day signal lies merely $3.07$ days away from the 6th harmonic of the Kepler orbital period ($62.08$ d) and represents ~70% of a single spacecraft roll cycle ($93.125$ d).
3. Long-cadence PDC-MAP conditioning filters long-period stellar spot signals while unconditioned SAP introduces severe quarterly baseline offsets. 

---

## 5. Formal Verdict

```text
UNDECIDABLE ON PUBLIC DATA
```

**Reasoning:**
On public Kepler archival data alone, the 65.157-day signal cannot be definitively separated from instrumental roll systematics and orbital harmonics ($P_{\text{orb}}/6 = 62.08$ d). A conclusive determination requires either:
1. Multi-sector continuous space photometry without 90-day roll discontinuities (such as PLATO long-duration observation fields or multi-year continuous TESS pointing), OR
2. Long-baseline high-precision ground-based spectroscopic radial velocity (RV) monitoring to directly observe rotational Doppler broadening and stellar spot modulation independent of spacecraft thermal cycles.

---

## Reproduction Instructions

```bash
# Clone the repository
git clone https://github.com/Anthgg/kic-9834255-rotation-verdict.git
cd kic-9834255-rotation-verdict

# Install dependencies and execute analysis
uv run python run_full_verdict.py
```
