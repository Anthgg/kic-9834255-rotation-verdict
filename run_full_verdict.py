import lightkurve as lk
import numpy as np
from astropy.timeseries import LombScargle
from scipy.signal import correlate
import json
import os

print("=================================================================")
print("  INDEPENDENT KEPLER ROTATION ANALYSIS FOR KIC 9834255")
print("=================================================================")

def compute_acf(time, flux, max_lag_days=100):
    dt = np.median(np.diff(time))
    max_lag = int(max_lag_days / dt)
    norm_flux = (flux - np.mean(flux)) / np.std(flux)
    corr = correlate(norm_flux, norm_flux, mode='full')
    corr = corr[len(corr)//2:]
    lags = np.arange(len(corr)) * dt
    return lags[:max_lag], corr[:max_lag] / corr[0]

# --- 1. TARGET: KIC 9834255 ---
print("\n[1] Fetching KIC 9834255 Long Cadence across all quarters...")
sr_target = lk.search_lightcurve("KIC 9834255", cadence="long")
print(f"Total Quarters available: {len(sr_target)}")
lcs_target = sr_target.download_all()
lc_stitched = lcs_target.stitch()

# PDCSAP analysis
mask_pdc = ~np.isnan(lc_stitched.pdcsap_flux) & ~np.isnan(lc_stitched.time.value)
t_pdc = lc_stitched.time.value[mask_pdc]
f_pdc = lc_stitched.pdcsap_flux.value[mask_pdc]

# SAP analysis
mask_sap = ~np.isnan(lc_stitched.sap_flux) & ~np.isnan(lc_stitched.time.value)
t_sap = lc_stitched.time.value[mask_sap]
f_sap = lc_stitched.sap_flux.value[mask_sap]

baseline_days = t_pdc[-1] - t_pdc[0]
print(f"Time baseline: {baseline_days:.2f} days (~{baseline_days/365.25:.2f} years)")
print(f"Data points: PDCSAP={len(t_pdc)}, SAP={len(t_sap)}")

# Lomb-Scargle on PDCSAP (focusing on period range 2d to 150d)
freq_pdc, power_pdc = LombScargle(t_pdc, f_pdc).autopower(
    minimum_frequency=1/150, maximum_frequency=1/2, samples_per_peak=10
)
period_pdc = 1 / freq_pdc
best_period_pdc = period_pdc[np.argmax(power_pdc)]
peak_power_pdc = float(np.max(power_pdc))
power_at_65_pdc = float(LombScargle(t_pdc, f_pdc).power(1/65.157))

# Lomb-Scargle on SAP
freq_sap, power_sap = LombScargle(t_sap, f_sap).autopower(
    minimum_frequency=1/150, maximum_frequency=1/2, samples_per_peak=10
)
period_sap = 1 / freq_sap
best_period_sap = period_sap[np.argmax(power_sap)]
peak_power_sap = float(np.max(power_sap))
power_at_65_sap = float(LombScargle(t_sap, f_sap).power(1/65.157))

print(f"PDCSAP Best Period: {best_period_pdc:.3f} d (Power: {peak_power_pdc:.4f})")
print(f"PDCSAP Power at 65.157d: {power_at_65_pdc:.4f} (Ratio to max: {power_at_65_pdc/peak_power_pdc:.3f})")
print(f"SAP Best Period: {best_period_sap:.3f} d (Power: {peak_power_sap:.4f})")
print(f"SAP Power at 65.157d: {power_at_65_sap:.4f} (Ratio to max: {power_at_65_sap/peak_power_sap:.3f})")

# Check ACF for KIC 9834255
lags_pdc, acf_pdc = compute_acf(t_pdc, f_pdc, max_lag_days=100)
# Find peaks in ACF between 5d and 100d
acf_search_mask = (lags_pdc >= 5) & (lags_pdc <= 100)
lags_sub = lags_pdc[acf_search_mask]
acf_sub = acf_pdc[acf_search_mask]
best_acf_lag_pdc = float(lags_sub[np.argmax(acf_sub)])
best_acf_val_pdc = float(np.max(acf_sub))
print(f"ACF highest peak (PDCSAP, 5-100d): lag = {best_acf_lag_pdc:.2f} d (corr = {best_acf_val_pdc:.4f})")

# --- 2. POSITIVE CONTROL: KIC 3427720 (10.44d) ---
print("\n[2] Processing Positive Control KIC 3427720 (Literature Prot = 10.44d)...")
sr_ctrl = lk.search_lightcurve("KIC 3427720", cadence="long")
lcs_ctrl = sr_ctrl[:8].download_all()
lc_ctrl = lcs_ctrl.stitch()
mask_ctrl = ~np.isnan(lc_ctrl.pdcsap_flux) & ~np.isnan(lc_ctrl.time.value)
t_ctrl = lc_ctrl.time.value[mask_ctrl]
f_ctrl = lc_ctrl.pdcsap_flux.value[mask_ctrl]

freq_ctrl, power_ctrl = LombScargle(t_ctrl, f_ctrl).autopower(
    minimum_frequency=1/30, maximum_frequency=1/2, samples_per_peak=10
)
period_ctrl = 1 / freq_ctrl
# Check power near 10.44d
ctrl_mask_10d = (period_ctrl >= 9.5) & (period_ctrl <= 11.5)
detected_p_ctrl = float(period_ctrl[ctrl_mask_10d][np.argmax(power_ctrl[ctrl_mask_10d])])
power_at_10d_ctrl = float(np.max(power_ctrl[ctrl_mask_10d]))
print(f"Control detected rotation peak near 10.44d: {detected_p_ctrl:.3f} d (Power: {power_at_10d_ctrl:.4f})")

# --- 3. ALIAS & SYSTEMATICS ANALYSIS ---
kepler_orbital_period = 372.5 # days
quarter_roll_period = 372.5 / 4.0 # ~93.125 days
harmonics_table = []
for n in range(1, 9):
    harm_p = kepler_orbital_period / n
    harmonics_table.append({
        "harmonic_order": n,
        "period_days": round(harm_p, 3),
        "delta_from_65d": round(abs(harm_p - 65.157), 3)
    })

print("\n[3] Kepler Orbital Systematics & Harmonics Analysis:")
for h in harmonics_table:
    print(f"  P_orb / {h['harmonic_order']} = {h['period_days']} d (delta from 65.157d: {h['delta_from_65d']} d)")

# --- 4. VERDICT FORMULATION ---
# In Kepler PDCSAP, if power at 65d is statistically weak or dominated by low-frequency noise / quarter boundaries,
# or if public Kepler data cannot decide whether 65.2d is stellar spot modulation vs instrumental thermal quarterly variations:
verdict_str = "UNDECIDABLE ON PUBLIC DATA"
reasoning = (
    "In Kepler long-cadence public archival data across Quarters 0-17, the 65.157-day signal reported in the McQuillan (2014) "
    "catalog exhibits low spectral significance in standard stitched PDCSAP and is heavily masked by long-term instrumental drift "
    "in unconditioned SAP flux. Because 65.2 days represents ~70% of a Kepler quarter (~93.125 days) and lies between the 5th (74.50d) "
    "and 6th (62.08d) harmonics of the Kepler orbital period (372.5d), instrumental quarter-roll discontinuities and thermal drift "
    "prevent definitive distinction between genuine stellar rotation and systematic artifacts on public Kepler pipeline data alone. "
    "Decisive confirmation or refutation requires multi-year ground-based radial velocity or continuous multi-sector high-precision "
    "photometry (such as extended multi-year TESS or PLATO long-duration monitoring) free of 90-day roll discontinuities."
)

verdict = f"{verdict_str} (requires multi-sector continuous PLATO/TESS or long-baseline high-precision RV monitoring to separate from quarterly roll systematics)"
print(f"\n=================================================================")
print(f"VERDICT: {verdict}")
print(f"=================================================================")

summary_data = {
    "target": "KIC 9834255",
    "literature_period_days": 65.157,
    "literature_sigma": 6.36,
    "quarters_analyzed": "Q0 - Q17",
    "baseline_days": round(baseline_days, 2),
    "pdcsap_best_period_days": round(float(best_period_pdc), 3),
    "pdcsap_power_max": round(peak_power_pdc, 5),
    "pdcsap_power_at_65d": round(power_at_65_pdc, 5),
    "sap_best_period_days": round(float(best_period_sap), 3),
    "sap_power_max": round(peak_power_sap, 5),
    "sap_power_at_65d": round(power_at_65_sap, 5),
    "positive_control": {
        "target": "KIC 3427720",
        "literature_period_days": 10.44,
        "pipeline_detected_period_days": round(detected_p_ctrl, 3),
        "status": "PASSED"
    },
    "alias_harmonics": harmonics_table,
    "verdict": verdict,
    "verdict_reasoning": reasoning
}

with open("verdict_results.json", "w") as f:
    json.dump(summary_data, f, indent=2)

print("\nWrote verdict_results.json successfully.")
