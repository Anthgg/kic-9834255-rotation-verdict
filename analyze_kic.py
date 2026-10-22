import lightkurve as lk
import numpy as np
from astropy.timeseries import LombScargle
import json
import os

print("--- 1. PROCESSING POSITIVE CONTROL: KIC 3427720 (Expected ~10.44d) ---")
sr_ctrl = lk.search_lightcurve("KIC 3427720", mission="Kepler", author="Kepler")
print(f"Positive control quarters found: {len(sr_ctrl)}")
# Download a subset of quarters to verify pipeline sensitivity (e.g. Q3, Q4)
lc_ctrl_collection = sr_ctrl[1:5].download_all()
lc_ctrl = lc_ctrl_collection.stitch().remove_nans().remove_outliers()

# Run Lomb-Scargle on control
freq_ctrl, power_ctrl = LombScargle(lc_ctrl.time.value, lc_ctrl.pdcsap_flux.value).autopower(
    minimum_frequency=1/50, maximum_frequency=1/1
)
best_period_ctrl = 1 / freq_ctrl[np.argmax(power_ctrl)]
print(f"Positive Control Best Period detected: {best_period_ctrl:.2f} days (Catalog: 10.44 days)")

print("\n--- 2. PROCESSING TARGET: KIC 9834255 ---")
sr_target = lk.search_lightcurve("KIC 9834255", mission="Kepler", author="Kepler")
print(f"Target quarters found: {len(sr_target)}")
lc_target_collection = sr_target.download_all()

# Extract both SAP and PDCSAP
lc_sap = lc_target_collection.stitch(flux_column="sap_flux").remove_nans().remove_outliers()
lc_pdcsap = lc_target_collection.stitch(flux_column="pdcsap_flux").remove_nans().remove_outliers()

print(f"Total SAP points: {len(lc_sap)}, Total PDCSAP points: {len(lc_pdcsap)}")
time_span = lc_pdcsap.time.value[-1] - lc_pdcsap.time.value[0]
print(f"Time baseline: {time_span:.1f} days")

# Periodogram on PDCSAP
freq_pdc, power_pdc = LombScargle(lc_pdcsap.time.value, lc_pdcsap.flux.value).autopower(
    minimum_frequency=1/200, maximum_frequency=1/5
)
best_p_pdc = 1 / freq_pdc[np.argmax(power_pdc)]
max_power_pdc = np.max(power_pdc)
print(f"PDCSAP Best Period: {best_p_pdc:.3f} days (Peak Power: {max_power_pdc:.4f})")

# Periodogram on SAP
freq_sap, power_sap = LombScargle(lc_sap.time.value, lc_sap.flux.value).autopower(
    minimum_frequency=1/200, maximum_frequency=1/5
)
best_p_sap = 1 / freq_sap[np.argmax(power_sap)]
max_power_sap = np.max(power_sap)
print(f"SAP Best Period: {best_p_sap:.3f} days (Peak Power: {max_power_sap:.4f})")

# Alias / Harmonics analysis
kepler_year = 372.5
quarter_period = 372.5 / 4.0 # ~93.125 days
harmonics = [kepler_year / i for i in range(1, 10)]
print(f"Kepler year: {kepler_year}d, Quarter roll: {quarter_period:.2f}d")
print("Harmonics of Kepler year:", [round(h, 2) for h in harmonics])

# Check power at 65.157d specifically
power_at_65_pdc = LombScargle(lc_pdcsap.time.value, lc_pdcsap.flux.value).power(1/65.157)
power_at_65_sap = LombScargle(lc_sap.time.value, lc_sap.flux.value).power(1/65.157)
print(f"Power at 65.157d (PDCSAP): {float(power_at_65_pdc):.5f}")
print(f"Power at 65.157d (SAP): {float(power_at_65_sap):.5f}")

# Save results summary to JSON
results = {
    "target": "KIC 9834255",
    "control": "KIC 3427720",
    "control_detected_period": round(float(best_period_ctrl), 3),
    "control_expected_period": 10.44,
    "target_baseline_days": round(float(time_span), 1),
    "target_best_period_pdcsap": round(float(best_p_pdc), 3),
    "target_best_period_sap": round(float(best_p_sap), 3),
    "target_power_at_65_pdc": round(float(power_at_65_pdc), 5),
    "target_power_at_65_sap": round(float(power_at_65_sap), 5),
    "kepler_quarter_roll_days": round(quarter_period, 3),
    "kepler_year_days": kepler_year
}
with open("analysis_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("Analysis complete. Saved to analysis_results.json")
