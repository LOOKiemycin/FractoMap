#!/usr/bin/env python3
"""
Example: Bioactivity-Chromatogram Overlay
=========================================

This example demonstrates how to overlay ABTS antioxidant activity data
with LC-MS/MS chromatogram from a Quercetin standard injection.

Data:
- Quercetin 2 mg/mL standard (Negative mode)
- ABTS assay after 6 min incubation
- Microfractionation: 86 fractions, 7 sec each, starting at 1 min
"""

import numpy as np
import matplotlib.pyplot as plt
from bioactivity_overlay import BioactivityOverlay

# =============================================================================
# STEP 1: Define your plate reader data
# =============================================================================

# Raw absorbance data from plate reader (8x12 matrix)
# This is the data you showed in the image (Quercetin ABTS 6 min)
plate_data = np.array([
    # Row A (fractions 1-12, left to right)
    [0.7652, 0.7666, 0.7676, 0.7605, 0.7601, 0.7554, 0.75,   0.77,   0.767,  0.7654, 0.7666, 0.7689],
    # Row B (fractions 24-13, right to left in serpentine)
    [0.7382, 0.7308, 0.7453, 0.7619, 0.76,   0.756,  0.7564, 0.7697, 0.7737, 0.7662, 0.7646, 0.7685],
    # Row C (fractions 25-36)
    [0.7386, 0.7352, 0.7396, 0.7293, 0.7234, 0.7265, 0.7394, 0.737,  0.7311, 0.7145, 0.6896, 0.6455],
    # Row D (fractions 48-37)
    [0.6478, 0.5982, 0.5473, 0.4867, 0.3537, 0.1848, 0.1932, 0.2508, 0.3457, 0.5255, 0.6003, 0.627],
    # Row E (fractions 49-60)
    [0.636,  0.6728, 0.6914, 0.7255, 0.7354, 0.7364, 0.7498, 0.7582, 0.7569, 0.7574, 0.751,  0.7654],
    # Row F (fractions 72-61)
    [0.7581, 0.7628, 0.7678, 0.768,  0.7687, 0.7658, 0.7668, 0.7664, 0.7622, 0.7645, 0.7509, 0.7691],
    # Row G (fractions 73-84)
    [0.7597, 0.7604, 0.7665, 0.7706, 0.7721, 0.766,  0.7649, 0.7711, 0.7644, 0.7765, 0.7659, 0.7723],
    # Row H (fractions 96-85, controls in 87-96)
    [0.7776, 0.7789, 0.7782, 0.777,  0.7746, 0.7762, 0.7709, 0.7765, 0.7761, 0.7727, 0.7641, 0.7456],
])

# =============================================================================
# STEP 2: Calculate inhibition from plate data
# =============================================================================

# Control wells are fractions 87-96 (serpentine: H10, H9, H8, H7, H6, H5, H4, H3, H2, H1)
# In the serpentine pattern, Row H goes from 96 (column 1) to 85 (column 12)
# So controls (87-96) are in H1-H10

# Control values (H1-H10, reading right-to-left for serpentine row H)
control_wells = plate_data[7, :10]  # H1 to H10
control_avg = np.mean(control_wells)
print(f"Control average (T₀): {control_avg:.4f}")

# Calculate inhibition for fractions 1-86 using serpentine pattern
def serpentine_to_inhibition(plate, control_avg, num_fractions=86):
    """Convert serpentine plate data to inhibition values."""
    inhibition = []
    
    for row in range(8):
        if row % 2 == 0:  # Even rows: left to right (A, C, E, G)
            for col in range(12):
                frac_num = row * 12 + col + 1
                if frac_num <= num_fractions:
                    abs_val = plate[row, col]
                    inh = (1 - abs_val / control_avg) * 100
                    inhibition.append(max(0, inh))  # Clamp negative to 0
        else:  # Odd rows: right to left (B, D, F, H)
            for col in range(11, -1, -1):
                frac_num = row * 12 + (11 - col) + 1
                if frac_num <= num_fractions:
                    abs_val = plate[row, col]
                    inh = (1 - abs_val / control_avg) * 100
                    inhibition.append(max(0, inh))
    
    return np.array(inhibition)

inhibition_data = serpentine_to_inhibition(plate_data, control_avg)

# Print summary
print(f"\nInhibition Summary:")
print(f"  Max: {inhibition_data.max():.1f}% at fraction {np.argmax(inhibition_data)+1}")
print(f"  Active fractions (>50%): {np.sum(inhibition_data > 50)}")

# =============================================================================
# STEP 3: Create the overlay plot
# =============================================================================

# Initialize overlay object
overlay = BioactivityOverlay(
    collection_start=1.0,       # Start collecting at 1 min
    collection_interval=7/60,   # 7 seconds per fraction
    num_fractions=86,
    fraction_offset=-1          # Correct for ~7 sec dead volume delay
)

# Load mzML file (use your actual path)
mzml_path = "/mnt/user-data/uploads/Quercetin_2mg-ml__Neg_FC_02.mzML"
overlay.load_mzml(mzml_path, chromatogram_type="TIC")

# Load inhibition data
overlay.inhibition = inhibition_data

# Create the overlay plot
fig = overlay.plot_overlay(
    output_path="/home/claude/quercetin_bioactivity_overlay.png",
    title="Quercetin Standard - ABTS Activity vs LC-MS/MS Chromatogram",
    rt_range=(1, 11),  # Collection window
    highlight_active=True,
    activity_threshold=50.0,
    show_fraction_lines=False
)

plt.show()

# =============================================================================
# STEP 4: Export results to CSV
# =============================================================================

overlay.export_results("/home/claude/quercetin_fraction_results.csv")

# =============================================================================
# STEP 5: Find compounds in active fractions
# =============================================================================

print("\n" + "="*60)
print("Compounds in Active Fractions")
print("="*60)

active_compounds = overlay.find_active_compounds(
    mzml_path,
    activity_threshold=50.0
)

if len(active_compounds) > 0:
    # Show top compounds per fraction
    for frac in active_compounds['Fraction'].unique():
        frac_data = active_compounds[active_compounds['Fraction'] == frac].head(5)
        print(f"\nFraction {frac} ({frac_data['% Inhibition'].iloc[0]}% inhibition):")
        print(f"  RT: {frac_data['RT (min)'].iloc[0]:.2f} min")
        print("  Top m/z values:")
        for _, row in frac_data.iterrows():
            print(f"    m/z {row['m/z']:.4f} (intensity: {row['Intensity']:.0f})")

# =============================================================================
# BONUS: Create EIC overlay for Quercetin identification
# =============================================================================

# Quercetin [M-H]- = 301.0348
quercetin_mz = 301.0348

# Create EIC overlay
fig2 = overlay.plot_eic_overlay(
    mzml_path,
    target_mz_list=[
        (301.0348, "Quercetin [M-H]⁻"),
        (300.0275, "Quercetin-H₂ [M-H]⁻"),  # Example fragment
    ],
    output_path="/home/claude/quercetin_eic_overlay.png",
    mz_tolerance=0.01,
    rt_range=(1, 11)
)

print("\n✅ Analysis complete!")
print("Output files:")
print("  - quercetin_bioactivity_overlay.png")
print("  - quercetin_eic_overlay.png") 
print("  - quercetin_fraction_results.csv")
