# FractoMap
Bioactivity-Chromatogram Overlay Tool

A Python tool for **bioactivity-guided fractionation** analysis. Overlays LC-MS/MS chromatograms with antioxidant activity data from microfractionation experiments.

## 🎯 Purpose

This tool helps identify bioactive compounds by correlating:
- **LC-MS/MS chromatogram** data (TIC, BPC, or EIC from mzML files)
- **Antioxidant activity** data (ABTS/DPPH % inhibition from microfractionation)

Perfect for **functional metabolomics** and **natural product discovery**!

## ✨ Features

- 📊 **Overlay plots**: TIC/BPC chromatogram with % inhibition bars
- 🎨 **Activity color-coding**: Strong (green) → Inactive (red)
- 🔍 **EIC extraction**: Create extracted ion chromatograms for compound identification
- 🐍 **Serpentine plate mapping**: Automatic conversion from 96-well plate reader data
- 📁 **Multiple input formats**: mzML, CSV, Excel, numpy arrays
- 📤 **Export results**: CSV/Excel with RT-fraction-activity mapping

### Python API

```python
from bioactivity_overlay import BioactivityOverlay

# Initialize
overlay = BioactivityOverlay(
    collection_start=1.0,       # Start collecting at 1 min
    collection_interval=7/60,   # 7 seconds per fraction
    num_fractions=86,
    fraction_offset=-1          # Correct for dead volume delay
)

# Load data
overlay.load_mzml("data.mzML", chromatogram_type="TIC")
overlay.load_inhibition_data("inhibition.csv")

# Create overlay plot
overlay.plot_overlay(
    output_path="overlay.png",
    title="My Bioactivity Analysis",
    rt_range=(1, 11),
    activity_threshold=50.0
)

# Export results
overlay.export_results("results.csv")
```

## ⚠️ Important: Fraction Offset Correction

Due to **dead volume** in the tubing between the splitter and fraction collector, there is typically a delay between MS detection and fraction collection. Use the `fraction_offset` parameter to correct this:

| Offset | Effect |
|--------|--------|
| 0 | No correction |
| -1 | Shift bars 1 fraction earlier (default, ~7 sec delay) |
| -2 | Shift bars 2 fractions earlier (~14 sec delay) |

**How to determine the correct offset:**
1. Inject a standard compound (e.g., Quercetin)
2. Compare the TIC peak RT with the bioactivity peak position
3. Adjust offset until they align

## 📋 Plate Layout

The tool uses **serpentine pattern** for 96-well plate mapping:

```
Row A: 1  → 2  → 3  → 4  → 5  → 6  → 7  → 8  → 9  → 10 → 11 → 12
Row B: 24 ← 23 ← 22 ← 21 ← 20 ← 19 ← 18 ← 17 ← 16 ← 15 ← 14 ← 13
Row C: 25 → 26 → 27 → 28 → 29 → 30 → 31 → 32 → 33 → 34 → 35 → 36
...
```

- **Wells 1-86**: Fractions (collected every 7 seconds from 1-11 min)
- **Wells 87-96**: Controls (solvent + reagent, 0% inhibition reference)

## 📊 Activity Classification

| % Inhibition | Classification | Color |
|--------------|----------------|-------|
| > 75% | **Strong** | 🟢 Green |
| 50-75% | **Moderate** | 🟡 Yellow |
| 25-50% | **Weak** | 🟠 Orange |
| < 25% | **Inactive** | 🔴 Red |

## 🔬 Workflow

```
┌─────────────────┐     ┌─────────────────┐
│  LC-MS/MS with  │     │   96-well plate │
│ microfractionation│     │  bioassay       │
│     (mzML)      │     │  (ABTS/DPPH)    │
└────────┬────────┘     └────────┬────────┘
         │                       │
         ▼                       ▼
┌─────────────────────────────────────────┐
│      Bioactivity-Chromatogram Overlay   │
│              (this tool)                │
└────────────────────┬────────────────────┘
                     │
         ┌───────────┴───────────┐
         ▼                       ▼
┌─────────────────┐     ┌─────────────────┐
│  Overlay Plot   │     │ Active Compound │
│                 │     │ Identification  │
└─────────────────┘     └─────────────────┘
```

## 📁 Example Files

- `bioactivity_overlay.py` - Main module
- `example_quercetin.py` - Example with Quercetin standard
- `requirements.txt` - Python dependencies

## 🧪 Example Output

Using Quercetin standard (2 mg/mL, ABTS assay 6 min):

```
Compounds in Active Fractions:

Fraction 43 (76.2% inhibition):
  RT: 6.01 min
  Top m/z values:
    m/z 301.0352 (Quercetin [M-H]⁻)
    m/z 603.0786 (Quercetin dimer [2M-H]⁻)
```

## 📚 References

1. Chaves N, et al. (2020) Quantification of the Antioxidant Activity of Plant Extracts. *Antioxidants* 9(1):76.
2. Re R, et al. (1999) Antioxidant Activity Applying an Improved ABTS Radical Cation Decolorization Assay. *Free Radical Biology and Medicine* 26(9):1231–37.

## Author
Thapanee Pruksatrakul
(Visiting scholar) Functional Metabolomics Laboratory  
University of California, Riverside

## 📄 License

MIT License - feel free to use and modify!
Made with ❤️ for metabolomics research
