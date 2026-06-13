# FractoMap WebApp

A Streamlit web application for bioactivity-guided microfractionation analysis.

## 📁 Structure

```
FractoMap-WebApp/
├── app.py                      # Main entry point
├── requirements.txt            # Python dependencies
├── .streamlit/
│   └── config.toml            # Streamlit configuration
├── pages/
│   ├── 1_📤_Upload_Data.py    # Data upload page
│   ├── 2_🧫_Plate_View.py     # Plate visualization
│   ├── 3_📊_Results.py        # Results analysis
│   └── 4_📈_Overlay.py        # Chromatogram overlay
└── src/
    ├── __init__.py
    └── common.py              # Shared functions
```

## 🚀 Local Installation

```bash
# Clone repository
git clone https://github.com/LOOKiemycin/FractoMap.git
cd FractoMap/FractoMap-WebApp

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

## ☁️ Deploy to Streamlit Cloud

1. Push this folder to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your GitHub repository
4. Set main file path: `FractoMap-WebApp/app.py`
5. Click Deploy

## 📋 Features

- **Upload Data**: Load 96-well plate Excel files and mzML chromatograms
- **Plate View**: Interactive heatmap visualization
- **Results**: Statistical analysis and activity classification
- **Overlay**: Bioactivity-chromatogram overlay plots

## 📖 Usage

1. Navigate to **Upload Data** page
2. Upload your plate data (Excel) and chromatogram (mzML)
3. Adjust parameters (collection start, interval, offset)
4. Click **Calculate Inhibition**
5. View results in **Results** and **Overlay** pages

## 👩‍🔬 Author

Thapanee Pruksatrakul  
Functional Metabolomics Lab, UC Riverside  
JGSEE, King Mongkut's University of Technology Thonburi
