import streamlit as st
from pathlib import Path

# Page configuration
st.set_page_config(
    page_title="FractoMap",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "About": "# FractoMap\nBioactivity-guided microfractionation analysis tool.\n\nDeveloped by Thapanee Pruksatrakul\n\nFunctional Metabolomics Lab, UC Riverside"
    }
)

# Custom CSS for FBMN-STATS style
st.markdown("""
<style>
    /* Main header styling */
    .main-header {
        font-size: 2.5rem;
        font-weight: 700;
        color: #1E3A5F;
        margin-bottom: 0.5rem;
    }
    
    /* Subheader styling */
    .sub-header {
        font-size: 1.1rem;
        color: #666;
        margin-bottom: 2rem;
    }
    
    /* Card styling */
    .info-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1.5rem;
        border-radius: 10px;
        color: white;
        margin-bottom: 1rem;
    }
    
    /* Feature box */
    .feature-box {
        background: #f8f9fa;
        padding: 1.2rem;
        border-radius: 8px;
        border-left: 4px solid #667eea;
        margin-bottom: 1rem;
    }
    
    /* Sidebar styling */
    [data-testid="stSidebar"] {
        background-color: #f8f9fa;
    }
    
    /* Hide Streamlit branding */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    
    /* Button styling */
    .stButton > button {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 0.5rem 2rem;
        font-weight: 500;
    }
    
    .stButton > button:hover {
        background: linear-gradient(135deg, #764ba2 0%, #667eea 100%);
    }
</style>
""", unsafe_allow_html=True)

# Sidebar logo and info
with st.sidebar:
    st.markdown("# 🧪 FractoMap")
    st.markdown("**Bioactivity-Chromatogram Overlay**")
    st.markdown("---")
    
    st.markdown("""
    ### 📚 Navigation
    Use the pages in the sidebar to:
    1. **Upload Data** - Load plate & mzML files
    2. **Plate View** - Visualize 96-well plate
    3. **Results** - View inhibition data
    4. **Overlay** - Chromatogram overlay
    """)
    
    st.markdown("---")
    st.markdown("""
    ### 🔗 Links
    - [GitHub Repository](https://github.com/LOOKiemycin/FractoMap)
    - [Documentation](https://github.com/LOOKiemycin/FractoMap#readme)
    """)
    
    st.markdown("---")
    st.markdown("""
    <small>
    Developed by<br>
    <b>Thapanee Pruksatrakul</b><br>
    Functional Metabolomics Lab<br>
    UC Riverside
    </small>
    """, unsafe_allow_html=True)

# Main page content
st.markdown('<p class="main-header">🧪 FractoMap</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Bioactivity-Guided Microfractionation Analysis Tool</p>', unsafe_allow_html=True)

# Introduction
st.markdown("""
<div class="info-card">
    <h3>Welcome to FractoMap!</h3>
    <p>Overlay bioactivity data (ABTS/DPPH antioxidant assays) with LC-MS chromatograms 
    for compound identification in natural product research.</p>
</div>
""", unsafe_allow_html=True)

# Features section
st.markdown("## ✨ Features")

col1, col2 = st.columns(2)

with col1:
    st.markdown("""
    <div class="feature-box">
        <h4>📋 Data Import</h4>
        <ul>
            <li>96-well plate Excel files</li>
            <li>mzML chromatogram files</li>
            <li>Automatic serpentine mapping</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="feature-box">
        <h4>🧮 Calculations</h4>
        <ul>
            <li>% Inhibition from absorbance</li>
            <li>Automatic control averaging</li>
            <li>Activity classification</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

with col2:
    st.markdown("""
    <div class="feature-box">
        <h4>📈 Visualization</h4>
        <ul>
            <li>96-well plate heatmap</li>
            <li>TIC/BPC chromatograms</li>
            <li>Bioactivity overlay plots</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="feature-box">
        <h4>📥 Export</h4>
        <ul>
            <li>Results as CSV</li>
            <li>Interactive plots</li>
            <li>Publication-ready figures</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

# Quick start
st.markdown("## 🚀 Quick Start")

st.markdown("""
1. **Upload Data** → Go to the Upload page and load your plate data (Excel) and chromatogram (mzML)
2. **Adjust Parameters** → Set collection start time, interval, and fraction offset
3. **View Results** → Check the Results page for inhibition calculations
4. **Generate Overlay** → Create bioactivity-chromatogram overlay plots
""")

# Workflow diagram
st.markdown("## 📊 Workflow")

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.markdown("""
    <div style="text-align: center; padding: 1rem; background: #e3f2fd; border-radius: 8px;">
        <h1 style="margin: 0;">📤</h1>
        <p><b>1. Upload</b></p>
        <small>Plate + mzML</small>
    </div>
    """, unsafe_allow_html=True)

with col2:
    st.markdown("""
    <div style="text-align: center; padding: 1rem; background: #e8f5e9; border-radius: 8px;">
        <h1 style="margin: 0;">⚙️</h1>
        <p><b>2. Configure</b></p>
        <small>Parameters</small>
    </div>
    """, unsafe_allow_html=True)

with col3:
    st.markdown("""
    <div style="text-align: center; padding: 1rem; background: #fff3e0; border-radius: 8px;">
        <h1 style="margin: 0;">🧮</h1>
        <p><b>3. Calculate</b></p>
        <small>% Inhibition</small>
    </div>
    """, unsafe_allow_html=True)

with col4:
    st.markdown("""
    <div style="text-align: center; padding: 1rem; background: #fce4ec; border-radius: 8px;">
        <h1 style="margin: 0;">📈</h1>
        <p><b>4. Overlay</b></p>
        <small>Visualize</small>
    </div>
    """, unsafe_allow_html=True)

# Citation
st.markdown("---")
st.markdown("## 📖 Citation")
st.code("""
Pruksatrakul, T. (2024). FractoMap: A tool for bioactivity-guided 
microfractionation analysis. GitHub. 
https://github.com/LOOKiemycin/FractoMap
""", language=None)

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; font-size: 0.9rem;">
    <p>FractoMap v1.0 | Developed by Thapanee Pruksatrakul</p>
    <p>Functional Metabolomics Lab, UC Riverside | JGSEE, KMUTT</p>
</div>
""", unsafe_allow_html=True)
