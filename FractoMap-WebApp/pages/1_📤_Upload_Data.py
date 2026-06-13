"""
FractoMap - Upload Data Page
"""

import streamlit as st
import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent.parent / "src"))
from common import (
    page_setup, init_session_state, parse_plate_excel, 
    parse_mzml, generate_demo_data, calculate_inhibition, ROWS
)

# Initialize
init_session_state()

# Page setup
st.set_page_config(page_title="FractoMap - Upload", page_icon="📤", layout="wide")

st.markdown("""
<style>
    .upload-box {
        background: #f8f9fa;
        padding: 2rem;
        border-radius: 10px;
        border: 2px dashed #dee2e6;
        text-align: center;
        margin-bottom: 1rem;
    }
    .success-msg {
        background: #d4edda;
        color: #155724;
        padding: 1rem;
        border-radius: 8px;
        margin-top: 1rem;
    }
    .info-box {
        background: #e3f2fd;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #2196f3;
        margin-bottom: 1rem;
    }
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("# 📤 Upload Data")
st.markdown("Upload your plate data and chromatogram files to begin analysis.")

st.markdown("---")

# Two columns for upload
col1, col2 = st.columns(2)

# Column 1: Plate Data
with col1:
    st.markdown("## 📋 Plate Data")
    
    st.markdown("""
    <div class="info-box">
        <b>Required format:</b> Excel file with 8×12 absorbance values (96-well plate)
    </div>
    """, unsafe_allow_html=True)
    
    plate_file = st.file_uploader(
        "Upload plate data",
        type=['xlsx', 'xls', 'csv'],
        key="plate_upload",
        help="96-well plate absorbance data"
    )
    
    if plate_file:
        plate_data = parse_plate_excel(plate_file)
        if plate_data:
            st.session_state.plate_data = plate_data
            st.success(f"✅ Successfully loaded 8×12 plate data")
            
            # Show preview
            with st.expander("Preview plate data"):
                import pandas as pd
                preview_df = pd.DataFrame(
                    plate_data,
                    index=ROWS,
                    columns=range(1, 13)
                )
                st.dataframe(preview_df.style.format("{:.4f}"), use_container_width=True)
        else:
            st.error("❌ Could not find 8×12 numeric matrix in file")
    
    # Status
    if st.session_state.plate_data:
        st.markdown("""
        <div class="success-msg">
            ✅ Plate data loaded
        </div>
        """, unsafe_allow_html=True)

# Column 2: Chromatogram
with col2:
    st.markdown("## 📈 Chromatogram Data")
    
    st.markdown("""
    <div class="info-box">
        <b>Supported formats:</b> mzML (TIC/BPC) or CSV (RT, Intensity)
    </div>
    """, unsafe_allow_html=True)
    
    chrom_file = st.file_uploader(
        "Upload chromatogram",
        type=['mzML', 'mzml', 'csv'],
        key="chrom_upload",
        help="mzML or CSV chromatogram file"
    )
    
    if chrom_file:
        with st.spinner("Parsing chromatogram..."):
            if chrom_file.name.lower().endswith('.mzml'):
                content = chrom_file.read().decode('utf-8')
                st.session_state.chrom_data = parse_mzml(content)
            else:
                import pandas as pd
                df = pd.read_csv(chrom_file)
                st.session_state.chrom_data['tic'] = list(zip(df.iloc[:, 0], df.iloc[:, 1]))
        
        tic_len = len(st.session_state.chrom_data.get('tic', []))
        bpc_len = len(st.session_state.chrom_data.get('bpc', []))
        
        if tic_len > 0 or bpc_len > 0:
            st.success(f"✅ TIC: {tic_len} points | BPC: {bpc_len} points")
        else:
            st.error("❌ No chromatogram data found")
    
    # Demo data button
    st.markdown("---")
    if st.button("🎲 Generate Demo Data", use_container_width=True):
        st.session_state.chrom_data = generate_demo_data()
        st.success("✅ Demo chromatogram generated")
    
    # Status
    if st.session_state.chrom_data.get('tic'):
        st.markdown("""
        <div class="success-msg">
            ✅ Chromatogram loaded
        </div>
        """, unsafe_allow_html=True)

st.markdown("---")

# Parameters section
st.markdown("## ⚙️ Parameters")

col1, col2, col3, col4 = st.columns(4)

with col1:
    collection_start = st.number_input(
        "Collection start (min)",
        value=st.session_state.params['collection_start'],
        min_value=0.0,
        max_value=30.0,
        step=0.1,
        help="When fraction collection begins"
    )
    st.session_state.params['collection_start'] = collection_start

with col2:
    collection_interval = st.number_input(
        "Interval (sec)",
        value=st.session_state.params['collection_interval'],
        min_value=1,
        max_value=60,
        step=1,
        help="Time per fraction"
    )
    st.session_state.params['collection_interval'] = collection_interval

with col3:
    fraction_offset = st.number_input(
        "Fraction offset",
        value=st.session_state.params['fraction_offset'],
        min_value=-10,
        max_value=10,
        step=1,
        help="Dead volume correction (-1 recommended)"
    )
    st.session_state.params['fraction_offset'] = fraction_offset

with col4:
    assay_type = st.selectbox(
        "Assay type",
        options=['ABTS (734 nm)', 'DPPH (515 nm)'],
        index=0 if st.session_state.params['assay_type'] == 'ABTS' else 1
    )
    st.session_state.params['assay_type'] = assay_type.split()[0]

# Calculate button
st.markdown("---")

if st.button("🧮 Calculate Inhibition", type="primary", use_container_width=True):
    if st.session_state.plate_data:
        with st.spinner("Calculating..."):
            df, control_avg = calculate_inhibition(
                st.session_state.plate_data,
                st.session_state.params
            )
            st.session_state.inhibition_df = df
            st.session_state.control_avg = control_avg
        
        st.success("✅ Calculation complete! Go to Results page to view.")
        
        # Quick stats
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Control Average", f"{control_avg:.4f}")
        col2.metric("Max Inhibition", f"{df['% Inhibition'].max():.1f}%")
        col3.metric("Active Fractions", len(df[df['% Inhibition'] > 50]))
        col4.metric("Best Fraction", f"F{df.loc[df['% Inhibition'].idxmax(), 'Fraction']}")
    else:
        st.error("❌ Please upload plate data first")

# Help section
with st.expander("ℹ️ Excel Template Help"):
    st.markdown("""
    ### Expected Format
    
    The tool reads **8 rows × 12 columns** of absorbance values:
    
    ```
          1       2       3      ...     12
    A   0.234   0.245   0.198   ...   0.312
    B   0.456   0.423   0.398   ...   0.289
    ...
    H   0.776   0.781   0.773   ...   0.778  ← Controls
    ```
    
    ### Serpentine Mapping
    
    - **Fractions 1-86**: Sample wells
    - **Fractions 87-96**: Control wells (used for 0% inhibition reference)
    
    ### Fraction Offset
    
    Due to dead volume in tubing (~7 sec delay), use `offset = -1` to align peaks correctly.
    """)
