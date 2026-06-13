"""
FractoMap - Plate View Page
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent.parent / "src"))
from common import init_session_state, ROWS, get_serpentine_mapping

# Initialize
init_session_state()
SERPENTINE = get_serpentine_mapping()

# Page setup
st.set_page_config(page_title="FractoMap - Plate", page_icon="🧫", layout="wide")

st.markdown("""
<style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .plate-legend {
        display: flex;
        gap: 20px;
        margin-top: 1rem;
        justify-content: center;
    }
    .legend-item {
        display: flex;
        align-items: center;
        gap: 8px;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("# 🧫 96-Well Plate View")
st.markdown("Visualize your plate data with serpentine fraction mapping.")

st.markdown("---")

if st.session_state.plate_data is None:
    st.warning("⚠️ No plate data loaded. Please go to **Upload Data** page first.")
    st.stop()

plate_data = st.session_state.plate_data
inhibition_df = st.session_state.inhibition_df

# Create plate heatmap
st.markdown("## Plate Heatmap")

# Build matrix for heatmap
plate_matrix = np.array(plate_data)

# Create hover text
hover_text = []
for r in range(8):
    row_text = []
    for c in range(12):
        # Find fraction number
        frac = None
        for f in range(1, 97):
            if SERPENTINE[f]['row'] == r and SERPENTINE[f]['col'] == c:
                frac = f
                break
        
        val = plate_data[r][c]
        well = f"{ROWS[r]}{c+1}"
        
        if frac and frac <= 86 and inhibition_df is not None:
            inh = inhibition_df[inhibition_df['Fraction'] == frac]['% Inhibition'].values
            if len(inh) > 0:
                text = f"Well: {well}<br>Fraction: F{frac}<br>Abs: {val:.4f}<br>Inhibition: {inh[0]:.1f}%"
            else:
                text = f"Well: {well}<br>Fraction: F{frac}<br>Abs: {val:.4f}"
        else:
            text = f"Well: {well}<br>Fraction: F{frac}<br>Abs: {val:.4f}<br>(Control)"
        
        row_text.append(text)
    hover_text.append(row_text)

# Create heatmap
fig = go.Figure(data=go.Heatmap(
    z=plate_matrix,
    x=list(range(1, 13)),
    y=ROWS,
    colorscale='RdYlGn_r',
    hovertext=hover_text,
    hovertemplate='%{hovertext}<extra></extra>',
    colorbar=dict(title="Absorbance")
))

fig.update_layout(
    title="96-Well Plate Absorbance",
    xaxis=dict(title="Column", tickmode='linear', dtick=1),
    yaxis=dict(title="Row", autorange='reversed'),
    height=400,
    width=800
)

st.plotly_chart(fig, use_container_width=True)

# Legend
st.markdown("""
<div class="plate-legend">
    <div class="legend-item">🟢 Low absorbance (high inhibition)</div>
    <div class="legend-item">🟡 Medium absorbance</div>
    <div class="legend-item">🔴 High absorbance (low inhibition)</div>
</div>
""", unsafe_allow_html=True)

st.markdown("---")

# Inhibition heatmap (if calculated)
if inhibition_df is not None:
    st.markdown("## Inhibition Heatmap")
    
    # Build inhibition matrix
    inh_matrix = np.zeros((8, 12))
    for _, row in inhibition_df.iterrows():
        r, c = row['Row'], row['Col']
        inh_matrix[r][c] = row['% Inhibition']
    
    # Mark controls
    for f in range(87, 97):
        r, c = SERPENTINE[f]['row'], SERPENTINE[f]['col']
        inh_matrix[r][c] = np.nan  # Mark as NaN for different color
    
    # Hover text for inhibition
    hover_text_inh = []
    for r in range(8):
        row_text = []
        for c in range(12):
            frac = None
            for f in range(1, 97):
                if SERPENTINE[f]['row'] == r and SERPENTINE[f]['col'] == c:
                    frac = f
                    break
            
            well = f"{ROWS[r]}{c+1}"
            
            if frac and frac <= 86:
                inh = inhibition_df[inhibition_df['Fraction'] == frac]['% Inhibition'].values
                if len(inh) > 0:
                    text = f"Well: {well}<br>Fraction: F{frac}<br>Inhibition: {inh[0]:.1f}%"
                else:
                    text = f"Well: {well}<br>Fraction: F{frac}"
            else:
                text = f"Well: {well}<br>(Control)"
            
            row_text.append(text)
        hover_text_inh.append(row_text)
    
    fig_inh = go.Figure(data=go.Heatmap(
        z=inh_matrix,
        x=list(range(1, 13)),
        y=ROWS,
        colorscale=[
            [0, '#E24B4A'],      # 0% - Red
            [0.25, '#F5A623'],   # 25% - Orange
            [0.5, '#EF9F27'],    # 50% - Yellow
            [0.75, '#7CB342'],   # 75% - Light green
            [1.0, '#1D9E75']     # 100% - Green
        ],
        zmin=0,
        zmax=100,
        hovertext=hover_text_inh,
        hovertemplate='%{hovertext}<extra></extra>',
        colorbar=dict(title="% Inhibition")
    ))
    
    fig_inh.update_layout(
        title="% Inhibition by Well",
        xaxis=dict(title="Column", tickmode='linear', dtick=1),
        yaxis=dict(title="Row", autorange='reversed'),
        height=400,
        width=800
    )
    
    st.plotly_chart(fig_inh, use_container_width=True)

st.markdown("---")

# Serpentine pattern reference
with st.expander("📋 Serpentine Pattern Reference"):
    st.markdown("""
    ### Fraction → Well Mapping
    
    ```
    Row A: →  F1   F2   F3   F4   F5   F6   F7   F8   F9  F10  F11  F12
    Row B: ← F24  F23  F22  F21  F20  F19  F18  F17  F16  F15  F14  F13
    Row C: → F25  F26  F27  F28  F29  F30  F31  F32  F33  F34  F35  F36
    Row D: ← F48  F47  F46  F45  F44  F43  F42  F41  F40  F39  F38  F37
    Row E: → F49  F50  F51  F52  F53  F54  F55  F56  F57  F58  F59  F60
    Row F: ← F72  F71  F70  F69  F68  F67  F66  F65  F64  F63  F62  F61
    Row G: → F73  F74  F75  F76  F77  F78  F79  F80  F81  F82  F83  F84
    Row H: ← F96  F95  F94  F93  F92  F91  F90  F89  F88  F87  F86  F85
                                        ↑________________________↑
                                             Controls (87-96)
    ```
    """)

# Raw data table
with st.expander("📊 Raw Plate Data"):
    raw_df = pd.DataFrame(
        plate_data,
        index=ROWS,
        columns=range(1, 13)
    )
    st.dataframe(raw_df.style.format("{:.4f}"), use_container_width=True)
    
    # Download button
    csv = raw_df.to_csv()
    st.download_button(
        "📥 Download Plate Data",
        csv,
        "plate_data.csv",
        "text/csv"
    )
