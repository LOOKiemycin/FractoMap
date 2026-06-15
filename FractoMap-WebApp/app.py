"""
FractoMap - Bioactivity-Guided Microfractionation Analysis
Styled Streamlit App
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import base64
import zlib
import struct
import xml.etree.ElementTree as ET

# Page Config
st.set_page_config(page_title="FractoMap", page_icon="🧪", layout="wide")

# Custom CSS
st.markdown("""
<style>
    /* Header banner */
    .header-banner {
        background: linear-gradient(135deg, #667eea 0%, #1DB984 100%);
        padding: 1.2rem 2rem;
        border-radius: 16px;
        margin-bottom: 1.5rem;
        color: white;
    }
    .header-banner h1 {
        font-size: 2.2rem;
        font-weight: 700;
        margin: 0 0 0.5rem 0;
    }
    .header-banner p {
        font-size: 1rem;
        opacity: 0.95;
        margin: 0;
    }
    
    /* Card styling */
    .upload-card {
        background: white;
        border: 1px solid #e5e7eb;
        border-radius: 12px;
        padding: 1.5rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.08);
    }
    .card-title {
        font-size: 1.1rem;
        font-weight: 600;
        color: #1f2937;
        margin-bottom: 0.5rem;
    }
    .card-subtitle {
        font-size: 0.9rem;
        color: #6b7280;
        margin-bottom: 1rem;
    }
    
    /* Tags */
    .tag {
        display: inline-block;
        padding: 4px 10px;
        border-radius: 6px;
        font-size: 0.8rem;
        font-weight: 500;
        margin-right: 6px;
    }
    .tag-teal { background: #d1fae5; color: #065f46; }
    .tag-blue { background: #dbeafe; color: #1e40af; }
    
    /* Hide default streamlit */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .stDeployButton {display: none;}
    
    /* Tab styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        border-radius: 8px;
        padding: 10px 20px;
        font-weight: 500;
    }
    
    /* Metric cards */
    .metric-row {
        display: flex;
        gap: 1rem;
        margin: 1rem 0;
    }
    .metric-card {
        background: #f9fafb;
        border-radius: 10px;
        padding: 1rem 1.25rem;
        flex: 1;
    }
    .metric-label {
        font-size: 0.8rem;
        color: #6b7280;
        margin-bottom: 4px;
    }
    .metric-value {
        font-size: 1.5rem;
        font-weight: 700;
        color: #111827;
    }
</style>
""", unsafe_allow_html=True)

# Constants
ROWS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

def get_serpentine_mapping():
    mapping = {}
    frac = 1
    for r in range(8):
        if r % 2 == 0:
            for c in range(12):
                mapping[frac] = {'row': r, 'col': c, 'well': f"{ROWS[r]}{c+1}"}
                frac += 1
        else:
            for c in range(11, -1, -1):
                mapping[frac] = {'row': r, 'col': c, 'well': f"{ROWS[r]}{c+1}"}
                frac += 1
    return mapping

SERPENTINE = get_serpentine_mapping()

# Helper Functions
def decode_binary(base64_string, is_64bit=True, is_zlib=False):
    decoded = base64.b64decode(base64_string)
    if is_zlib:
        try: decoded = zlib.decompress(decoded)
        except: pass
    if is_64bit:
        return struct.unpack(f'<{len(decoded)//8}d', decoded)
    return struct.unpack(f'<{len(decoded)//4}f', decoded)

def parse_mzml(content):
    try:
        root = ET.fromstring(content)
        result = {'tic': [], 'bpc': []}
        for elem in root.iter():
            if 'chromatogram' in elem.tag.lower():
                chrom_id = elem.get('id', '').lower()
                time_array, intensity_array = None, None
                is_64bit, is_zlib = True, False
                for ba in elem.iter():
                    if 'binaryDataArray' in ba.tag:
                        is_time, is_int = False, False
                        for cv in ba.iter():
                            if 'cvParam' in cv.tag:
                                n = cv.get('name', '').lower()
                                if '64-bit' in n: is_64bit = True
                                elif '32-bit' in n: is_64bit = False
                                if 'zlib' in n: is_zlib = True
                                if 'time array' in n: is_time = True
                                if 'intensity array' in n: is_int = True
                        for b in ba.iter():
                            if 'binary' in b.tag.lower() and b.text:
                                d = decode_binary(b.text.strip(), is_64bit, is_zlib)
                                if is_time: time_array = d
                                elif is_int: intensity_array = d
                if time_array and intensity_array:
                    times = [t/60 if t > 100 else t for t in time_array]
                    data = list(zip(times, intensity_array))
                    if 'tic' in chrom_id or 'total' in chrom_id: result['tic'] = data
                    elif 'bpc' in chrom_id or 'base' in chrom_id: result['bpc'] = data
                    elif not result['tic']: result['tic'] = data
        return result
    except Exception as e:
        st.error(f"Error parsing mzML: {e}")
        return {'tic': [], 'bpc': []}

def parse_plate_excel(file):
    try:
        df = pd.read_csv(file, header=None) if file.name.endswith('.csv') else pd.read_excel(file, header=None)
        matrix = []
        for i in range(len(df)):
            row = df.iloc[i].values
            nums = []
            start = 1 if (len(row) > 0 and isinstance(row[0], str) and row[0].strip().upper() in ROWS) else 0
            for j in range(start, min(start+12, len(row))):
                if pd.notna(row[j]):
                    try: nums.append(float(row[j]))
                    except: pass
            if len(nums) >= 10:
                while len(nums) < 12: nums.append(0.0)
                matrix.append(nums[:12])
                if len(matrix) == 8: break
        return matrix if len(matrix) == 8 else None
    except Exception as e:
        st.error(f"Error: {e}")
        return None

def calculate_inhibition(plate_data, params):
    start, interval, offset = params['start'], params['interval'], params['offset']
    ctrl = [plate_data[SERPENTINE[f]['row']][SERPENTINE[f]['col']] for f in range(87, 97)]
    ctrl_avg = np.mean(ctrl)
    results = []
    for f in range(1, 87):
        r, c = SERPENTINE[f]['row'], SERPENTINE[f]['col']
        abs_val = plate_data[r][c]
        inh = max(0, (1 - abs_val / ctrl_avg) * 100)
        rt = start + (f - 1 + 0.5 + offset) * (interval / 60)
        act = 'Strong' if inh > 75 else 'Moderate' if inh > 50 else 'Weak' if inh > 25 else 'Inactive'
        results.append({'Fraction': f, 'Well': SERPENTINE[f]['well'], 'RT (min)': round(rt, 3),
                       'Absorbance': round(abs_val, 4), '% Inhibition': round(inh, 1), 'Activity': act})
    return pd.DataFrame(results), ctrl_avg

def generate_demo():
    tic = [(rt, 1e6 + np.random.random()*1e5 + 8e6*np.exp(-((rt-5.8)/0.15)**2) + 3e6*np.exp(-((rt-4.5)/0.2)**2)) 
           for rt in np.arange(0, 12, 0.02)]
    return {'tic': tic, 'bpc': [(t, i*0.5) for t, i in tic]}

def get_color(inh):
    if inh > 75: return '#2E8B57'
    elif inh > 50: return '#FFD700'
    elif inh > 25: return '#FF8C00'
    return '#DC143C'

# Session State
for k, v in [('plate', None), ('chrom', {'tic': [], 'bpc': []}), ('df', None), ('ctrl', None), ('gnps', None)]:
    if k not in st.session_state: st.session_state[k] = v

def parse_gnps_excel(file):
    """Parse GNPS2 annotation Excel file"""
    try:
        df = pd.read_excel(file) if file.name.endswith(('.xlsx', '.xls')) else pd.read_csv(file)
        
        # Find RT column
        rt_cols = [c for c in df.columns if any(x in c.lower() for x in ['rt', 'retention', 'time'])]
        rt_col = rt_cols[0] if rt_cols else None
        
        # Find compound name column
        name_cols = [c for c in df.columns if any(x in c.lower() for x in ['compound', 'name', 'library', 'annotation'])]
        name_col = name_cols[0] if name_cols else None
        
        # Find m/z column
        mz_cols = [c for c in df.columns if any(x in c.lower() for x in ['mz', 'mass', 'precursor'])]
        mz_col = mz_cols[0] if mz_cols else None
        
        # Find score column
        score_cols = [c for c in df.columns if any(x in c.lower() for x in ['score', 'cosine', 'mq'])]
        score_col = score_cols[0] if score_cols else None
        
        if not rt_col or not name_col:
            st.error(f"Cannot find RT or Compound columns. Found: {list(df.columns)[:10]}")
            return None
        
        results = []
        for _, row in df.iterrows():
            rt = row[rt_col]
            name = row[name_col]
            if pd.notna(rt) and pd.notna(name) and str(name).strip():
                # Convert RT to minutes if in seconds
                rt_val = float(rt)
                if rt_val > 100:  # likely in seconds
                    rt_val = rt_val / 60
                
                results.append({
                    'RT': rt_val,
                    'Compound': str(name).strip(),
                    'mz': float(row[mz_col]) if mz_col and pd.notna(row[mz_col]) else None,
                    'Score': float(row[score_col]) if score_col and pd.notna(row[score_col]) else None
                })
        
        return pd.DataFrame(results) if results else None
    except Exception as e:
        st.error(f"Error parsing GNPS file: {e}")
        return None

# Header Banner with custom icon
st.markdown("""
<div class="header-banner">
    <div style="display: flex; align-items: center; gap: 16px;">
        <div style="background: rgba(255,255,255,0.2); border-radius: 12px; padding: 10px; display: flex; align-items: center; justify-content: center;">
            <svg width="36" height="36" viewBox="0 0 100 100" fill="none" xmlns="http://www.w3.org/2000/svg">
                <!-- Flask body -->
                <path d="M35 15 L35 40 L15 85 Q12 95 22 95 L78 95 Q88 95 85 85 L65 40 L65 15" 
                      fill="rgba(255,255,255,0.9)" stroke="white" stroke-width="3"/>
                <!-- Flask neck -->
                <rect x="32" y="8" width="36" height="12" rx="3" fill="white"/>
                <!-- Liquid gradient -->
                <path d="M20 75 Q25 65 35 70 Q50 78 65 68 Q75 62 80 75 L78 90 Q78 92 76 92 L24 92 Q22 92 22 90 Z" 
                      fill="url(#liquidGrad)"/>
                <!-- Bubbles -->
                <circle cx="35" cy="78" r="4" fill="rgba(255,255,255,0.6)"/>
                <circle cx="55" cy="82" r="3" fill="rgba(255,255,255,0.5)"/>
                <circle cx="45" cy="72" r="2.5" fill="rgba(255,255,255,0.4)"/>
                <circle cx="62" cy="76" r="2" fill="rgba(255,255,255,0.5)"/>
                <!-- Peak lines (chromatogram) -->
                <path d="M25 55 L35 55 L40 35 L45 50 L55 25 L60 50 L65 45 L75 55" 
                      stroke="#FFD700" stroke-width="3" fill="none" stroke-linecap="round" stroke-linejoin="round"/>
                <defs>
                    <linearGradient id="liquidGrad" x1="0%" y1="0%" x2="100%" y2="100%">
                        <stop offset="0%" style="stop-color:#34D399"/>
                        <stop offset="100%" style="stop-color:#059669"/>
                    </linearGradient>
                </defs>
            </svg>
        </div>
        <div>
            <h1 style="margin:0; font-size:1.8rem; font-weight:700;">FractoMap</h1>
            <p style="margin:0; opacity:0.9; font-size:0.95rem;">Bioactivity-guided microfractionation analysis</p>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.markdown("### ⚙️ Parameters")
    start = st.number_input("Collection start (min)", value=1.0, step=0.1)
    interval = st.number_input("Interval (sec)", value=7, step=1)
    offset = st.number_input("Fraction offset", value=-1, step=1)
    params = {'start': start, 'interval': interval, 'offset': offset}
    
    st.markdown("---")
    st.markdown("### 📖 About")
    st.markdown("FractoMap overlays bioactivity data with LC-MS chromatograms for compound identification.")
    st.markdown("[GitHub](https://github.com/LOOKiemycin/FractoMap)")

# Main Tabs
tab1, tab2, tab3, tab4 = st.tabs(["📤 Upload", "🧫 Plate", "📊 Results", "📈 Overlay"])

# Tab 1: Upload
with tab1:
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        <div class="upload-card">
            <div class="card-title">📁 1. Plate Data (Excel)</div>
            <div class="card-subtitle">96-well plate absorbance (8×12)</div>
        </div>
        """, unsafe_allow_html=True)
        pf = st.file_uploader("Upload plate data", type=['xlsx', 'xls', 'csv'], label_visibility="collapsed")
        if pf:
            p = parse_plate_excel(pf)
            if p: 
                st.session_state.plate = p
                st.success("✅ Plate data loaded")
            else: 
                st.error("❌ Invalid format")
    
    with col2:
        st.markdown("""
        <div class="upload-card">
            <div class="card-title">📈 2. MS Data (mzML/CSV)</div>
            <div class="card-subtitle"><span class="tag tag-teal">TIC</span><span class="tag tag-blue">BPC</span></div>
        </div>
        """, unsafe_allow_html=True)
        cf = st.file_uploader("Upload chromatogram", type=['mzML', 'mzml', 'csv'], label_visibility="collapsed")
        if cf:
            if cf.name.lower().endswith('.mzml'):
                st.session_state.chrom = parse_mzml(cf.read().decode('utf-8'))
            else:
                d = pd.read_csv(cf)
                st.session_state.chrom['tic'] = list(zip(d.iloc[:,0], d.iloc[:,1]))
            st.success(f"✅ Chromatogram loaded")
        
        if st.button("🎲 Demo Data", use_container_width=True):
            st.session_state.chrom = generate_demo()
            st.success("✅ Demo loaded")
    
    with col3:
        st.markdown("""
        <div class="upload-card">
            <div class="card-title">🏷️ 3. GNPS2 Annotation</div>
            <div class="card-subtitle">Compound IDs (optional)</div>
        </div>
        """, unsafe_allow_html=True)
        gf = st.file_uploader("Upload GNPS2 results", type=['xlsx', 'xls', 'csv', 'tsv'], label_visibility="collapsed", key="gnps_upload")
        if gf:
            gnps_df = parse_gnps_excel(gf)
            if gnps_df is not None:
                st.session_state.gnps = gnps_df
                st.success(f"✅ {len(gnps_df)} annotations loaded")
            else:
                st.error("❌ Cannot parse GNPS file")
    
    st.markdown("---")
    
    if st.button("🧮 Calculate Inhibition", type="primary", use_container_width=True):
        if st.session_state.plate:
            df, ctrl = calculate_inhibition(st.session_state.plate, params)
            st.session_state.df, st.session_state.ctrl = df, ctrl
            st.success("✅ Calculation complete!")
            
            # Metrics
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Control Avg", f"{ctrl:.4f}")
            c2.metric("Max Inhibition", f"{df['% Inhibition'].max():.1f}%")
            c3.metric("Active (>50%)", len(df[df['% Inhibition'] > 50]))
            c4.metric("Best Fraction", f"F{df.loc[df['% Inhibition'].idxmax(), 'Fraction']}")
        else:
            st.error("❌ Please upload plate data first")

# Tab 2: Plate
with tab2:
    st.markdown("### 🧫 96-Well Plate Heatmap")
    if st.session_state.plate:
        fig = go.Figure(go.Heatmap(
            z=np.array(st.session_state.plate), 
            x=list(range(1,13)), 
            y=ROWS, 
            colorscale='RdYlGn_r',
            colorbar=dict(title="Abs")
        ))
        fig.update_layout(
            yaxis=dict(autorange='reversed'), 
            height=400,
            xaxis_title="Column",
            yaxis_title="Row"
        )
        st.plotly_chart(fig, use_container_width=True)
        
        st.info("🟢 Low absorbance = High inhibition | 🔴 High absorbance = Low inhibition")
    else:
        st.warning("⚠️ Upload plate data first")

# Tab 3: Results
with tab3:
    st.markdown("### 📊 Inhibition Results")
    if st.session_state.df is not None:
        df = st.session_state.df
        
        # Metrics row
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Control Avg", f"{st.session_state.ctrl:.4f}")
        c2.metric("Max Inhibition", f"{df['% Inhibition'].max():.1f}%")
        c3.metric("Active (>50%)", len(df[df['% Inhibition'] > 50]))
        c4.metric("Best Fraction", f"F{df.loc[df['% Inhibition'].idxmax(), 'Fraction']}")
        
        # Bar chart
        fig = go.Figure(go.Bar(
            x=df['Fraction'], 
            y=df['% Inhibition'], 
            marker_color=[get_color(i) for i in df['% Inhibition']]
        ))
        fig.update_layout(
            yaxis_range=[0, 100], 
            height=350,
            xaxis_title="Fraction",
            yaxis_title="% Inhibition"
        )
        fig.add_hline(y=50, line_dash="dash", line_color="crimson")
        st.plotly_chart(fig, use_container_width=True)
        
        # Data table
        st.dataframe(df, use_container_width=True, height=300)
        st.download_button("📥 Download CSV", df.to_csv(index=False), "fractomap_results.csv")
    else:
        st.warning("⚠️ Calculate inhibition first")

# Tab 4: Overlay
with tab4:
    st.markdown("### 📈 Bioactivity-Chromatogram Overlay")
    
    c1, c2, c3 = st.columns([1, 1, 2])
    with c1: 
        chrom_type = st.selectbox("Chromatogram", ["TIC", "BPC"])
    with c2: 
        bar_width = st.slider("Bar width", 0.05, 0.30, 0.12, 0.01)
    
    df, chrom = st.session_state.df, st.session_state.chrom.get(chrom_type.lower(), [])
    gnps_df = st.session_state.gnps
    
    if df is not None and chrom:
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        
        # TIC - light purple fill
        fig.add_trace(go.Scatter(
            x=[p[0] for p in chrom], y=[p[1] for p in chrom], 
            name=chrom_type, 
            line=dict(color='rgba(100,100,180,0.8)', width=1.5), 
            fill='tozeroy', 
            fillcolor='rgba(180,180,220,0.4)'
        ), secondary_y=True)
        
        # Separate traces for each activity level
        for act, color, label in [
            ('Strong', '#2E8B57', 'Strong (>75%)'),
            ('Moderate', '#FFD700', 'Moderate (50-75%)'),
            ('Weak', '#FF8C00', 'Weak (25-50%)'),
            ('Inactive', '#DC143C', 'Inactive (<25%)')
        ]:
            subset = df[df['Activity'] == act]
            if len(subset) > 0:
                fig.add_trace(go.Bar(
                    x=subset['RT (min)'], y=subset['% Inhibition'],
                    name=label, marker_color=color, opacity=0.85, width=bar_width,
                    text=[f"F{f}" if inh > 50 else "" for f, inh in zip(subset['Fraction'], subset['% Inhibition'])],
                    textposition='outside', textfont=dict(size=9, color='#333')
                ), secondary_y=False)
        
        # Add GNPS compound labels
        if gnps_df is not None and len(gnps_df) > 0:
            # Get max intensity for positioning
            max_int = max([p[1] for p in chrom]) if chrom else 1
            
            for _, row in gnps_df.iterrows():
                rt = row['RT']
                name = row['Compound']
                # Shorten long names
                if len(name) > 20:
                    name = name[:18] + "..."
                
                # Add annotation
                fig.add_annotation(
                    x=rt, y=max_int * 0.95,
                    text=f"<b>{name}</b>",
                    showarrow=True,
                    arrowhead=2,
                    arrowsize=0.8,
                    arrowwidth=1,
                    arrowcolor="#6366F1",
                    ax=0, ay=-30,
                    font=dict(size=10, color="#4F46E5"),
                    bgcolor="rgba(255,255,255,0.9)",
                    bordercolor="#6366F1",
                    borderwidth=1,
                    borderpad=3,
                    xref="x", yref="y2"
                )
        
        # 50% threshold line
        fig.add_hline(y=50, line_dash="dash", line_color="crimson", line_width=1.5, secondary_y=False)
        
        fig.update_layout(
            title=dict(
                text=f"<b>Bioactivity Overlay with {chrom_type}</b><br><sup>Offset = {params['offset']} fraction (~{abs(params['offset'])*7} sec dead volume)</sup>", 
                x=0.5, xanchor='center'
            ),
            xaxis_title="Retention Time (min)",
            height=600, 
            bargap=0, 
            barmode='overlay',
            legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99, bgcolor='rgba(255,255,255,0.9)', bordercolor='#ccc', borderwidth=1),
            plot_bgcolor='white'
        )
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='#eee')
        fig.update_yaxes(title_text="% Inhibition", range=[0, 105], showgrid=True, gridcolor='#eee', 
                        title_font=dict(color='#2E8B57'), tickfont=dict(color='#2E8B57'), secondary_y=False)
        fig.update_yaxes(title_text=f"{chrom_type} Intensity", showgrid=False,
                        title_font=dict(color='#6464B4'), tickfont=dict(color='#6464B4'), secondary_y=True)
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Active fractions summary
        active = df[df['% Inhibition'] > 50]
        if len(active) > 0:
            st.success(f"**Active fractions:** {', '.join([f'F{f}' for f in active['Fraction'].values])} (RT {active['RT (min)'].min():.2f} - {active['RT (min)'].max():.2f} min)")
        
        # Show GNPS annotations in active region
        if gnps_df is not None and len(active) > 0:
            rt_min, rt_max = active['RT (min)'].min() - 0.5, active['RT (min)'].max() + 0.5
            active_compounds = gnps_df[(gnps_df['RT'] >= rt_min) & (gnps_df['RT'] <= rt_max)]
            if len(active_compounds) > 0:
                st.info(f"**🏷️ Compounds in active region:** {', '.join(active_compounds['Compound'].tolist())}")
    else:
        st.warning("⚠️ Upload data and calculate inhibition first")

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #6b7280; font-size: 0.9rem;">
    <b>FractoMap</b> • Bioactivity-Guided Microfractionation Analysis<br>
    Developed by Thapanee Pruksatrakul • <a href="https://github.com/LOOKiemycin/FractoMap">GitHub</a>
</div>
""", unsafe_allow_html=True)
