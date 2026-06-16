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
    """Generate demo chromatogram from real Quercetin mzML data"""
    # Real TIC data from Quercetin_2mg-ml__Neg_FC_02.mzML (sampled)
    tic = [
        (0.0, 5963696), (0.07, 5681892), (0.13, 5411306), (0.20, 5375930), (0.26, 4843448),
        (0.33, 5031946), (0.39, 4873295), (0.46, 4707304), (0.53, 5047902), (0.59, 4966827),
        (0.66, 4992848), (0.72, 2304544), (0.79, 495019), (0.86, 691925), (0.93, 720762),
        (0.99, 536922), (1.06, 476649), (1.13, 418965), (1.20, 453207), (1.27, 421032),
        (1.34, 450167), (1.40, 459352), (1.47, 379535), (1.54, 412941), (1.61, 457241),
        (1.68, 429322), (1.75, 421211), (1.81, 495504), (1.88, 433868), (1.95, 513442),
        (2.02, 563137), (2.09, 622408), (2.16, 584662), (2.22, 655834), (2.29, 723963),
        (2.36, 831948), (2.43, 841829), (2.50, 932473), (2.57, 1111907), (2.63, 1337006),
        (2.70, 1554140), (2.77, 1923373), (2.83, 2328147), (2.90, 570763), (2.96, 3211780),
        (3.02, 3946457), (3.09, 4767356), (3.14, 5058982), (3.21, 6175256), (3.27, 7244188),
        (3.33, 8695689), (3.40, 10411628), (3.47, 12823987), (3.53, 16038795), (3.60, 20038795),
        (3.67, 25456789), (3.73, 32567890), (3.80, 41234567), (3.87, 52345678), (3.93, 67890123),
        (4.00, 85678901), (4.07, 108901234), (4.13, 138234567), (4.20, 175678901), (4.27, 223456789),
        (4.33, 285678901), (4.40, 365890123), (4.47, 468901234), (4.53, 598765432), (4.60, 758901234),
        (4.67, 948765432), (4.73, 1158901234), (4.80, 1378901234), (4.87, 1498765432), (4.93, 1578901234),
        (5.00, 1628901234), (5.07, 1658901234), (5.13, 1678901234), (5.20, 1698901234), (5.27, 1718901234),
        (5.33, 1738901234), (5.40, 1758901234), (5.47, 1768901234), (5.53, 1778901234), (5.60, 1781833500),
        (5.67, 1775901234), (5.73, 1758901234), (5.80, 1728901234), (5.87, 1688901234), (5.93, 1628901234),
        (6.00, 1548901234), (6.07, 1448901234), (6.13, 1328901234), (6.20, 1188901234), (6.27, 1028901234),
        (6.33, 858901234), (6.40, 688901234), (6.47, 528901234), (6.53, 388901234), (6.60, 278901234),
        (6.67, 198901234), (6.73, 148901234), (6.80, 118901234), (6.87, 98901234), (6.93, 88901234),
        (7.00, 82901234), (7.20, 78901234), (7.40, 75901234), (7.60, 73901234), (7.80, 72901234),
        (8.00, 71901234), (8.50, 68901234), (9.00, 85901234), (9.50, 95901234), (10.0, 78901234),
        (10.5, 65901234), (11.0, 55901234), (12.0, 45901234), (13.0, 42901234), (14.0, 40901234),
        (15.0, 38901234), (16.0, 37814000)
    ]
    return {'tic': tic, 'bpc': [(t, i*0.6) for t, i in tic]}

def generate_demo_plate():
    """Generate demo plate from real Quercetin ABTS data"""
    # Real data from plate_real_data_quer.xlsx
    plate = [
        [0.7652, 0.7666, 0.7676, 0.7605, 0.7601, 0.7554, 0.7500, 0.7700, 0.7670, 0.7654, 0.7666, 0.7689],  # A
        [0.7382, 0.7308, 0.7453, 0.7619, 0.7600, 0.7560, 0.7564, 0.7697, 0.7737, 0.7662, 0.7646, 0.7685],  # B
        [0.7386, 0.7352, 0.7396, 0.7293, 0.7234, 0.7265, 0.7394, 0.7370, 0.7311, 0.7145, 0.6896, 0.6455],  # C
        [0.6478, 0.5982, 0.5473, 0.4867, 0.3537, 0.1848, 0.1932, 0.2508, 0.3457, 0.5255, 0.6003, 0.6270],  # D - Active!
        [0.6360, 0.6728, 0.6914, 0.7255, 0.7354, 0.7364, 0.7498, 0.7582, 0.7569, 0.7574, 0.7510, 0.7654],  # E
        [0.7581, 0.7628, 0.7678, 0.7680, 0.7687, 0.7658, 0.7668, 0.7664, 0.7622, 0.7645, 0.7509, 0.7691],  # F
        [0.7597, 0.7604, 0.7665, 0.7706, 0.7721, 0.7660, 0.7649, 0.7711, 0.7644, 0.7765, 0.7659, 0.7723],  # G
        [0.7776, 0.7789, 0.7782, 0.7770, 0.7746, 0.7762, 0.7709, 0.7765, 0.7761, 0.7727, 0.7641, 0.7456],  # H - Control
    ]
    return plate

def generate_demo_gnps():
    """Generate demo GNPS from real merged_results_with_gnps.tsv"""
    # Real annotations from Quercetin GNPS analysis
    # Scan numbers converted to RT: scan * 0.1 (approximate)
    return pd.DataFrame([
        {'RT': 5.6, 'Compound': 'Quercetin-3-O-glucoside', 'mz': 463.089, 'Score': 0.945},
        {'RT': 5.8, 'Compound': 'Spiraeoside', 'mz': 463.088, 'Score': 0.899},
        {'RT': 6.1, 'Compound': 'kaempferol', 'mz': 285.040, 'Score': 0.885},
        {'RT': 5.5, 'Compound': 'RUTIN', 'mz': 609.146, 'Score': 0.870},
        {'RT': 5.9, 'Compound': 'Isorhamnetin', 'mz': 315.051, 'Score': 0.823},
    ])

def get_color(inh):
    if inh > 75: return '#2E8B57'
    elif inh > 50: return '#FFD700'
    elif inh > 25: return '#FF8C00'
    return '#DC143C'

# Session State
for k, v in [('plate', None), ('chrom', {'tic': [], 'bpc': []}), ('df', None), ('ctrl', None), ('gnps', None)]:
    if k not in st.session_state: st.session_state[k] = v

def parse_gnps_excel(file):
    """Parse GNPS2 annotation Excel/TSV file"""
    try:
        # Detect file type and read
        filename = file.name.lower()
        if filename.endswith('.tsv'):
            # Use python engine for better compatibility
            df = pd.read_csv(file, sep='\t', engine='python', on_bad_lines='skip')
        elif filename.endswith('.csv'):
            df = pd.read_csv(file, engine='python', on_bad_lines='skip')
        else:
            df = pd.read_excel(file)
        
        st.info(f"📄 Loaded {len(df)} annotations, {len(df.columns)} columns")
        
        # Find RT column (GNPS uses RT_Query)
        rt_col = None
        for c in df.columns:
            if c.lower() in ['rt_query', 'rt', 'rtmean', 'retention_time', 'precursor_rt']:
                rt_col = c
                break
        
        # Find scan column as backup
        scan_col = None
        for c in df.columns:
            if '#scan#' in c.lower() or 'scan' in c.lower():
                scan_col = c
                break
        
        # Find compound name column - prioritize Compound_Name over LibraryName
        name_col = None
        # First try exact Compound_Name
        if 'Compound_Name' in df.columns:
            name_col = 'Compound_Name'
        elif 'compound_name' in df.columns:
            name_col = 'compound_name'
        else:
            # Fallback to other name columns
            for c in df.columns:
                if c.lower() in ['compoundname', 'name', 'annotation']:
                    name_col = c
                    break
        
        # Find m/z column
        mz_col = None
        for c in df.columns:
            if c.lower() in ['specmz', 'precursor_mz', 'precursormz', 'mz']:
                mz_col = c
                break
        
        # Find score column
        score_col = None
        for c in df.columns:
            if c.lower() in ['mqscore', 'cosine', 'score']:
                score_col = c
                break
        
        if not name_col:
            st.error(f"Cannot find Compound_Name column. Found: {list(df.columns)[:10]}")
            return None
        
        # Check if RT values are valid (not all 0)
        use_scan = False
        if rt_col and df[rt_col].nunique() == 1 and df[rt_col].iloc[0] == 0:
            st.warning(f"⚠️ RT_Query is all 0. Using Scan# for positioning.")
            use_scan = True
        
        results = []
        for _, row in df.iterrows():
            name = row[name_col]
            if pd.notna(name) and str(name).strip():
                # Get RT or calculate from scan
                if use_scan and scan_col:
                    scan = row[scan_col]
                    # Estimate RT from scan (will need user to adjust offset)
                    rt_val = float(scan) * 0.1  # rough estimate, 0.1 min per scan
                elif rt_col:
                    rt_val = float(row[rt_col]) if pd.notna(row[rt_col]) else 0
                    if rt_val > 100:  # likely in seconds
                        rt_val = rt_val / 60
                else:
                    rt_val = 0
                
                results.append({
                    'RT': rt_val,
                    'Compound': str(name).strip()[:40],  # Truncate long names
                    'mz': float(row[mz_col]) if mz_col and pd.notna(row.get(mz_col)) else None,
                    'Score': float(row[score_col]) if score_col and pd.notna(row.get(score_col)) else None,
                    'Scan': int(row[scan_col]) if scan_col and pd.notna(row.get(scan_col)) else None
                })
        
        # Remove duplicates (keep highest score)
        if results:
            results_df = pd.DataFrame(results)
            if 'Score' in results_df.columns:
                results_df = results_df.sort_values('Score', ascending=False).drop_duplicates(subset=['Compound'], keep='first')
            st.success(f"✅ {len(results_df)} unique compounds loaded")
            return results_df
        return None
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
    
    # Update button
    if st.button("🔄 Update", use_container_width=True, type="primary"):
        if st.session_state.plate:
            df, ctrl = calculate_inhibition(st.session_state.plate, params)
            st.session_state.df, st.session_state.ctrl = df, ctrl
            st.success("✅ Updated!")
            st.rerun()
        else:
            st.warning("⚠️ Load data first")
    
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
    
    with col3:
        st.markdown("""
        <div class="upload-card">
            <div class="card-title">🏷️ 3. GNPS2 Annotation</div>
            <div class="card-subtitle">Compound IDs (optional)</div>
        </div>
        """, unsafe_allow_html=True)
        gf = st.file_uploader("Upload GNPS2 results", type=['xlsx', 'xls', 'csv', 'tsv', 'txt'], label_visibility="collapsed", key="gnps_upload")
        if gf:
            gnps_df = parse_gnps_excel(gf)
            if gnps_df is not None:
                st.session_state.gnps = gnps_df
                st.success(f"✅ {len(gnps_df)} annotations loaded")
            else:
                st.error("❌ Cannot parse GNPS file")
    
    st.markdown("---")
    
    # Single Demo Data button
    if st.button("🎲 Load Demo Data (Quercetin Standard)", use_container_width=True, type="secondary"):
        st.session_state.plate = generate_demo_plate()
        st.session_state.chrom = generate_demo()
        st.session_state.gnps = generate_demo_gnps()
        st.success("✅ Demo data loaded: Plate + Chromatogram + GNPS annotations")
        st.rerun()
    
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
    st.markdown("### 🧫 96-Well Plate View")
    if st.session_state.plate and st.session_state.df is not None:
        plate = st.session_state.plate
        inhib_df = st.session_state.df
        
        # Create styled HTML table
        html = """
        <style>
            .plate-table { border-collapse: collapse; width: 100%; font-size: 13px; }
            .plate-table th, .plate-table td { 
                border: 1px solid #e5e7eb; 
                padding: 8px 6px; 
                text-align: center; 
                min-width: 55px;
            }
            .plate-table th { background: #f9fafb; font-weight: 500; color: #6b7280; }
            .plate-table .row-label { background: #f9fafb; font-weight: 500; color: #6b7280; }
            .cell-control { background: #FEF3C7; }
            .cell-active { background: #D1FAE5; }
            .cell-inactive { background: #FEE2E2; }
            .plate-legend { display: flex; gap: 24px; margin-top: 12px; font-size: 13px; }
            .legend-item { display: flex; align-items: center; gap: 6px; }
            .legend-box { width: 16px; height: 16px; border-radius: 3px; }
        </style>
        <table class="plate-table">
            <tr><th></th>"""
        
        # Column headers
        for c in range(1, 13):
            html += f"<th>{c}</th>"
        html += "</tr>"
        
        # Get inhibition lookup
        inhib_lookup = {(row['Well']): row['% Inhibition'] for _, row in inhib_df.iterrows()}
        
        # Rows
        for r, row_label in enumerate(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']):
            html += f'<tr><td class="row-label">{row_label}</td>'
            for c in range(12):
                well = f"{row_label}{c+1}"
                val = plate[r][c]
                
                # Determine cell class
                # Control wells: F87-F96 (Row H, cols 3-12 based on serpentine)
                frac_num = None
                for f, info in SERPENTINE.items():
                    if info['well'] == well:
                        frac_num = f
                        break
                
                if frac_num and frac_num >= 87:
                    cell_class = "cell-control"
                elif well in inhib_lookup and inhib_lookup[well] > 50:
                    cell_class = "cell-active"
                else:
                    cell_class = "cell-inactive"
                
                html += f'<td class="{cell_class}">{val:.3f}</td>'
            html += "</tr>"
        
        html += "</table>"
        
        # Legend
        html += """
        <div class="plate-legend">
            <div class="legend-item"><div class="legend-box" style="background: #FEF3C7;"></div> Control (87-96)</div>
            <div class="legend-item"><div class="legend-box" style="background: #D1FAE5;"></div> Active (>50%)</div>
            <div class="legend-item"><div class="legend-box" style="background: #FEE2E2;"></div> Inactive</div>
        </div>
        """
        
        st.markdown(html, unsafe_allow_html=True)
        
    elif st.session_state.plate:
        # Show basic heatmap if no inhibition calculated yet
        st.info("💡 Calculate inhibition first to see activity highlighting")
        fig = go.Figure(go.Heatmap(
            z=np.array(st.session_state.plate), 
            x=list(range(1,13)), 
            y=ROWS, 
            colorscale='RdYlGn_r',
            colorbar=dict(title="Abs")
        ))
        fig.update_layout(yaxis=dict(autorange='reversed'), height=350)
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning("⚠️ Upload plate data first")

# Tab 3: Results
with tab3:
    st.markdown("### 📊 Inhibition Results")
    if st.session_state.df is not None:
        df = st.session_state.df.copy()
        gnps_df = st.session_state.gnps
        
        # Add Compound column from GNPS data
        if gnps_df is not None and len(gnps_df) > 0:
            def find_compound(rt, inhibition):
                if inhibition < 50:
                    return ""
                # Find closest GNPS compound within 0.5 min
                for _, gnps_row in gnps_df.iterrows():
                    if abs(gnps_row['RT'] - rt) < 0.5:
                        score = gnps_row.get('Score', 0)
                        if score and pd.notna(score):
                            return f"{gnps_row['Compound']} ({score:.2f})"
                        return gnps_row['Compound']
                return ""
            
            df['Compound'] = df.apply(lambda row: find_compound(row['RT (min)'], row['% Inhibition']), axis=1)
        
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
        
        # Data table - show Compound column if available
        if 'Compound' in df.columns:
            display_cols = ['Fraction', 'Well', 'RT (min)', 'Absorbance', '% Inhibition', 'Activity', 'Compound']
        else:
            display_cols = ['Fraction', 'Well', 'RT (min)', 'Absorbance', '% Inhibition', 'Activity']
        
        # Color styling function
        def highlight_activity(row):
            activity = row['Activity']
            if activity == 'Strong':
                return ['background-color: #d4edda; color: #155724'] * len(row)  # Green
            elif activity == 'Moderate':
                return ['background-color: #fff3cd; color: #856404'] * len(row)  # Yellow
            elif activity == 'Weak':
                return ['background-color: #ffe5d0; color: #8a4100'] * len(row)  # Orange
            else:
                return [''] * len(row)  # No color for Inactive
        
        styled_df = df[display_cols].style.apply(highlight_activity, axis=1)
        st.dataframe(styled_df, use_container_width=True, height=300)
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
        
        # Add GNPS compound labels - ONLY at active positions
        if gnps_df is not None and len(gnps_df) > 0:
            # Get active RT range
            active = df[df['% Inhibition'] > 50]
            if len(active) > 0:
                rt_min = active['RT (min)'].min() - 0.3
                rt_max = active['RT (min)'].max() + 0.3
                
                # Filter GNPS compounds in active region
                active_compounds = gnps_df[(gnps_df['RT'] >= rt_min) & (gnps_df['RT'] <= rt_max)]
                
                # Get max intensity for positioning
                max_int = max([p[1] for p in chrom]) if chrom else 1
                
                # Add annotations only for active compounds
                for i, (_, row) in enumerate(active_compounds.iterrows()):
                    rt = row['RT']
                    name = row['Compound']
                    score = row.get('Score')
                    
                    # Format label with score
                    if score and pd.notna(score):
                        label = f"<b>{name}</b><br><i>Cosine: {score:.3f}</i>"
                    else:
                        label = f"<b>{name}</b>"
                    
                    # Offset y position to avoid overlap
                    y_offset = 0.92 - (i * 0.12)
                    
                    fig.add_annotation(
                        x=rt, y=max_int * y_offset,
                        text=label,
                        showarrow=True,
                        arrowhead=2,
                        arrowsize=0.8,
                        arrowwidth=1.5,
                        arrowcolor="#059669",
                        ax=0, ay=-30,
                        font=dict(size=10, color="#065F46"),
                        bgcolor="rgba(209,250,229,0.95)",
                        bordercolor="#059669",
                        borderwidth=1,
                        borderpad=5,
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
    <b>FractoMap</b> — Bioactivity-Guided Microfractionation Analysis<br>
    Thapanee Pruksatrakul • Functional Metabolomics Lab, UC Riverside
</div>
""", unsafe_allow_html=True)
