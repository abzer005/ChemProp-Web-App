import streamlit as st
from src.common import *
import pandas as pd
from streamlit.components.v1 import html

page_setup()

# Add a tracking token
html('<script async defer data-website-id="74bc9983-13c4-4da0-89ae-b78209c13aaf" src="https://analytics.gnps2.org/umami.js"></script>', width=0, height=0)

try:
    st.image("assets/chemprop2_logo.png", caption="Chemprop2 Logo", use_container_width=True)
except TypeError:
    st.image("assets/chemprop2_logo.png", caption="Chemprop2 Logo", use_column_width=True)

# Introduction
st.markdown("""
<div style="color: #05577C; font-size: 25px;">
<b>Welcome to ChemProp2!</b> Before proceeding with your analysis, please take a moment to read this homepage.
</div>
""", unsafe_allow_html=True)

st.write(' ')

st.subheader("What is ChemProp2?")

st.write("""
**ChemProp2** is a network-based metabolomics tool for identifying and prioritizing 
**putative biotransformations** in non-targeted LC–MS/MS data. It is particularly suited for **time-series and multi-condition experiments**, such as:
- Microbiome-mediated drug metabolism
- Xenobiotic and environmental biotransformation
- Natural product transformation studies

ChemProp2 builds on **Feature-based Molecular Networking (FBMN)** and uses 
**correlation-based scoring** to detect anti-correlating metabolite pairs and 
infer transformation relationships across networks.

### ChemProp1 vs ChemProp2 (Quick Guide)

- **ChemProp1**  
  - Log-ratio–based approach  
  - Best for **two-point comparisons**  
  - [Published method](https://doi.org/10.1021/acs.analchem.1c01520). Available in this app for convenience

- **ChemProp2**  
  - Correlation-based approach  
  - Best for **three or more sequential data points**  
""")

# Input Files
st.subheader("Input Data Options")

st.write("""
ChemProp2 supports **three ways** to provide input data. Choose the option that best fits your workflow:
""")

with st.expander("🧪 Example Dataset (Quick Start)"):
    st.markdown("""
This option loads a **lightweight example dataset** to help you quickly explore the app.

- Contains a **subset** of an FBMN task used in the ChemProp2 manuscript
- Optimized for **fast execution**
- The drug-related feature to inspect is **`drug1`** (see annotation table)

Recommended if you want to **try the app immediately** without preparing input files.
""")

with st.expander("🧬 FBMN Task ID (GNPS2)"):
    st.markdown("""
Provide an **FBMN Task ID from GNPS2** to automatically retrieve all required inputs.

- The default Task ID corresponds to the **ChemProp2 manuscript dataset**  
  (link will be added soon)
- Dataset description:
  - Synthetic gut bacterial community (SynCom)
  - Treated with **12 drugs**
  - Samples collected across **multiple timepoints**

You may use this dataset as a **fully worked example** to explore ChemProp scoring (both ChemProp1 and ChemProp2 modes),
cascade analysis, and network-level visualization.
""")

with st.expander("📂 Manual Upload"):
    st.markdown("""
Upload your own data in **CSV**, **TXT**, or **XLSX** format.

Required / optional inputs include:

- **Feature Quantification Table** (required)
- **Metadata Table** (required)
- **Edge Table** from FBMN or user-generated networks  
  (required for 'Provided Edge Table' Mode and 'Cascade Analysis' Mode)
- **Annotation Tables** from FBMN (optional but recommended)

This option provides full flexibility for **custom datasets** and non-GNPS workflows.
""")

st.markdown("""
### 🔗 External Interpretation Tools

Integration with **Spectrum Resolver** and **ModiFinder** is available **only when using the FBMN Task ID option**.

These tools require **Universal Spectrum Identifiers (USIs)**, which are automatically available
when data are retrieved directly from GNPS2.

- Spectrum Resolver: https://metabolomics-usi.ucsd.edu
- ModiFinder: https://modifinder.gnps2.org
""")

st.warning(
    "Manual uploads **do not support** Spectrum Resolver or ModiFinder "
    "due to the absence of Universal Spectrum Identifiers (USIs)."
)

st.write("""
### Data Preparation Essentials

To ensure smooth and error-free processing, please follow these guidelines:

- Input files must include the `.mzML` extension in both the feature quantification table and the metadata table.
- **Metadata Table (required):**
  - Must include a column named **filename**, matching the feature table exactly.
  - May include additional attribute columns such as **replicates** or **sample type** (e.g., control, treatment).
  - A **Time column is mandatory** and **must contain numeric values** representing time points.  
    (e.g., `2`, `2h`, `2hr`, `2min`; text will be stripped automatically).
""")

with st.expander("📊 Example Input Tables"):
    st.markdown("""

**Example Feature Quantification Table**

| feature | sample1.mzML | sample2.mzML | blank.mzML |
|--------|---------------|---------------|------------|
| 1 | 1000 | 1100 | 100 |
| 2 | 2000 | 2200 | 200 |
""")

    st.markdown("""
**Example Metadata Table**

| filename | Sample_Type | Time |
|----------|-------------|------|
| sample1.mzML | Sample | 1h |
| sample2.mzML | Sample | 2h |
| blank.mzML | Blank | N/A |
""")

# Define the data
data = {
    "File Type": [
        "Feature Quantification Table",
        "Metadata",
        "Edge File",
        "Annotation Files",
        "Annotation Files"
    ],
    "Folder Location": [
        "`quantification_table` folder",
        "`metadata_table` folder",
        "`networking_pairs_results_file_filtered` folder",
        "`DB_result` folder",
        "`DB_analogresult` folder"
    ],
    "Description": [
        "",
        "",
        "",
        "Includes library annotations for features.",
        "Includes both library annotations and analog hits for features."
    ]
}

# Convert to a pandas DataFrame
df = pd.DataFrame(data)

data = {
    "ChemProp Input File": [
        "Feature Quantification Table",
        "Metadata Table",
        "Edge Table",
        "Annotation Table (Library Matches)",
        "Annotation Table (Library + Analog Matches)"
    ],
    "GNPS1 Folder Location": [
        "`quantification_table`",
        "`metadata_table`",
        "`networking_pairs_results_file_filtered`",
        "`DB_result`",
        "`DB_analogresult`"
    ],
    "Notes": [
        "Required for all ChemProp analyses.",
        "Must include filename and time information.",
        "Required for network-based and cascade analyses.",
        "Optional. Provides GNPS library annotations.",
        "Optional. Includes both library annotations and analog hits."
    ]
}

df = pd.DataFrame(data)


# Display the table
st.write(' ')
with st.expander("📌 For GNPS1 Users: Uploading Data"):
    st.markdown("""
⚠️ **GNPS1 Task IDs are not supported in this app.**

- This applies to **both ChemProp modules** available here (ChemProp1 and ChemProp2).
- Only **FBMN Task IDs from GNPS2** are supported. 
- If your data were generated using **GNPS1**, please do **not** enter a GNPS1 Task ID.

**What should GNPS1 users do instead?**

1. Navigate to your **GNPS1 FBMN task** on GNPS.
2. Go to the **Job Status** page.
3. Under **Export / Download Network Files**, click **Download Cytoscape Data**.
4. Unzip the downloaded folder.

Within the extracted folder, you will find the files required for **Manual Upload**, as summarized below:
""")

    st.dataframe(df, hide_index=True, use_container_width=True)


st.subheader("Choosing the Right Mode for Your Analysis")
with st.expander("🧭 **DECISION GUIDE**: Which mode should I use?"):
    st.markdown("""
Use this quick guide to choose the best mode for your question:

**✅ Use _Provided Edge Table mode_ if you want to…**
- Score transformations **across a whole network**
- Reproduce or extend results from **FBMN**
- Filter / prioritize edges based on ChemProp scores

**✅ Use _Cascade edges mode_ if you want to…**
- Start from **one feature (e.g., a drug node)** and explore downstream / upstream relationships
- Identify **multi-step transformation chains**
- Prioritize **higher-order** candidates beyond first-degree neighbors

**✅ Use _User defined edge mode_ if you want to…**
- Test a **specific hypothesis** (Feature A → Feature B)
- Compare two selected nodes even if they are not strongly connected in the provided network
- Quickly validate a suspected transformation between two features

**Rule of thumb:**  
- *Network-wide scoring* → **Provided Edge Table**  
- *Transformation chains from a starting node* → **Cascade Edges**  
- *One targeted comparison* → **User defined edge**
""")

st.subheader("Output File Information")

st.write("""
ChemProp2 generates a single **CSV output file** containing **ChemProp2 scores for each node pair** in the molecular network.
The output extends the input edge (node-pair) table by adding ChemProp2 scoring information for each node pair.

**Key details:**
- **Score range:** −1 to +1 for each node pair
- **Score magnitude:** indicates the strength of the putative transformation
- **Score sign:** + or - sign indicates transformation directionality  
  (e.g., whether A → B or B → A within a node pair)

Each row in the CSV corresponds to one edge in the molecular network.
""")

st.write("""
### Integration with Cytoscape
- You can download the results as a **zip file** containing the required **GraphML** and **style files** for Cytoscape visualization.
- To use the output in Cytoscape:
  1. Open the **GraphML** file in Cytoscape.
  2. Import the **style file**:
     - Navigate to **File > Import > Styles from File**.
     - Select and upload the `.xml` style file.
  3. Apply the imported style:
     - Go to the **Styles** panel in Cytoscape.  
     - By default, the style will be set to **default**.  
     - Use the drop-down menu to select the uploaded style (it may appear as **default_0** or another variation).  
     - The visualization will update based on the selected style.
        
         """)
  
st.write("""
### Additional Visualizations and Features
- **Global Transformation Patterns**: Visualize overall transformation trends in your data with a score plot showing ChemProp scores against various mass differences.
- **Filtering Options**: Use filters to refine scores by range, mass, annotation name, or cluster ID, displaying only relevant edges in a list. Selecting an edge shows:
  - The intensity trends for node pairs belonging to that edge.
  - The corresponding node pair within the molecular network as a subnetwork. (Note: Full network visualization is minimized for simplified rendering, showing only the relevant cluster and its ChemProp2 scores with directional indicators.)
  - **Export Visualizations**: Save these visualizations as PNG files for easy sharing and documentation.

""")

# Subheader and Interactive Features
st.subheader('About the App Elements')
st.markdown("""
🔍 **How to Know If the App Is Running?**  
If you're performing a calculation or generating a figure and don't see any updates on the main page, 
check the **top-right corner** of the page. If you see the message **'RUNNING'**, the app is active and processing your request.  
            
💡 **All plots are interactive!**  
- Use your mouse to select areas and zoom in on specific regions.  
- Double-click on the plot to reset the zoom and return to the default view.  
- Save plots using the **camera icon** in the top-right corner of the plot. You can specify the image format (e.g., PNG) in the settings panel before saving.
""")



st.subheader("⬇️ Download the ChemProp2 App")
# ----------------------------
# Windows expander
# ----------------------------
with st.expander("🪟 Windows users", expanded=False):
    st.markdown("""
### ✅ Option 1 — Website download (recommended)
1. Go to: https://www.functional-metabolomics.com/resources  
2. Under **Software**, click **Download** next to **ChemProp2 Web App**

### ⚙️ Option 2 — GitHub Actions build (advanced)
> Requires a GitHub account

1. Go to: https://github.com/Functional-Metabolomics-Lab/ChemProp-Web-App/actions  
2. Select the **latest successful workflow run**
3. Scroll to **Artifacts** and **download** the ZIP file  
4. **Extract** the ZIP locally  
5. Run the **ChemProp2 app executable** inside the extracted folder
""")
# ----------------------------
# macOS expander
# ----------------------------
with st.expander("🍎 macOS users", expanded=False):
    st.markdown("""
🚫 A pre-built macOS version is **not available**.

### ▶️ Run locally from source
For **macOS** and **Windows** users who want full control over the app.

**Step-by-step**
1. Clone the repository:  
   ``git clone https://github.com/Functional-Metabolomics-Lab/ChemProp-Web-App.git``

2. Navigate into the project folder:  
   ``cd ChemProp-Web-App``

3. Install required dependencies:  
   ``pip install -r requirements.txt``

4. Launch the app:  
   ``streamlit run Home.py``

The app will open automatically in your web browser at a local address  
(for example: ``http://localhost:8501``).
                """)

# ----------------------------
# Optional note expander
# ----------------------------
with st.expander("ℹ️ Notes & tips", expanded=False):
 st.markdown("""
- The **GitHub Actions build** gives more control and transparency  
- Running from source requires **Python + Streamlit** installed locally
- Make sure you have Python 3.11 installed (same version used in the Windows .exe build).
- You can run Online (Recommended for smaller networks). If your input data is relatively large, please consider running it locally.

""")
                

# Citation and Resources
st.subheader('Citation and Further Resources')
st.write('If you use ChemProp2 in your research, please cite:')
st.markdown("""
            * [FBMN-STATS](https://fbmn-statsguide.gnps2.org/) - A statistical pipeline for downstream processing of FBMN results.
            * Pakkir Shah, A.K., Walter, A., Ottosson, F. et al. Statistical analysis of feature-based molecular networking results from non-targeted metabolomics data. Nat Protoc (2024). https://doi.org/10.1038/s41596-024-01046-3
            """
            )
            

# Add more links as needed

# Feedback Section
st.subheader("We Value Your Feedback")
st.markdown("""
            We welcome your feedback and suggestions to improve ChemProp2. Please feel free to create an issue on our GitHub repository to share your thoughts or report any issues you encounter. 
            Your input is invaluable in making ChemProp2 better for everyone.

            [Create an Issue on GitHub](https://github.com/abzer005/ChemProp2/issues/new)
""")

# Contribution and Follow Us
st.subheader("Contribute and Follow Us")
st.markdown("""
- Interested in contributing? Check out the [GitHub page](https://github.com/abzer005).
- For more about our work, visit our [lab's GitHub page](https://github.com/Functional-Metabolomics-Lab).
""")

# Optional: Footer
st.markdown("---")
st.text("ChemProp2 © 2023")