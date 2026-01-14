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
st.subheader('What is ChemProp2 Used For?')

st.write("""
         ChemProp2 is a tool developed to address a key challenge in non-targeted metabolomics using liquid chromatography-tandem mass spectrometry (LC–MS/MS), 
         particularly in the study of biotransformation studies like drug metabolism as well as xenobiotic and natural product biotransformation in the environment. 
        
         Analyzing and annotating the vast data from metabolomic studies still remains a challenge. Various in silico methods and spectral similarity metrics have been developed to tackle this issue. 
         Tools like GNPS (now GNPS2) use Feature-based Molecular Networking (FBMN) to create molecular networks by connecting metabolites with similar MS/MS spectral profiles. 
         ChemProp2 builds on this, identifying potential biotransformations within these networks. It detects anti-correlating metabolites and putatibve reaction pairs, scoring their correlation over time or space. 
         This helps in prioritizing and visualizing biochemical alterations within the network.
        
         ChemProp2 is particularly useful when dealing with more than two sequential data points. Go to the module ChemProp1 for studies with only two data points. [To read more about this](https://doi.org/10.1021/acs.analchem.1c01520)""")

# Input Files
st.subheader('Input File Requirements')
st.write(""" 
         The accepted formats for the manually uploading the input files are **CSV**, **TXT**, or **XLSX**. The necessary files include:
         1. **Feature Quantification Table**
         2. **Metadata** 
         3. **Edge File** from FBMN - This is essential for analyzing connections between node pairs (metabolites).
         4. **Annotation Files** (if available) from FBMN - These files provide additional context and annotation for the features in your dataset.
         
         **Instead of manually uploading these files, you can also provide your FBMN Job ID from GNPS2.**
         This will allow ChemProp2 to directly retrieve and process the necessary data. To get an idea of how these tables should be structured, you can use the test data available on the ‘Load Input Data’ page. \n
         **Note:** FBMN jobs from GNPS1 are not supported. 
         """)

st.write("""
### Data Preparation Essentials
To ensure smooth processing, follow these guidelines:

- Input files must include the `.mzML` extension in the feature quantification table and metadata table.
- **Metadata Table**:
  - Must include a column **filename**.
  - Can include attribute columns such as **replicates**, **sample type** (e.g., control, treatment), etc.
  - **Time**: A specific column for time is mandatory.
""")

st.markdown("""          
Example feature table:  
 
|feature|sample1.mzML|sample2.mzML|blank.mzML|
|---|---|---|---|
|1|1000|1100|100|
|2|2000|2200|200|
""")

st.write(' ')

st.markdown("""        
Example meta data table:
            
|filename|Sample_Type|Time_Point|
|---|---|---|
|sample1.mzML|Sample|1h|
|sample2.mzML|Sample|2h|
|blank.mzML|Blank|N/A| 
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

# Display the table
st.write(' ')
st.write("""
### For GNPS1 Users: Uploading Data
While ChemProp2 doesn’t support FBMN jobs from GNPS1 directly, you can manually upload the required files to the app.""")
st.write("""
         Once your FBMN job is complete, go to the **Job Status** page and press **Download Cytoscape Data**
         under **Export/Download Network Files**. 
         In the downloaded folder, you can find the required files as follows:

         """)
st.table(df)


# Output Files
st.subheader('Output File Information')
st.write("""
     
Upon processing your data, ChemProp2 generates an output in the form of a CSV file. 
This file is an enhanced version of the node-pair information, now including ChemProp2 scores for each pair. 
The scores range from -1 to +1, providing a comprehensive score for every edge (or node pair) within the molecular network.
         
Key aspects of the output CSV include:
- **Score Range**: Each node pair gets a score between -1 and +1.
- **Score Interpretation**: 
  - The magnitude of the score indicates the strength of the potential transformation. 
  - The sign of the score (positive or negative) reveals the directionality of the transformation. 
  - For example, in a node pair A & B, the sign of the score will indicate whether the transformation is from A to B or vice versa.
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