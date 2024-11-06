import streamlit as st
from src.common import *

from streamlit.components.v1 import html

page_setup()

# Add a tracking token
html('<script async defer data-website-id="74bc9983-13c4-4da0-89ae-b78209c13aaf" src="https://analytics.gnps2.org/umami.js"></script>', width=0, height=0)

# Page Header
#st.title("ChemProp2")

st.image("assets/ChemProp2.png", use_column_width=True)

# Introduction
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
         The accepted formats for the input files are CSV, TXT, and XLSX. The necessary files include:
         1. **Feature Quantification Table**
         2. **Metadata**
         3. **Node-Pair Information** from FBMN - This is crucial for understanding the connections between different metabolite pairs.
         4. **Annotation Files** (if available) from FBMN - These files provide additional context and annotation for the features in your dataset.
         
         Instead of manually uploading these files, you can also provide your FBMN Job ID from GNPS (or GNPS2). 
         This will allow ChemProp2 to directly retrieve and process the necessary data. To get an idea of how these tables should be structured, you can use the test data available on the ‘Load Input Data’ page.
         """)

# Output Files
st.subheader('Output File Information')
st.write("""
         Upon processing your data, ChemProp2 generates an output in the form of a CSV file. 
         This file is an enhanced version of the node-pair information, now including ChemProp2 scores for each pair. 
         The scores range from -1 to +1, providing a comprehensive score for every node pair within the molecular network.
         
         Key aspects of the output include:
         - **Score Range**: Each node pair gets a score between -1 and +1.
         - **Score Interpretation**: The magnitude of the score indicates the strength of the potential transformation. The sign of the score (positive or negative) reveals the directionality of the transformation. For example, in a node pair A & B, the sign of the score will indicate whether the transformation is from A to B or vice versa.
         - **Integration with Cytoscape**: Download results as a zip file (containing GraphML and style files) to open in Cytoscape, allowing for further analysis and detailed visualization.
         
        ### Additional Visualizations and Features
        - **Global Transformation Patterns**: Visualize overall transformation trends in your data with a score plot showing ChemProp scores against various mass differences.
        - **Filtering Options**: Use filters to refine scores by range, mass, annotation name, or cluster ID, displaying only relevant edges in a list. Selecting an edge shows:
          - The intensity trends for node pairs belonging to that edge.
          - The corresponding node pair within the molecular network as a subnetwork. (Note: Full network visualization is minimized for simplified rendering, showing only the relevant cluster and its ChemProp2 scores with directional indicators.)
        - **Export Visualizations**: Save these visualizations as PNG files for easy sharing and documentation.

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