import streamlit as st

# page_setup()

# Set page title
st.set_page_config(page_title="ChemProp2")

# Page Header
st.title("ChemProp2")

st.image("assets/ChemProp2.png", use_column_width=True)

# Under Construction Notice
st.subheader("🚧 Under Construction 🚧")

# Introduction
st.subheader('What is ChemProp2 Used For?')
st.write("""
         ChemProp2 is a tool developed to address a key challenge in non-targeted metabolomics using liquid chromatography-tandem mass spectrometry (LC–MS/MS), 
         particularly in the study of biotransformation sudies like drug metabolism as well as xenobiotic and natural product biotransformation in the environment. 
        
         Analyzing and annotating the vast data from metabolomic studies still remains a challenge. Various in silico methods and spectral similarity metrics have been developed to tackle this issue. 
         Tools like GNPS (now GNPS2) use Feature-based Molecular Networking (FBMN) to create molecular networks by connecting metabolites with similar MS/MS spectral profiles. 
         ChemProp2 builds on this, identifying potential biotransformations within these networks. It detects anti-correlating metabolites and putatibve reaction pairs, scoring their correlation over time or space. 
         This helps in prioritizing and visualizing biochemical alterations within the network.
        
         ChemProp2 is particularly useful when dealing with more than two sequential data points. Go to the module ChemProp1 for studies with only two data points. [To read more about this:](https://doi.org/10.1021/acs.analchem.1c01520)""")

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
         - **Integration with Cytoscape**: You can download the CSV file and merge it with Cytoscape for further analysis and visualization.
         
         Additionally, ChemProp2 offers a 'Smart View' feature:
         - **Visualize in Molecular Network**: Use this tab to view each node pair within a molecular network. This visualization includes the node pairs' first and/or second neighbors, along with the ChemProp2 scores and directionality indicators. You can save these visualizations as PNG images for your records or presentations.
         """)

# Citation and Resources
st.subheader('Citation and Further Resources')
st.write('If you use ChemProp2 in your research, please cite: ....')
st.write('* [FBMN-STATS](https://fbmn-statsguide.gnps2.org/) - A statistical pipeline for downstream processing of FBMN results.')
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
- Interested in contributing? Check out our [GitHub personal page](https://github.com/abzer005).
- For more about our work, visit our [lab's GitHub page](https://github.com/Functional-Metabolomics-Lab).
- Follow us on [Twitter](https://twitter.com/Functional-Metabolomics-Lab) for the latest updates.
""")

# Optional: Footer
st.markdown("---")
st.text("ChemProp2 © 2023")