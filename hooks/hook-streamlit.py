from PyInstaller.utils.hooks import copy_metadata

datas = []

# Include metadata for all required packages
datas += copy_metadata("streamlit")
datas += copy_metadata("plotly")
datas += copy_metadata("pingouin")
datas += copy_metadata("openpyxl")
datas += copy_metadata("kaleido")
datas += copy_metadata("scikit_posthocs")
datas += copy_metadata("gnpsdata")
datas += copy_metadata("scikit_learn")
datas += copy_metadata("tabulate")
datas += copy_metadata("networkx")
datas += copy_metadata("pandas_flavor")
datas += copy_metadata("numpy")
datas += copy_metadata("scipy")  # Added scipy
datas += copy_metadata("scikit_bio")  # Added scikit-bio
datas += copy_metadata("streamlit_agraph")  # Added scipy
datas += copy_metadata("pygraphviz")

