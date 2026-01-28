import streamlit as st
import pandas as pd
import io
import uuid
import base64

def clear_cache_button():
   if st.button("Clear Cache"):
        # Clear cache for both newer and older Streamlit versions
        if hasattr(st, "cache_data"):
            st.cache_data.clear()
        if hasattr(st, "cache_resource"):
            st.cache_resource.clear()
        st.success("Cache cleared!")

    # initialize global session state variables if not already present
    # DataFrames

def v_space(n, col=None):
    for _ in range(n):
        if col:
            col.write("")
        else:
            st.write("")

def page_setup():
    # streamlit configs
    st.set_page_config(
        page_title="ChemProp2",
        page_icon="assets/chemprop2.ico", 
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )

    # initialize global session state variables if not already present
    # DataFrames
    for key in dataframe_names:
        if key not in st.session_state:
            st.session_state[key] = pd.DataFrame()
    if "data_preparation_done" not in st.session_state:
        st.session_state["data_preparation_done"] = False

    with st.sidebar:
        with st.expander("⚙️ Settings", expanded=True):
            st.selectbox(
                "image export format",
                ["svg", "png", "jpeg", "webp"],
                key="image_format",
            )

            v_space(1)
            # Add the clear cache button
            clear_cache_button()
        
        v_space(1)
        
        col1, col2 = st.columns([1.5, 1])
        
        with col1:
            st.markdown(
                f'''
                <a href="https://www.functional-metabolomics.com/resources" target="_blank">
                    <img src="data:image/png;base64,{base64.b64encode(open("assets/chemprop2_logo.png", "rb").read()).decode()}"
                       style="max-width:100%; height:auto;">
                </a>
                ''',
                unsafe_allow_html=True,
                )
        with col2:
             st.markdown(
                f'''
                <a href="https://gnps2.org/homepage" target="_blank">
                    <img src="data:image/png;base64,{base64.b64encode(open("assets/GNPS2_logo.png", "rb").read()).decode()}"
                       style="max-width:100%; height:auto;">
                </a>
                ''',
                unsafe_allow_html=True,
            )
        v_space(1)

        st.markdown("## Functional-Metabolomics-Lab")
        c1, c2, c3 = st.columns(3)
        c1.markdown(
            """<a href="https://github.com/Functional-Metabolomics-Lab/ChemProp-Web-App">
            <img src="data:image/png;base64,{}" width="50">
            </a>""".format(
                base64.b64encode(open("./assets/github-mark.png", "rb").read()).decode()
            ),
            unsafe_allow_html=True,
        )
        c2.markdown(
            """<a href="https://www.functional-metabolomics.com/">
            <img src="data:image/png;base64,{}" width="50">
            </a>""".format(
                base64.b64encode(open("./assets/fmlab_logo_colored.png", "rb").read()).decode()
            ),
            unsafe_allow_html=True,
        )
        c3.markdown(
            """<a href="https://www.youtube.com/@functionalmetabolomics">
            <img src="data:image/png;base64,{}" width="50">
            </a>""".format(
                base64.b64encode(open("./assets/youtube-logo.png", "rb").read()).decode()
            ),
            unsafe_allow_html=True,
        )
        

dataframe_names = ("md",
                   "ft",
                   "nw",
                   "an_gnps",
                   "an_analog")


def show_input_tables_in_tabs():
    """
    Display Metadata, Metabolomics FT, and Other Omics FT
    in three tabs with shape info and full table view.
    """
    tab_defs = [
        ("ft", "Feature Quantification Table"),
        ("md", "Metadata"),
        ("nw", "Edge Table"),
        ("an_gnps", "Annotation Table from GNPS"),
    ]

    tabs = st.tabs([label for _, label in tab_defs])

    for tab, (key, label) in zip(tabs, tab_defs):
        with tab:
            df = st.session_state.get(key)

            if isinstance(df, pd.DataFrame) and not df.empty:
                # Show shape info
                num_rows, num_cols = df.shape
                st.caption(f"{num_rows} rows × {num_cols} columns")

                # Show full dataframe
                st.dataframe(df, use_container_width=True)
            else:
                st.info(f"{label} not loaded yet.")


def display_tables_in_tabs(tab_defs):
    """
    Display tables with their shape info and full table view.
    """
    #tab_defs = [
    #    ("ft", "Feature Quantification Table"),
    #    ("md", "Metadata"),
    #]

    tabs = st.tabs([label for _, label in tab_defs])

    for tab, (key, label) in zip(tabs, tab_defs):
        with tab:
            df = st.session_state.get(key)

            if isinstance(df, pd.DataFrame) and not df.empty:
                # Show shape info
                num_rows, num_cols = df.shape
                st.caption(f"{num_rows} rows × {num_cols} columns")

                # Show full dataframe
                st.dataframe(df, use_container_width=True)
            else:
                st.info(f"{label} not loaded yet.")

def reset_dataframes():
    for key in dataframe_names:
        st.session_state[key] = pd.DataFrame()

def open_df(file):
    separators = {"txt": "\t", "tsv": "\t", "csv": ","}
    try:
        if type(file) == str:
            ext = file.split(".")[-1]
            if ext != "xlsx":
                df = pd.read_csv(file, sep=separators[ext])
            else:
                df = pd.read_excel(file)
        else:
            ext = file.name.split(".")[-1]
            if ext != "xlsx":
                df = pd.read_csv(file, sep=separators[ext])
            else:
                df = pd.read_excel(file)
        
        # sometimes dataframes get saved with unnamed index, that needs to be removed
        if "Unnamed: 0" in df.columns:
            df.drop("Unnamed: 0", inplace=True, axis=1)
        return df
    except:
        return pd.DataFrame()

def show_table(df, title="", col="", download=True):
    if col:
        col = col
    else:
        col = st
    if download:
        col.download_button(
            f"Download Table",
            df.to_csv(sep="\t").encode("utf-8"),
            title.replace(" ", "-") + ".tsv",
            key=uuid.uuid1(),
        )
    col.dataframe(df, use_container_width=True)


def show_fig(fig, download_name, container_width=True):

    # Set default image format to 'svg' if not specified in session state
    image_format = st.session_state.get('image_format', 'svg')

    st.plotly_chart(
        fig,
        use_container_width=container_width,
        config={
            "displaylogo": False,
            "modeBarButtonsToRemove": [
                "zoom",
                "pan",
                "select",
                "lasso",
                "zoomin",
                "autoscale",
                "zoomout",
                "resetscale",
            ],
            "toImageButtonOptions": {
                "filename": download_name,
                "format": image_format,
            },
        },
    )


def download_plotly_figure(fig, filename="", col=""):
    buffer = io.BytesIO()
    fig.write_image(file=buffer, format="png")

    if col:
        col.download_button(
            label=f"Download Figure",
            data=buffer,
            file_name=filename,
            mime="application/png",
        )
    else:
        st.download_button(
            label=f"Download Figure",
            data=buffer,
            file_name=filename,
            mime="application/png",
        )
