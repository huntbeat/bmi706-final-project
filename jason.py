import altair as alt
import pandas as pd
import streamlit as st
#Disable the row limit

@st.cache_data
def load_data():
    df = pd.read_table("data/msk_met_2021_clinical_data.tsv", sep='\t')
    return df