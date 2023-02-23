import altair as alt
import pandas as pd
import streamlit as st


@st.cache_data
def load_mutation_data():
    df = pd.read_csv("data/data_mutations.txt", sep='\t')
    return df


@st.cache_data
def load_patient_data():
    df = pd.read_csv("data/data_clinical_patient.txt", sep='\t')
    return df


@st.cache_data
def load_sample_data():
    df = pd.read_csv(
        "data/data_clinical_sample.txt",
        sep='\t',
        skiprows=4)
    return df


@st.cache_data
def load_cna_data():
    df = pd.read_csv("data/data_cna.txt", sep='\t')
    return df


sample_df = load_sample_data()
cna_df = load_cna_data()

all_cancer_types = sample_df["CANCER_TYPE"].unique()
all_cancer_genes = cna_df['Hugo_Symbol'].unique()

# Figure 3
st.write("## Genomic characterization of metastatic tumors for different distant organ sites")

# Allow the user to select cancer types of interest.
selected_cancer_types = st.multiselect(
    'Cancer type',
    options=all_cancer_types,
    default=["Non-Small Cell Lung Cancer"])

# First, we want to filter for samples associated with the cancer type.
sample_df_filtered_cancer_types = sample_df[sample_df["CANCER_TYPE"].isin(selected_cancer_types)]

# Next, we'll merge copy number alternation data for each of the filtered samples.
# 1. Transpose the CNA dataframe to turn the sample ID into indices. This will ease the merge.
cna_df_T = cna_df.set_index('Hugo_Symbol').T
cna_df_T.reset_index(inplace=True)
cna_df_T = cna_df_T.rename(columns={'index': 'SAMPLE_ID'})
# 2. Merge.
sample_cna_df_filtered_cancer_types = sample_df_filtered_cancer_types.merge(cna_df_T,
                                                                            on='SAMPLE_ID',
                                                                            # Left is filtered, right is not
                                                                            how='left')

st.write(sample_cna_df_filtered_cancer_types.head())

cancer_genes = st.multiselect(
    'Cancer gene',
    options=all_cancer_genes,
    default=["TAP1"])