from seunghun import get_seunghun_charts
from hunter import *
from jason import *

import altair as alt
import pandas as pd
import streamlit as st

# Disable the row limit


st.set_page_config(layout="wide")  # Use Wide mode as default

alt.data_transformers.disable_max_rows()

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
mut_df = load_mutation_data()
cli_sample_df = pd.read_table("data/msk_met_2021_clinical_data.tsv",sep='\t')

all_cancer_types = sample_df["CANCER_TYPE"].unique()
all_cancer_genes = cna_df['Hugo_Symbol'].unique()

wanted_mutation_types = ['frameshift_variant',
                         'missense_variant',
                         'splice_acceptor_variant',
                         'splice_donor_variant',
                         'start_lost',
                         'stop_gained',
                         'stop_lost']
all_mutations = list(mut_df['Consequence'].dropna().unique())

valid_mutation_types = []
for mutation in all_mutations:
    if not str(mutation):
        continue
    for wanted_mutation in wanted_mutation_types:
        if wanted_mutation in mutation:
            valid_mutation_types.append(mutation)
            break
valid_mutation_types = all_mutations

all_organs = list(filter(lambda column_name: "DMETS_DX_" in column_name, sample_df.columns))

filtered_mut_df = mut_df[mut_df['Consequence'].isin(valid_mutation_types)]
filtered_cna_df = cna_df
cna_genes = set(filtered_cna_df['Hugo_Symbol'].dropna().unique())
mut_genes = set(filtered_mut_df['Hugo_Symbol'].dropna().unique())
valid_genes = list(cna_genes.intersection(mut_genes))

# Add sidebar features
with st.sidebar:
    figure_radio = st.radio(
        "Choose Panel",
        ("Clinical information", "Primary vs Metastasis", "Difference in organ sites")
    )
    st.markdown("#")  # Adding whitespace
    st.markdown("#")
    # Hunter: Commenting below out for now because it's causing column render issues
    # st.altair_chart(imgs, use_container_width=True)
    list_for_dropdown = all_cancer_types
    if figure_radio == "Clinical information":
        list_for_dropdown = ["Pan Cancer", *all_cancer_types]
    selected_cancer = st.selectbox(
        'Choose a cancer type',
        (list_for_dropdown))
    st.markdown("#")
    st.markdown("#")
    st.image("hms_logo.png", width=300)


@st.cache_data
def get_merged_sample_df(cancer_type, filtered_mut_df, filtered_cna_df):
    sample_df_filtered_cancer_types = sample_df[sample_df["CANCER_TYPE"].isin([cancer_type])]

    # Next, we'll merge copy number alternation data for each of the filtered samples.
    # 1. Transpose the CNA dataframe to turn the sample ID into indices. This will ease the merge.
    filtered_cna_df_T = filtered_cna_df.set_index('Hugo_Symbol').T
    filtered_cna_df_T.reset_index(inplace=True)
    filtered_cna_df_T = filtered_cna_df_T.add_prefix("cna_")
    filtered_cna_df_T = filtered_cna_df_T.rename(columns={'cna_index': 'SAMPLE_ID'})
    # 2. Merge CNA into sample table.
    sample_cna_df_filtered_cancer_types = sample_df_filtered_cancer_types.merge(filtered_cna_df_T,
                                                                                on='SAMPLE_ID',
                                                                                # Left is filtered, right is not
                                                                                how='left')
    # 3. Wrangle with mutation dataframe.
    filtered_mut_df = filtered_mut_df[["Hugo_Symbol", "Tumor_Sample_Barcode", "Consequence"]]
    filtered_mut_df = filtered_mut_df.groupby(["Tumor_Sample_Barcode", "Hugo_Symbol"]).size().reset_index(name='Count')
    filtered_mut_df = pd.pivot(filtered_mut_df, columns=['Hugo_Symbol'], index=["Tumor_Sample_Barcode"],
                               values="Count").add_prefix("mut_").reset_index()
    filtered_mut_df = filtered_mut_df.rename(columns={'Tumor_Sample_Barcode': 'SAMPLE_ID'})
    # 4. Merge mutation into sample table.
    _merged_sample_df = sample_cna_df_filtered_cancer_types.merge(filtered_mut_df,
                                                                 on='SAMPLE_ID',
                                                                 # Left is filtered by cancer type, right is not
                                                                 how='left')
    _merged_sample_df = _merged_sample_df.fillna(0)
    return _merged_sample_df

merged_sample_df = get_merged_sample_df(selected_cancer, filtered_mut_df, filtered_cna_df)

if figure_radio == "Clinical information":
    seunghun_chart = get_seunghun_charts(selected_cancer, all_cancer_types)
    st.altair_chart(seunghun_chart, use_container_width=True)
elif figure_radio == "Primary vs Metastasis":
    default_genes = valid_genes[:10]
    selected_genes = st.sidebar.multiselect("Choose genes", valid_genes, default_genes)
    jason_chart = get_jason_charts(selected_cancer, selected_genes, cna_df, mut_df, cli_sample_df)
    st.altair_chart(jason_chart, use_container_width=True)
elif figure_radio == "Difference in organ sites":
    with st.spinner("Loading... It's a lot of genes, may take a minute or two."):
        heatmap_df, selector = build_heatmap_df_and_selector(merged_sample_df, valid_genes, all_organs)
        final_chart = alt.hconcat(
            build_heatmap(heatmap_df, selector),
            build_chart(heatmap_df, selector)
        )
        st.altair_chart(final_chart, use_container_width=True)