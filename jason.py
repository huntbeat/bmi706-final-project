import altair as alt
import pandas as pd
import streamlit as st
#Disable the row limit

@st.cache_data
def load_data():
    df = pd.read_table("data/msk_met_2021_clinical_data.tsv", sep='\t')
    return df

def build_heatmap_types(df, genes, sample_types, max_num_top_genes=20):
    _heatmap_df = []
    top_genes = {}
    for gene in genes:
        total_cna = df[f'cna_{gene}'].sum()
        total_mut = df[f'mut_{gene}'].sum()
        if total_cna == 0 or total_mut == 0:
            continue
        top_genes[gene] = total_mut
        for sample_type in sample_types:
            gene_data = {'Gene': gene, 'Organ': organ}
            filtered_df = df[df[sample_type] == 'Yes']

            _amplification = filtered_df[df[f'cna_{gene}'] > 0][f'cna_{gene}'].sum()
            gene_data['Amplification'] = _amplification
            gene_data['Amplification_Fraction'] = _amplification / float(total_cna)
            _deletion = filtered_df[df[f'cna_{gene}'] < 0][f'cna_{gene}'].sum() / float(total_cna)
            gene_data['Deletion'] = _deletion
            gene_data['Deletion_Fraction'] = abs(_deletion) / float(total_cna)
            _mutation = filtered_df[f'mut_{gene}'].sum()
            gene_data['Mutation'] = _mutation
            gene_data['Mutation_Fraction'] = _mutation / float(total_mut)

            gene_data['Total MUT'] = total_mut
            gene_data['Total CNA'] = total_cna

            _heatmap_df.append(gene_data)