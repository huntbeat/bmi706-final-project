import altair as alt
import pandas as pd
import streamlit as st


def build_heatmap_df(df, genes, organs, max_num_top_genes=20):
    _heatmap_df = []
    top_genes = {}
    for gene in genes:
        total_cna = df[f'cna_{gene}'].sum()
        total_mut = df[f'mut_{gene}'].sum()
        if total_cna == 0 or total_mut == 0:
            continue
        top_genes[gene] = total_mut
        for organ in organs:
            gene_data = {'Gene': gene, 'Organ': organ}
            filtered_df = df[df[organ] == 'Yes']

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

    if len(top_genes) > max_num_top_genes:
        top_genes = sorted(top_genes.items(), key=lambda x: x[1], reverse=True)[:max_num_top_genes]
    top_genes = map(lambda x: x[0], top_genes)
    _heatmap_df = pd.DataFrame(_heatmap_df)
    _heatmap_df = _heatmap_df[_heatmap_df['Gene'].isin(top_genes)]
    return _heatmap_df

def build_heatmap(heatmap_df):
    # Base chart
    base = alt.Chart(heatmap_df).encode(

    )

    # Create Amplification heatmap
    amp_heatmap = alt.Chart(heatmap_df).mark_rect().encode(
        alt.X('Gene:N'),
        alt.Y('Organ:N'),
        alt.Color('Amplification_Fraction:Q', scale=alt.Scale(scheme='greenblue')),
        opacity=alt.Opacity('Amplification_Fraction:Q', scale=alt.Scale(range=[0.2, 1]), title='fraction'),
    ).properties(title='Amplification Heatmap')

    # Create Deletion heatmap
    del_heatmap = alt.Chart(heatmap_df).mark_rect().encode(
        alt.X('Gene:N'),
        alt.Y('Organ:N', axis=alt.Axis(labels=False)),
        alt.Color('Deletion_Fraction:Q', scale=alt.Scale(scheme='greenblue')),
        opacity=alt.Opacity('Deletion_Fraction:Q', scale=alt.Scale(range=[0.2, 1]), title='fraction'),
    ).properties(title='Deletion Heatmap')

    # Create Mutation heatmap
    mut_heatmap = alt.Chart(heatmap_df).mark_rect().encode(
        alt.X('Gene:N'),
        alt.Y('Organ:N', axis=alt.Axis(labels=False)),
        alt.Color('Mutation_Fraction:Q', scale=alt.Scale(scheme='greenblue')),
        opacity=alt.Opacity('Mutation_Fraction:Q', scale=alt.Scale(range=[0.2, 1]), title='fraction'),
    ).properties(title='Mutation Heatmap')

    final_heatmap = alt.hconcat(
        amp_heatmap,
        del_heatmap,
        mut_heatmap
    ).resolve_scale(
        y='shared'
    )

    return final_heatmap

def build_chart(heatmap_df, selected_gene):
    filtered_df = heatmap_df[heatmap_df['Gene'] == selected_gene]
    cna_df = pd.melt(filtered_df,
                     id_vars=['Gene', 'Organ'],
                     value_vars=['Amplification', 'Deletion'],
                     var_name='CNA type',
                     value_name='Count')

    cna_bar = alt.Chart(cna_df).mark_bar().encode(
        x='CNA type:N',
        y='Count:Q',
        color='CNA type:N',
        column='Organ:N'
    )

    mut_bar = alt.Chart(filtered_df).mark_bar().encode(
        x='Organ:N',
        y='Mutation:Q'
    )
    return cna_bar & mut_bar

