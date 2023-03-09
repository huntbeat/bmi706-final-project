import altair as alt
import pandas as pd
import streamlit as st

@st.cache_data
def build_heatmap_df_and_selector(df, genes, organs, max_num_top_genes=20):
    _heatmap_df = []
    top_genes = {}
    for gene in genes:
        total_amp = df[df[f'cna_{gene}'] > 0][f'cna_{gene}'].sum()
        total_del = df[df[f'cna_{gene}'] < 0][f'cna_{gene}'].sum()
        total_cna = total_amp + abs(total_del)
        total_mut = df[f'mut_{gene}'].sum()
        if total_cna == 0 or total_mut == 0:
            continue
        top_genes[gene] = total_mut
        for organ in organs:
            gene_data = {'Gene': gene, 'Organ': organ.removeprefix("DMETS_DX_")}
            filtered_df = df[df[organ] == 'Yes']

            _amplification = filtered_df[filtered_df[f'cna_{gene}'] > 0][f'cna_{gene}'].sum()
            gene_data['Amplification'] = _amplification
            gene_data['Amplification_Fraction'] = _amplification / float(total_cna)
            _deletion = filtered_df[filtered_df[f'cna_{gene}'] < 0][f'cna_{gene}'].sum() / float(total_cna)
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
    top_genes = list(map(lambda x: x[0], top_genes))
    _heatmap_df = pd.DataFrame(_heatmap_df)
    _heatmap_df = _heatmap_df[_heatmap_df['Gene'].isin(top_genes)]

    selector = alt.selection_single(
        fields=["Gene"],
        init={'Gene': top_genes[0]}
    )

    return _heatmap_df, selector

def build_heatmap(heatmap_df, selector):
    # Create Amplification heatmap
    amp_heatmap = alt.Chart(heatmap_df).mark_rect().encode(
        alt.X('Gene:N'),
        alt.Y('Organ:N'),
        alt.Color('Amplification_Fraction:Q', scale=alt.Scale(scheme='reds')),
    ).add_selection(selector).properties(title='Amplification Heatmap')

    # Create Deletion heatmap
    del_heatmap = alt.Chart(heatmap_df).mark_rect().encode(
        alt.X('Gene:N'),
        alt.Y('Organ:N'),
        alt.Color('Deletion_Fraction:Q', scale=alt.Scale(scheme='blues')),
    ).add_selection(selector).properties(title='Deletion Heatmap')

    # Create Mutation heatmap
    mut_heatmap = alt.Chart(heatmap_df).mark_rect().encode(
        alt.X('Gene:N'),
        alt.Y('Organ:N'),
        alt.Color('Mutation_Fraction:Q', scale=alt.Scale(scheme='greens')),
    ).add_selection(selector).properties(title='Mutation Heatmap')

    return alt.vconcat(
        amp_heatmap,
        del_heatmap,
        mut_heatmap
    ).resolve_scale(
        color='independent'
    )

def build_chart(heatmap_df, selector):
    cna_df = pd.melt(heatmap_df,
                     id_vars=['Gene', 'Organ'],
                     value_vars=['Amplification', 'Deletion'],
                     var_name='CNA type',
                     value_name='Count')
    cna_df['Count'] = cna_df['Count'].abs()

    cna_bar = alt.Chart(cna_df).mark_bar().encode(
        x='Count:Q',
        y='CNA type:N',
        color=alt.Color(
            'CNA type:N',
            scale=alt.Scale(domain=['Amplification', 'Deletion'], range=['red', 'blue'])),
        row=alt.Row('Organ:N', header=alt.Header(labelAngle=0)),
        tooltip=['Gene', 'Organ', 'Count']
    ).transform_filter(selector)

    mut_bar = alt.Chart(heatmap_df).mark_bar().encode(
        x='Mutation:Q',
        y='Organ:N',
        color=alt.value('green'),
        tooltip=['Gene', 'Organ', 'Mutation']
    ).transform_filter(selector)

    return alt.vconcat(cna_bar, mut_bar)
