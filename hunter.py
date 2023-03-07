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
mut_df = load_mutation_data()

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

if True:
    sample_df_filtered_cancer_types = sample_df[sample_df["CANCER_TYPE"].isin(["Breast Cancer"])]
else:
    sample_df_filtered_cancer_types = sample_df

# Next, we'll merge copy number alternation data for each of the filtered samples.
# 1. Transpose the CNA dataframe to turn the sample ID into indices. This will ease the merge.
cna_genes = set(filtered_cna_df['Hugo_Symbol'].dropna().unique())
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
mut_genes = set(filtered_mut_df['Hugo_Symbol'].dropna().unique())
filtered_mut_df = filtered_mut_df[["Hugo_Symbol", "Tumor_Sample_Barcode", "Consequence"]]
filtered_mut_df = filtered_mut_df.groupby(["Tumor_Sample_Barcode", "Hugo_Symbol"]).size().reset_index(name='Count')
filtered_mut_df = pd.pivot(filtered_mut_df, columns=['Hugo_Symbol'], index=["Tumor_Sample_Barcode"],
                           values="Count").add_prefix("mut_").reset_index()
filtered_mut_df = filtered_mut_df.rename(columns={'Tumor_Sample_Barcode': 'SAMPLE_ID'})
# 4. Merge mutation into sample table.
merged_sample_df = sample_cna_df_filtered_cancer_types.merge(filtered_mut_df,
                                                             on='SAMPLE_ID',
                                                             # Left is filtered by cancer type, right is not
                                                             how='left')
merged_sample_df = merged_sample_df.fillna(0)

st.write(merged_sample_df)


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
        alt.X('Gene:N'),
        alt.Y('Organ:N'),
    )

    # Create Amplification heatmap
    amp_heatmap = base.mark_rect().encode(
        alt.Color('Amplification_Fraction:Q', scale=alt.Scale(scheme='greenblue')),
        opacity=alt.Opacity('Amplification_Fraction:Q', scale=alt.Scale(range=[0.2, 1]), title='fraction'),
    ).properties(title='Amplification Heatmap')

    # Create Deletion heatmap
    del_heatmap = base.mark_rect().encode(
        alt.Color('Deletion_Fraction:Q', scale=alt.Scale(scheme='greenblue')),
        opacity=alt.Opacity('Deletion_Fraction:Q', scale=alt.Scale(range=[0.2, 1]), title='fraction'),
    ).properties(title='Deletion Heatmap')

    # Create Mutation heatmap
    mut_heatmap = base.mark_rect().encode(
        alt.Color('Mutation_Fraction:Q', scale=alt.Scale(scheme='greenblue')),
        opacity=alt.Opacity('Mutation_Fraction:Q', scale=alt.Scale(range=[0.2, 1]), title='fraction'),
    ).properties(title='Mutation Heatmap')

    return amp_heatmap & del_heatmap & mut_heatmap

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

valid_genes = list(cna_genes.intersection(mut_genes))
heatmap_df = build_heatmap_df(merged_sample_df, valid_genes, all_organs)
st.write(heatmap_df)
st.write(build_heatmap(heatmap_df))
st.write(build_chart(heatmap_df, 'NCOR1'))


