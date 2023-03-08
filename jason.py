import altair as alt
import pandas as pd
import streamlit as st
#Disable the row limit

@st.cache_data

def count_amplification(df, gene, sample_type):
    return df[(df[f'cna_{gene}'] > 0) & (df['Sample_Types'] == sample_type)]['SAMPLE_ID'].nunique()

def count_deletion(df, gene, sample_type):
    return df[(df[f'cna_{gene}'] < 0) & (df['Sample_Types'] == sample_type)]['SAMPLE_ID'].nunique()

def count_mutation(df, gene, sample_type):
    return df[(df[f'mut_{gene}'] != 0) & (df['Sample_Types'] == sample_type)]['SAMPLE_ID'].nunique()

def count_samples(df, sample_type):
    return df[df['Sample_Types'] == sample_type]['SAMPLE_ID'].nunique()

def fraction(row):
    total = row['Total']
    if total == 0:
        return 0
    else:
        return int(row[row.name]) / int(total)

def get_jason_charts(selected_cancer, valid_genes):
    sample_df = pd.read_csv("data/data_clinical_sample.txt", sep="\t",skiprows=4)
    cna_df = pd.read_csv("data/data_cna.txt", sep="\t")
    mut_df = pd.read_csv("data/data_mutations.txt", sep="\t")

    sample_df = pd.read_table("data/msk_met_2021_clinical_data.tsv",sep='\t')

    sample_df['Sample_Types'] = "None"
    sample_df.loc[(sample_df['Sample Type'] == 'Primary') & (sample_df['Metastatic patient'] == True), 
                  'Sample_Types'] = "Primary_from_Met"
    sample_df.loc[(sample_df['Sample Type'] == 'Primary') & (sample_df['Metastatic patient'] == False), 
                  'Sample_Types'] = "Primary_from_NoMet"
    sample_df.loc[(sample_df['Sample Type'] == 'Metastasis'), 'Sample_Types'] = "Metastasis"

    sample_df['patient_count'] = 1
    sample_df = sample_df.rename(columns={'Sample ID': 'SAMPLE_ID'})

    filtered_mut_df = mut_df[mut_df["Hugo_Symbol"].isin(valid_genes)]
    filtered_cna_df = cna_df[cna_df["Hugo_Symbol"].isin(valid_genes)]

    filtered_cna_df_T = filtered_cna_df.set_index('Hugo_Symbol').T
    filtered_cna_df_T.reset_index(inplace=True)
    filtered_cna_df_T = filtered_cna_df_T.add_prefix("cna_")
    filtered_cna_df_T = filtered_cna_df_T.rename(columns={'cna_index': 'SAMPLE_ID'})

    sample_df_filtered_cancer_types = None
    if selected_cancer != "Pan Cancer":
        sample_df_filtered_cancer_types = sample_df[sample_df["Cancer Type"].isin([selected_cancer])]
    else:
        sample_df_filtered_cancer_types = sample_df

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
    filtered_mut_heat = filtered_mut_df[["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"]]
    filtered_mut_heat = filtered_mut_heat.groupby(["Tumor_Sample_Barcode", "Hugo_Symbol"]).size().reset_index(name='Count')
    filtered_mut_heat = pd.pivot(filtered_mut_heat, columns=['Hugo_Symbol'], index=["Tumor_Sample_Barcode"], values="Count").add_prefix("mut_").reset_index()
    filtered_mut_heat = filtered_mut_heat.rename(columns={'Tumor_Sample_Barcode': 'SAMPLE_ID'})
    # 4. Merge mutation into sample table.
    merged_sample_df = sample_cna_df_filtered_cancer_types.merge(filtered_mut_heat,
                                                                on='SAMPLE_ID',
                                                                # Left is filtered by cancer type, right is not
                                                                how='left')
    merged_sample_df = merged_sample_df.fillna(0)

    counts = []
    Pri_Met = ['Primary_from_Met', 'Primary_from_NoMet', 'Metastasis']
    primary_from_met_count = count_samples(merged_sample_df, 'Primary_from_Met')
    primary_from_no_met_count = count_samples(merged_sample_df, 'Primary_from_NoMet')
    metastasis_count = count_samples(merged_sample_df, 'Metastasis')
    for gene in valid_genes:
        for pm in Pri_Met:
            amplification_count = count_amplification(merged_sample_df, gene,pm)
            deletion_count = count_deletion(merged_sample_df, gene,pm)
            mutation_count = count_mutation(merged_sample_df, gene,pm)
            counts.append({
                'gene': gene,
                'Sample_Types': pm,
                'Amplification': amplification_count,
                'Deletion': deletion_count,
                'Mutation': mutation_count,
                'Total': primary_from_met_count,
                'Amplification_Fraction': amplification_count / primary_from_met_count,
                'Deletion_Fraction': deletion_count / primary_from_met_count,
                'Mutation_Fraction': mutation_count / primary_from_met_count
            })
            
    counts_df = pd.DataFrame(counts)

    # Base chart
    base = alt.Chart(counts_df).encode(
        alt.X('gene', sort=alt.EncodingSortField(field='gene', op='count', order='ascending')),
        alt.Y('Sample_Types', sort=alt.EncodingSortField(field='Sample_Types', op='count', order='ascending')),
    )

    # Create Amplification heatmap
    amp_heatmap = base.mark_rect().encode(
        alt.Color('Amplification_Fraction:Q', scale=alt.Scale(scheme='greenblue')),
        opacity=alt.Opacity('Amplification_Fraction:Q', scale=alt.Scale(range=[0.2, 1]), title ='fraction'),
    ).properties(width=500,title='Amplification Heatmap')

    # Create Deletion heatmap
    del_heatmap = base.mark_rect().encode(
        alt.Color('Deletion_Fraction:Q', scale=alt.Scale(scheme='greenblue')),
        opacity=alt.Opacity('Deletion_Fraction:Q', scale=alt.Scale(range=[0.2, 1]),title ='fraction'),
    ).properties(width=500,title='Deletion Heatmap')

    # Create Mutation heatmap
    mut_heatmap = base.mark_rect().encode(
        alt.Color('Mutation_Fraction:Q', scale=alt.Scale(scheme='greenblue')),
        opacity=alt.Opacity('Mutation_Fraction:Q', scale=alt.Scale(range=[0.2, 1]),title ='fraction'),
    ).properties(width=500,title='Mutation Heatmap')

    
    #cna plot
    bar_cna_df_T = filtered_cna_df.set_index('Hugo_Symbol').T
    bar_cna_df_T.reset_index(inplace=True)
    bar_cna_df_T = bar_cna_df_T.rename(columns={'index': 'SAMPLE_ID'})
    merged_bar_cna_df_T = sample_df_filtered_cancer_types.merge(bar_cna_df_T, on='SAMPLE_ID', how='left')   
    index = merged_bar_cna_df_T.columns.get_loc('patient_count')
    cna_id_vars = merged_bar_cna_df_T.columns[:index+1]
    df_tall = pd.melt(merged_bar_cna_df_T, id_vars=cna_id_vars, var_name='Gene', value_name='Copy_Number')
    df_tall['Copy_Number_Status'] = df_tall['Copy_Number'].apply(
        lambda x: 'Amplification' if x > 0 else ('Deletion' if x < 0 else None)
    )
    cna_counts = df_tall.groupby(['Gene', 'Sample_Types', 'Copy_Number_Status']).size().reset_index(name='Count')

    cna_chart = alt.Chart(cna_counts).mark_bar().encode(
        x='Gene:N',
        y='Count:Q',
        color=alt.Color('Copy_Number_Status:N', scale=alt.Scale(domain=['Deletion', 'Amplification'])),
        column= alt.Column('Sample_Types:N', header=alt.Header(title=None)) 
    ).properties(
        width=180,
        title='Copy Number Count'  
    )

    vc_mut_df = filtered_mut_df[["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"]]
    vc_mut_df = vc_mut_df.groupby(["Tumor_Sample_Barcode", "Hugo_Symbol","Variant_Classification"]).size().reset_index(name='Count')
    vc_mut_df = vc_mut_df.rename(columns={'Tumor_Sample_Barcode': 'SAMPLE_ID'})

    merged_vc_mut_df = sample_df_filtered_cancer_types.merge(vc_mut_df, on='SAMPLE_ID', how='inner')
    grouped_df = merged_vc_mut_df.groupby(['Hugo_Symbol', 'Sample_Types', 'Variant_Classification']).size().reset_index(name='count')

    vc_chart = alt.Chart(grouped_df).mark_bar().encode(
        x=alt.X('Hugo_Symbol:N', title='Gene'),
        y='count:Q',
        color='Variant_Classification:N',
        column=alt.Column('Sample_Types:N', header=alt.Header(title=None)) 
    ).properties(width=180, title='Variant Classification Count')

    combined_heatmap = amp_heatmap & del_heatmap & mut_heatmap & cna_chart & vc_chart

    return combined_heatmap
