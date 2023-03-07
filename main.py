from functools import reduce

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


@st.cache_data
def load_data():
    df = pd.read_table("data/msk_met_2021_clinical_data.tsv", sep='\t')
    return df


sample_df = load_sample_data()
cna_df = load_cna_data()
mut_df = load_mutation_data()

all_cancer_types = sample_df["CANCER_TYPE"].unique()
all_cancer_genes = cna_df['Hugo_Symbol'].unique()

valid_mutation_types = ['frameshift_variant',
                        'missense_variant',
                        'splice_acceptor_variant',
                        'splice_donor_variant',
                        'start_lost',
                        'stop_gained',
                        'stop_lost']

valid_genes = ["TAP1",
               "STK19",
               "SCG5",
               "STK11",
               "CRKL",
               "MEN1",
               "B2M",
               "ERRFI1",
               "TAP2",
               "CDC73"]

clin = load_data()
# Adding Sample_Types column

clin['Sample_Types'] = "None"
clin.loc[(clin['Sample Type'] == 'Primary') & (clin['Metastatic patient'] == True), 'Sample_Types'] = "Primary_from_Met"
clin.loc[
    (clin['Sample Type'] == 'Primary') & (clin['Metastatic patient'] == False), 'Sample_Types'] = "Primary_from_NoMet"
clin.loc[(clin['Sample Type'] == 'Metastasis'), 'Sample_Types'] = "Metastasis"

# Adding patient_count column so that we can use sum(patient_count) for altair charts
clin['patient_count'] = 1

#clin = clin.head(1000)  # limiting for faster loading

# Make list of cancer types for sidebar dropdown
a0 = ['Pan Cancer']
all_cancer_list = list(set(clin['Cancer Type']))
list_for_dropdown = a0 + all_cancer_list


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
    selected_cancer = st.selectbox(
        'Choose a cancer type',
        (list_for_dropdown))
    selected_cancer_genes = st.multiselect(
        'Cancer gene',
        options=all_cancer_genes,
        default=["TAP1"],
        disabled=figure_radio in ["Pan-cancer Clinical information"])
    st.markdown("#")
    st.markdown("#")
    st.image("./hms_logo.png", width=300)

filtered_mut_df = mut_df[mut_df["Hugo_Symbol"].isin(valid_genes)]
filtered_cna_df = cna_df[cna_df["Hugo_Symbol"].isin(valid_genes)]

# Get counts of met
all_cols = (clin.columns)
met_cols = [x for x in all_cols if "Distant" in x]


def get_counts(df, cancer_type, total):
    df = df[met_cols]
    met_all = []
    counts_all = []
    for mets in met_cols:
        yes_count = sum(df[mets] == 'Yes')
        met_all.append(mets)
        counts_all.append(yes_count)
    final = pd.DataFrame(list(zip(met_all, counts_all)), columns=['organ', 'count'])
    final['primary_cancer'] = cancer_type
    final['prop'] = final['count'] / total * 100
    return final


dfs = []
for cancers in all_cancer_list:
    df_temp = clin[clin['Cancer Type'] == cancers]
    total = len(df_temp)
    dfs.append(get_counts(df_temp, cancers, total))
final = pd.concat(dfs, axis=0, ignore_index=True)

cancer_counts = clin['Cancer Type'].value_counts().rename_axis('Cancer Type').reset_index(name='counts')
primary_cancer_ordered = list(cancer_counts['Cancer Type'])

# Heatmap

##get order of organ sites total
organs = []
organ_count = []
final['organ'] = final['organ'].str.replace("Distant Mets: ", "")
for org in list(set(final.organ)):
    ocount = sum(final[final.organ == org]['count'])
    organs.append(org)
    organ_count.append(ocount)

sum_organ_count = pd.DataFrame(list(zip(organs, organ_count)), columns=['Organ', 'total'])
sum_organ_count = sum_organ_count.sort_values(by='total', ascending=False)

organ_order = list(sum_organ_count.Organ)

pan_heat = alt.Chart(final).mark_rect().encode(
    x=alt.X("primary_cancer:N", sort=primary_cancer_ordered),
    y=alt.Y("organ:N", title='', axis=alt.Axis(ticks=False, labels=False), sort=organ_order),
    color=alt.Color("prop", legend=alt.Legend(orient='none', legendX=1150, legendY=150, title="% of Met samples"))
).properties(
    width=800
)

# Total Cancer Count bar.

if selected_cancer == "Pan Cancer":
    cancers = clin['Cancer Type'].unique()

    cancer_bar = alt.Chart(cancer_counts).mark_bar(color='#17becf').encode(
        y=alt.Y('counts:Q', title='total patients'),
        x=alt.X('Cancer Type:N', sort='-y', title='', axis=alt.Axis(ticks=False, labels=False))
    ).properties(width=800, height=200)

else:
    cancer_bar = alt.Chart(cancer_counts).mark_bar(color='#17becf').encode(
        y=alt.Y('counts:Q', title='total patients'),
        x=alt.X('Cancer Type:N', sort='-y', title='', axis=alt.Axis(ticks=False, labels=False)),
        opacity=alt.condition(alt.expr.datum['Cancer Type'] == selected_cancer, alt.value(1.0), alt.value(0.2))
    ).properties(
        width=800, height=200
    )

# Metastasis count bar
final2 = final.rename(columns={'primary_cancer': "Cancer Type"})
if selected_cancer != "Pan Cancer":
    final2 = final2[final2['Cancer Type'] == selected_cancer]
met_count_bar = alt.Chart(final2).mark_bar().encode(
    y=alt.Y('organ:N', sort=organ_order, title='Organ Sites of Mets'),
    # color='Cancer Type',
    x=alt.X('sum(count):Q', title='Number of total Mets', scale=alt.Scale(reverse=True)),
    tooltip=['sum(count)', 'organ:N']
).properties(
    width=250
)

# Met age boxplot

metage = clin[clin['Age at First Mets Dx'].notnull()]  # Not all samples have metastasis

if selected_cancer == "Pan Cancer":
    met_box = alt.Chart(metage).mark_boxplot(extent='min-max', color='thistle').encode(
        x=alt.X("Cancer Type:N", title='', sort=primary_cancer_ordered, axis=alt.Axis(ticks=False, labels=False)),
        y=alt.Y("Age at First Mets Dx:Q", title="Age at First Metastasis Diagnosis")
    ).properties(
        width=800, height=200
    )

else:
    met_box = alt.Chart(metage).mark_boxplot(extent='min-max', color='thistle').encode(
        x=alt.X("Cancer Type:N", title='', sort=primary_cancer_ordered, axis=alt.Axis(ticks=False, labels=False)),
        y=alt.Y("Age at First Mets Dx:Q", title="Age at First Metastasis Diagnosis"),
        opacity=alt.condition(alt.expr.datum['Cancer Type'] == selected_cancer, alt.value(1.0), alt.value(0.2))
    ).properties(
        width=800, height=200
    )

# Sex piechart
sex_colorscale = alt.Scale(domain=['Female', 'Male'], range=['pink', 'navy'])
gender = clin[clin.Sex.notnull()]
gender_count = gender.groupby(['Cancer Type', 'Sex'], as_index=False).count()
a = gender_count[['Cancer Type', 'Sex', 'Study ID']]
a = a.rename(columns={"Study ID": "patients"})

if selected_cancer != "Pan Cancer":
    a = a[a['Cancer Type'] == selected_cancer]

donut1 = alt.Chart(a).mark_arc(innerRadius=50, outerRadius=90).encode(
    theta=alt.Theta(field="patients", aggregate="sum", type="quantitative"),
    color=alt.Color('Sex:N', scale=sex_colorscale, legend=alt.Legend(orient='none', legendX=-130))  # ,legendY=200))
    , tooltip=['Sex', 'sum(patients):Q']
).properties(
    width=250,
    height=200
)

# Sample type piechart
sampletype_count = clin.groupby(['Cancer Type', 'Sample_Types'], as_index=False).count()
sampletype_count = sampletype_count[['Cancer Type', 'Sample_Types', 'Study ID']]
sampletype_count = sampletype_count.rename(columns={"Study ID": "patients"})

type_colorscale = alt.Scale(domain=['Metastasis', 'Primary_from_Met', 'Primary_from_NoMet'],
                            range=['lightcoral', 'mediumpurple', 'cornflowerblue', ''])

if selected_cancer != "Pan Cancer":
    sampletype_count = sampletype_count[sampletype_count['Cancer Type'] == selected_cancer]
donut2 = alt.Chart(sampletype_count).mark_arc(innerRadius=50, outerRadius=90).encode(
    theta=alt.Theta(field="patients", aggregate="sum", type="quantitative"),
    color=alt.Color('Sample_Types:N', scale=type_colorscale, legend=alt.Legend(orient='none', legendX=-130))
    # ,legendY=200))
    , tooltip=['Sample_Types:N', 'sum(patients):Q']
).properties(
    width=250,
    height=200
)

# TMB

if selected_cancer != "Pan Cancer":
    clin2 = clin[clin['Cancer Type'] == selected_cancer]
else:
    clin2 = clin.copy()
tmb_box = alt.Chart(clin2).mark_boxplot(extent=1.5, clip=True, outliers=False).encode(
    x=alt.X("Sample_Types:N", title='', axis=alt.Axis(ticks=False, labels=False)),
    y=alt.Y("TMB (nonsynonymous):Q"),
    color=alt.Color('Sample_Types:N')
).properties(
    width=250, height=200
)

if selected_cancer == "Pan Cancer":
    tmb_dots = alt.Chart(clin).mark_point(size=100).encode(
        x=alt.X("Cancer Type:N", sort=primary_cancer_ordered, title='', axis=alt.Axis(ticks=False, labels=False)),
        y=alt.Y("median(TMB (nonsynonymous)):Q", title='Median TMB'),
        color=alt.Color('Sample_Types:N', scale=type_colorscale),
        tooltip=['Sample_Types', 'Cancer Type', 'median(TMB (nonsynonymous)):Q']
    ).properties(
        width=800,
        height=200
    )
else:
    tmb_dots = alt.Chart(clin).mark_point(size=100).encode(
        x=alt.X("Cancer Type:N", sort=primary_cancer_ordered, title='', axis=alt.Axis(ticks=False, labels=False)),
        y=alt.Y("median(TMB (nonsynonymous)):Q", title='Median TMB'),
        color=alt.Color('Sample_Types:N', scale=type_colorscale),
        tooltip=['Sample_Types', 'Cancer Type', 'median(TMB (nonsynonymous)):Q'],
        opacity=alt.condition(alt.expr.datum['Cancer Type'] == selected_cancer, alt.value(1.0), alt.value(0.2))
    ).properties(
        width=800,
        height=200
    )

# FGA
if selected_cancer != "Pan Cancer":
    clin2 = clin[clin['Cancer Type'] == selected_cancer]
else:
    clin2 = clin.copy()
FGA_box = alt.Chart(clin2).mark_boxplot(extent=1.5, clip=True, outliers=False).encode(
    x=alt.X("Sample_Types:N", title='', axis=alt.Axis(ticks=False, labels=False)),
    y=alt.Y("FGA:Q", title='Fraction of Genome Altered (FGA)'),
    color=alt.Color('Sample_Types:N')
).properties(
    width=250, height=200
)

if selected_cancer == "Pan Cancer":
    FGA_dots = alt.Chart(clin).mark_point(size=100).encode(
        x=alt.X("Cancer Type:N", sort=primary_cancer_ordered, title='', axis=alt.Axis(ticks=False, labels=False)),
        y=alt.Y("median(FGA):Q", title='Median FGA'),
        color=alt.Color('Sample_Types:N', scale=type_colorscale),
        tooltip=['Sample_Types', 'Cancer Type', 'median(FGA):Q']
    ).properties(
        width=800,
        height=200
    )
else:
    FGA_dots = alt.Chart(clin).mark_point(size=100).encode(
        x=alt.X("Cancer Type:N", sort=primary_cancer_ordered, title='', axis=alt.Axis(ticks=False, labels=False)),
        y=alt.Y("median(FGA):Q", title='Median FGA'),
        color=alt.Color('Sample_Types:N', scale=type_colorscale),
        tooltip=['Sample_Types', 'Cancer Type', 'median(FGA):Q'],
        opacity=alt.condition(alt.expr.datum['Cancer Type'] == selected_cancer, alt.value(1.0), alt.value(0.2))
    ).properties(
        width=800,
        height=200
    )

# MSI
if selected_cancer != "Pan Cancer":
    clin2 = clin[clin['Cancer Type'] == selected_cancer]
else:
    clin2 = clin.copy()
MSI_box = alt.Chart(clin2).mark_boxplot(extent=1.5, clip=True, outliers=False).encode(
    x=alt.X("Sample_Types:N", title='', axis=alt.Axis(ticks=False, labels=False)),
    y=alt.Y("MSI Score:Q", title='Microsatellite Instability Score (MSI))', scale=alt.Scale(zero=False)),
    color=alt.Color('Sample_Types:N')
).properties(
    width=250, height=200
)

if selected_cancer == "Pan Cancer":
    MSI_dots = alt.Chart(clin).mark_point(size=100).encode(
        x=alt.X("Cancer Type:N", sort=primary_cancer_ordered),
        y=alt.Y("median(MSI Score):Q", title='Median MSI'),
        color=alt.Color('Sample_Types:N', scale=type_colorscale),
        tooltip=['Sample_Types', 'Cancer Type', 'median(MSI Score):Q']

    ).properties(
        width=800,
        height=200
    )
else:
    MSI_dots = alt.Chart(clin).mark_point(size=100).encode(
        x=alt.X("Cancer Type:N", sort=primary_cancer_ordered),
        y=alt.Y("median(MSI Score):Q", title='Median MSI'),
        color=alt.Color('Sample_Types:N', scale=type_colorscale),
        tooltip=['Sample_Types', 'Cancer Type', 'median(MSI Score):Q'],
        opacity=alt.condition(alt.expr.datum['Cancer Type'] == selected_cancer, alt.value(1.0), alt.value(0.2))

    ).properties(
        width=800,
        height=200
    )

row3 = alt.hconcat(donut2, met_box, spacing=41)
row2 = alt.hconcat(donut1, cancer_bar, spacing=35)
row1 = alt.hconcat(met_count_bar, pan_heat, spacing=100)
row0 = alt.hconcat(tmb_box, tmb_dots)
row4 = alt.hconcat(FGA_box, FGA_dots)
row5 = alt.hconcat(MSI_box, MSI_dots)
row00 = alt.vconcat(row0, row4, row5, spacing=100)
# chart = row3 & row2 & row1 & row0 & row4
chart = row3 & row2 & row1 & row00

# col1, col2, col3 = st.beta_columns([1,60,1])
st.altair_chart(chart, use_container_width=True)
# with col1:
#    st.write("")

# with col2:

# Figure 3 - Hunter

# First, we want to filter for samples associated with the cancer type.
if selected_cancer != "Pan Cancer":
    sample_df_filtered_cancer_types = sample_df[sample_df["CANCER_TYPE"].isin([selected_cancer])]
else:
    sample_df_filtered_cancer_types = sample_df

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
filtered_mut_df = pd.pivot(filtered_mut_df, columns=['Hugo_Symbol'], index=["Tumor_Sample_Barcode"], values="Count").add_prefix("mut_").reset_index()
filtered_mut_df = filtered_mut_df.rename(columns={'Tumor_Sample_Barcode': 'SAMPLE_ID'})
# 4. Merge mutation into sample table.
merged_sample_df = sample_cna_df_filtered_cancer_types.merge(filtered_mut_df,
                                                             on='SAMPLE_ID',
                                                             # Left is filtered by cancer type, right is not
                                                             how='left')
merged_sample_df = merged_sample_df.fillna(0)

if figure_radio == "Pan-cancer Clinical information":
    st.altair_chart(chart, use_container_width=True)
elif figure_radio == "Primary vs Metastasis":
    st.write("TODO by Jason")
elif figure_radio == "Difference in organ sites":
    st.write(merged_sample_df.head())

# with col3:
#    st.write("")
