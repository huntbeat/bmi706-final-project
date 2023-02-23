import altair as alt
import pandas as pd
import streamlit as st

### P1.2 ###

@st.cache_data
def load_data():
    # Move this code into `load_data` function {{
    cancer_df = pd.read_csv(
        "https://raw.githubusercontent.com/hms-dbmi/bmi706-2022/main/cancer_data/cancer_ICD10.csv").melt(
        # type: ignore
        id_vars=["Country", "Year", "Cancer", "Sex"],
        var_name="Age",
        value_name="Deaths",
    )

    pop_df = pd.read_csv("https://raw.githubusercontent.com/hms-dbmi/bmi706-2022/main/cancer_data/population.csv").melt(
        # type: ignore
        id_vars=["Country", "Year", "Sex"],
        var_name="Age",
        value_name="Pop",
    )

    df = pd.merge(left=cancer_df, right=pop_df, how="left")
    df["Pop"] = df.groupby(["Country", "Sex", "Age"])["Pop"].fillna(method="bfill")
    df.dropna(inplace=True)

    df = df.groupby(["Country", "Year", "Cancer", "Age", "Sex"]).sum().reset_index()
    df["Rate"] = df["Deaths"] / df["Pop"] * 100_000
    # }}
    return df

# Uncomment the next line when finished
df = load_data()

### End of P1.2 ###

st.write("## Age-specific cancer mortality rates")

### P2.1 ###
YEAR_COLUMN = "Year"
year = st.slider(
    'Select the year to filter by.',
    min_value=int(df[YEAR_COLUMN].min()),
    max_value=int(df[YEAR_COLUMN].max()))
subset = df[df[YEAR_COLUMN] == year]
### End of P2.1 ###

### P2.2 ###
# replace with st.radio
SEX_COLUMN = "Sex"
sex = st.radio(
    "Select to sex to filter by.",
    ('M', 'F'))
subset = subset[subset[SEX_COLUMN] == sex]
### End of P2.2 ###

### P2.3 ###
# replace with st.multiselect
# (hint: can use current hard-coded values below as as `default` for selector)
default_countries = [
    "Austria",
    "Germany",
    "Iceland",
    "Spain",
    "Sweden",
    "Thailand",
    "Turkey",
]
countries = st.multiselect(
    'Countries',
    options=df["Country"].unique(),
    default=default_countries)
subset = subset[subset["Country"].isin(countries)]
### End of P2.3 ###

### P2.4 ###
# replace with st.selectbox
cancer = st.selectbox(
    'Cancer',
    options=df["Cancer"].unique())
subset = subset[subset["Cancer"] == cancer]
### End of P2.4 ###

### P2.5 ###
ages = [
    "Age <5",
    "Age 5-14",
    "Age 15-24",
    "Age 25-34",
    "Age 35-44",
    "Age 45-54",
    "Age 55-64",
    "Age >64",
]

# Removing demo chart given with template
# chart = alt.Chart(subset).mark_bar().encode(
#     x=alt.X("Age", sort=ages),
#     y=alt.Y("Rate", title="Mortality rate per 100k"),
#     color="Country",
#     tooltip=["Rate"],
# ).properties(
#     title=f"{cancer} mortality rates for {'males' if sex == 'M' else 'females'} in {year}",
# )

interval_brush = alt.selection(type='interval', encodings=['x'], fields=['Age'])

pop_chart = alt.Chart().mark_bar().encode(
    x="sum(Pop)",
    y='Country:N'
).transform_filter(interval_brush).properties(
    width=600,
)

heatmap = alt.Chart().mark_rect().encode(
    x=alt.X('Age:O', sort=ages),
    y='Country:N',
    color=alt.Color('Rate:Q', scale=alt.Scale(type='log', domain=(0.01, 1000), clamp=True))
).add_selection(
  interval_brush
).properties(
    title=f"{cancer} mortality rates for {'males' if sex == 'M' else 'females'} in {year}",
    width=600,
)

### End of P2.5 ###

# st.altair_chart(chart, use_container_width=True)
all_charts = alt.vconcat(heatmap, pop_chart, data=subset)
st.altair_chart(all_charts, use_container_width=True)

countries_in_subset = subset["Country"].unique()
if len(countries_in_subset) != len(countries):
    if len(countries_in_subset) == 0:
        st.write("No data avaiable for given subset.")
    else:
        missing = set(countries) - set(countries_in_subset)
        st.write("No data available for " + ", ".join(missing) + ".")
