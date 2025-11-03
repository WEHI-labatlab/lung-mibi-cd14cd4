# TNF⍺-mediated myeloid-instructed CD14+CD4+ T cells in NSCLC tumor microenvironment

Code and analysis supporting the paper: "TNF⍺-mediated myeloid-instructed CD14+CD4+ T cells within the tumor microenvironment are associated with poor survival in non-small cell lung cancer"

## Summary of the paper

In this work, we profiled a clinically annotated cohort of 64 non-small cell lung cancer (NSCLC) patient samples using a spatial proteomics technology called multiplex ion beam imaging (MIBIscope). We designed a bespoke antibody panel of 38 proteins to interrogate discrete immune cell populations and their spatial organisation in the TME and their association with clinical outcomes.
We demonstrated the existence of a subset of CD4+ T cells expressing the myeloid marker CD14. Tumors with high infiltration of CD14+CD4+ T cells were associated with reduced patient survival. Importantly we integrated spatial proteomics analysis with spatial transcriptomic approaches (GeoMx) and discovered the activation of TNF⍺ signaling in tumors with high infiltration of CD14+CD4+ T cells in the tumor core.
In vitro, TNF⍺ increased the transfer of membrane from myeloid cells to T cells through trogocytosis.

## Other repositories used in this project

You can find details of cell classification pipelines that was adapted to Nexflow here:

* [Data preprocessing](https://github.com/BioimageAnalysisCoreWEHI/MIBI-preprocess-data)
* [Training of XGBoost model](https://github.com/BioimageAnalysisCoreWEHI/MIBI-train-model)
* [Applying the model](https://github.com/BioimageAnalysisCoreWEHI/MIBI-apply-model)

## Data

Associated with this repo data is available on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17051520.svg)](https://doi.org/10.5281/zenodo.17051520).

The raw image data from MIBI is available by request.
