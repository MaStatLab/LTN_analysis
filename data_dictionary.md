The data used in the reproducible analysis is aggregated into a phyloseq-class object and saved as "./cache/ps_otu50.RData", which contains the following parts:

- `otu_table()`: the OTU table consisting of the counts of the OTUs
- `tax_table()`: the taxonomy table of the OTUs
- `phy_tree()`: the phylogenetic tree of the OTUs
- `sample_data()`: covariate information aggregated from the file "./res/mmc2.xlsx" and "./res/diabimmune_t1d_16s_metadata.rdata". The variables used in the reproducible analysis are described as below. 
  - `Case_Control`: factor with two levels "case" and "control", categorized by the T1D status of the subjects
  - `Age_at_Collection`: numeric, age (day) of subject at sample collection, maximum is 1233
  - `sero_age`: numeric, age (day) at seroconversion for cases; for controls, this variable is set to 1234
  - `post_seroconversion`: logical, defined as `post_seroconversion=(Age_at_Collection>sero_age)`
  - `BF`: factor, "true" if breastfeeding has not ceased, otherwise "false"
  - `Solid_Food`: factor, "true" if solid food is introduced to the subject's diet, otherwise "false"
  - `Eggs`: factor, "true" if eggs is introduced to the subject's diet, otherwise "false"
  - `Fish`: factor, "true" if fish is introduced to the subject's diet, otherwise "false"
  - `Soy_Prod`: factor, "true" if soy product is introduced to the subject's diet, otherwise "false"
  - `Rye`: factor, "true" if rye is introduced to the subject's diet, otherwise "false"
  - `Barley`: factor, "true" if barley is introduced to the subject's diet, otherwise "false"
  - `Buckwheat_Millet`: factor, "true" if buckwheat and millet is introduced to the subject's diet, otherwise "false"
  - `Country`: factor with levels "Estonia" and "Finland", country of origin of the subject
  - `Gender`: factor with levels "male" and "female", gender of the subject
  - `Subject_ID`: factor with 33 levels, subject ID 



