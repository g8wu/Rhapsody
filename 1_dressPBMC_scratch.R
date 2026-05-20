# Add patient metadata
rds$condition <- case_when(
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_hs"~ "act",
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_hs"~ "trm",
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_hs"~ "rec",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag08_hs"~ "act",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag09_hs"~ "rec",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_hs"~ "act",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_hs"~ "trm",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_hs"~ "rec",
  rds$orig.ident == "DRS03-act-BF" ~ "act-BF",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

rds$patient <- case_when(
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_hs"~ "DRS01",
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_hs"~ "DRS01",
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_hs"~ "DRS01",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag08_hs"~ "DRS02",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag09_hs"~ "DRS02",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_hs"~ "DRS03",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_hs"~ "DRS03",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_hs"~ "DRS03",
  rds$orig.ident == "DRS03-act-BF" ~ "DRS03",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

table(rds$patient)
table(rds$condition)

rds$pat.cond <- paste(rds$patient, rds$condition)
Idents(rds) <- rds$pat.cond
Idents(rds) <- factor(Idents(rds), levels = sort(levels(rds)))
rds$pat.cond <- Idents(rds)
write.csv(table(rds$patient, rds$condition),  paste0(rds@project.name, "-pat.cond.csv"))
table(rds$pat.cond)
