From raw tab-separated QIIME and taxa files, the following data manipulations should be performed prior to uploading the data:

## Import files

```{r}
#| label: import-files
#| message: false

v1_v3_taxa <- read_tsv(params$v1_v3_taxa)
v1_v3_qiime <- read_tsv(params$v1_v3_qiime)
v3_v4_taxa <- read_tsv(params$v3_v4_taxa)
v3_v4_qiime <- read_tsv(params$v3_v4_qiime)

```

## Clean data

```{r}
#| label: clean-data
#| message: false

# remove non-breast cancer study-related samples
v1_v3_qiime_bc <- v1_v3_qiime %>%
  select("Taxa", matches("^(BC|MBC)"))

v3_v4_qiime_bc <- v3_v4_qiime %>%
  select("Taxa", matches("^(BC|MBC)"))
```

## Join taxa

```{r}
#| label: join-taxa
#| message: false

# Total count across samples
v1_v3_qiime_bc_taxa <- v1_v3_qiime_bc_taxa %>%
  mutate("Total Count" = rowSums(across(where(is.numeric))))

v3_v4_qiime_bc_taxa <- v3_v4_qiime_bc_taxa %>%
  mutate("Total Count" = rowSums(across(where(is.numeric))))
```