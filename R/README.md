## Functions in the Package and Their Descriptions

---

### `change_gene_names`

Converts Ensembl IDs to gene symbols in a read counts file where Ensembl IDs
are in the columns, preserving any Ensembl IDs that do not have a corresponding symbol.

**Usage**
```r
data_transformed <- change_gene_names(data, transpose = TRUE, ensembl_id_col = 1, remove_cols = c(2))
```

---

### `combinatorial_gain`

Calculates the combinatorial gain of a rule, and helps generalize rules by
eliminating those with excessive objects from other classes in boundary regions.
This metric is particularly useful for analyzing rules that exhibit relatively
high accuracy but low support.

**Usage**
```r
rules_gain <- combinatorial_gain(data = autcon, rules)
```
Where `data` is a decision table with samples in rows and the outcome as the last column.

---

### `consensus_run`

Runs multiple iterations of RBM to reduce the stochasticity of the genetic reducer.

**Usage**
```r
rules <- consensus_run(data)
```
Where `data` is a decision table with samples in rows and the outcome as the last column.

---

### `findCodingMarkers`

For RNA-seq data containing pseudo and long non-coding RNAs, this function
removes them and returns a decision table with samples in rows and the outcome
as the last column, containing only protein-coding sequences.

**Usage**
```r
data_proteincoded <- findCodingMarkers(data, keytype = "SYMBOL")
```
Where `data` is a decision table with samples in rows and the outcome as the last column.

---

### `permutation_test`

Runs a permutation test to calculate the statistical significance of a rule-based model.

**Usage**
```r
result <- permutation_test(data, reducer = 'Johnson', discrete = FALSE)
```
Where `data` is a decision table with samples in rows and the outcome as the last column.

---

### `readFeatureCounts`

Converts single RNA-seq read count files (e.g., output from HTSeq) into a
decision table with samples in rows and the outcome as the last column, which
can also be used for differential expression analysis.

**Usage**
```r
path <- c("../data/data_all_old/")
df <- readFeatureCounts(path, makedt = TRUE)
```

---

### `runMCFS`

A direct pipeline that allows you to input and save files without having to
manually manage the process.

**Usage**
```r
runMCFS(data, 'mcfs_run1', permutations = 20, set_response_var = FALSE, response_var = NULL, plot_res = TRUE)
```
Where `data` is a decision table with samples in rows and the outcome as the last column.

---

### `mcfs_to_rosetta`

A pipeline that takes output from MCFS and creates a rule-based model with a
built-in stepwise selection step to select the best model.

**Usage**
```r
mcfs_to_rosetta(mcfs_df, outname = 'mcfs_run1')
```
Where `mcfs_df` is a decision table with samples in rows and the outcome as the last column.

---

### `step_wise`

A feature selection method that selects the most accurate classifier using the
top genes from the MCFS step.

**Usage**
```r
result <- step_wise(df, reducer = 'Johnson', discrete = FALSE)
```
Where `df` is a decision table with samples in rows and the outcome as the last column.
