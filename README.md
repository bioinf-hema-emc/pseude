
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pseuDE

<!-- badges: start -->
<!-- badges: end -->

The PseuDE package is designed as a suite of functions for the
aggregation of single-cell level gene-counts and application of
differential expression (DE) analysis through the DESeq2 (Love *et al.*,
2014) package. <br>

A paper regarding the abundance of false positives in single-cell DE
analyses (Squair *et al.*, 2021), clearly states that the DE algorithms
applied to single-cell level gene counts are severely lacking in
statistical strength and biased toward variable genes. Several times
over they indicate that a “pseudobulk” approach to DE analysis is proven
to be the better option. This package was made following the
recommendations given in that paper.

## Installation

You can install the development version of pseuDE from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("weversMJW/pseude")
```

Then you can load the functions into memory by running:

``` r
library(pseude)
```

PseuDE makes use of two packages that are only distributed by
BioConductor: DESeq2 and BiocParallel. If you do not have these packages
installed locally, you can run the following code:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BiocParallel")
BiocManager::install("DESeq2")
```

## Examples

The three main functions from pseuDE are `scagg`, `sciff` and `scout`.
These functions aggregate the gene-counts, perform DE analysis and give
the output file, respectively. <br>

------------------------------------------------------------------------

#### Function 1: scagg

<br>

##### Input

1.  `sc_counts`: A cell-level gene-count matrix of form *m* rows x *n*
    columns, with *m* being the genes and *n* the cells. If you’re
    working with Seurat (Hao *et al.*, 2021), such an object can be
    found by calling the following:

``` r
seurat_counts <- seurat_object[["RNA"]]@counts
```

Which would look like this sparse matrix:

|               | AACAAGCGTCAT | GCAAGTAACACT | GAATCTCTGTCC |  …  |
|---------------|:------------:|:------------:|:------------:|:---:|
| RP5-857K21.4  |      4       |      .       |      .       |  …  |
| FO538757.2    |      6       |      3       |      5       |  …  |
| RP4-669L17.10 |      1       |      .       |      .       |  …  |
| RP11-34P13.7  |      .       |      .       |      8       |  …  |
| AP006222.2    |      .       |      .       |      9       |  …  |
| …             |      …       |      …       |      …       |  …  |

2.  `sc_table`: A look-up table that contains the groups you wish to
    aggregate your samples to. Take this random sampling as an example:
    <br>

|              | sample | state |
|--------------|:------:|:-----:|
| AACAAGCGTCAT |   y4   | young |
| GCAAGTAACACT |   o7   |  old  |
| GAATCTCTGTCC |   o7   |  old  |
| GTGCCTTATTCC |   y3   | young |
| GCTTGCCGGGTG |   y6   | young |
| …            |   …    |   …   |

3.  `sc_group`: A character string that indicate which column(s) in
    `sc_table` contains the groups you wish to use for aggregation.

<br>

##### Running the code

The most rudimentary aggregation is from cell to sample. If we take the
examples from the previous section, the call would look like this: <br>

``` r
aggregated_counts <- scagg(sc_counts = cell_counts, sc_table = cell_table,
                           sc_group = "sample")
```

Say you have clustered the cells as well and you want to aggregate to
clusters as you consider them the different cell types. Then, something
like the following would do:

``` r
aggregated_counts <- scagg(sc_counts = cell_counts, sc_table = cell_table,
                           sc_group = "clusters")
```

You can also use multiple groupings. For example, if you want
sample-specific clusters:

``` r
aggregated_counts <- scagg(sc_counts = cell_counts, sc_table = cell_table,
                           sc_group = c("sample", "clusters"))
```

<br>

##### Output

the `aggregated_counts` object will contain a *m* x *n’* sparse matrix,
with *n’* being the cells aggregated to whichever level you set it to.

------------------------------------------------------------------------

#### Function 2: sciff

<br>

##### Input

1.  `agg_counts`: The *aggregated_counts* object.
2.  `agg_table`: A look-up table, similar to the table used in scagg,
    but now on the aggregated level; linking the aggregated samples to
    the variable(s) you wish to compare. For example, you’ve aggregated
    your cell-level counts to sample-level counts:

| sample | state |
|:------:|:-----:|
|   y3   | young |
|   y4   | young |
|   y6   | young |
|   o6   |  old  |
|   o7   |  old  |
|   …    |   …   |

3.  `design`: A design that describes how the counts for each gene
    depend on the variable(s) you wish to compare. For example, if you
    are only interested in the *state* of your samples, you set design:

<!-- -->

    design <- "state"`

Interested in comparing two cell types in a patient-specific paired
test? Simply run:

    design <- "type + patient"

4.  `comparison`: This is basically a character string that indicates
    which column in the look-up table it needs to use and which
    variables you which to compare against each other. For example:

<!-- -->

    comparison <- "state,young,old"

Be mindful that this follows young VS old, meaning that everything with
a **positive Fold Change (FC)** is up-regulated in **young** and
everything with a **negative FC** is up-regulated in **old**. <br> 5.
`do_prefilter` (optional): Defaults to TRUE, this removes genes that
only have a counts in a few samples. <br> 6. `threads` (optional):
Defaults to 10 and sets the number of BiocParallel (Morgan *et al.*,
2022) workers for the DE analysis.

**Always** keep in mind that the character strings you assign to the
design and comparison variables are case-sensitive.

<br>

##### Running the code

The most basic function call looks like this:

``` r
aggregated_dds <- sciff(agg_counts = aggregated_counts,
                          agg_table = sample_table,
                          design = "state",
                          comparison = "state,young,old")
```

If you do not want to prefilter the genes and speed up the process a
little, you can alter the function call like so:

``` r
aggregated_dds <- sciff(agg_counts = aggregated_counts,
                          agg_table = sample_table,
                          design = "state",
                          comparison = "state,young,old",
                          do_prefilter = FALSE,
                          threads = 12)
```

<br>

##### Output

The `aggregated_dds` object consists out of a DESeqDataSet object with
the results of the DE analysis attached as metadata columns.

------------------------------------------------------------------------

#### Function 3: scout

<br>

##### Input

1.  `dds`: A DESeqDataSet object with attached results, like the one
    generated by `sciff`.
2.  `comparison`: See Input point 4 for the `sciff` function for more
    details.
3.  `out_path` (optional): Defaults to work directory. Indicates to
    which directory the output file needs to be written.
4.  `do_write` (optional): This argument defaults to TRUE, which means
    that it will generate an .xlsx file with the DESeq2 DE genes.
5.  `do_return` (optional): This argument default to FALSE. If you wish
    to assign the output within your rsession, set this to TRUE.

<br>

##### Running the code

Running the following call will generate an output xlsx file:

``` r
aggregated_output <- scout(dds = aggregated_dds, comparison = "state,young,old")
```

If you do not wish to generate an output file:

``` r
aggregated_output <- scout(dds = aggregated_dds, comparison = "state,young,old",
                           do_write = FALSE, do_return = TRUE)
```

<br>

##### Output

`aggregated_output` will contain a table with the genes ordered by
adjusted p-value, the corresponding statistics and the counts for your
aggregated samples/clusters/objects. <br> Likewise, the xlsx file that’s
been generated will contain three sheets, with the last corresponding to
the table found in the `aggregated_output` object. The first sheet will
only contain the statistics of those genes with a adjusted p-value \<
0.05 and the second sheet will have the counts for those genes attached.

------------------------------------------------------------------------

#### Additional function: scuality

The function `scuality` generates a document with basic plots for visual
inspection. Right now it will generate PCA plots for the regularized log
and variance stabilized transform, a MAplot and a volcano plot of the DE
results. It’s also possible to assign these results as a global variable
without writing to a file.

<br>

##### Input

1.  `dds`: A DESeqDataSet object with attached results, like the one
    generated by `sciff`.
2.  `comparison`: See Input point 4 for the `sciff` function for more
    details.
3.  `out_path` (optional): Defaults to work directory. Indicates to
    which directory the output file needs to be written.
4.  `do_write` (optional): This argument defaults to TRUE, which means
    that it will generate a file called quality_metrics.pdf with the
    various plots.
5.  `do_labels` (optional): This argument defaults to TRUE, which means
    that it will give the points in the PCA plots the corresponding
    sample lables.
6.  `do_rlog` (optional): This argument defaults to TRUE, which means
    that it will do both regularized log and variance stabilized
    transfrom on the DE results. You can turn of the rlog, as it might
    take quite long with larger data sets.
7.  `do_return` (optional): This argument defaults to FALSE, which means
    that it will not assign any results to global variables.

<br>

##### Running the code

Running the following call will generate the quality_metrics file; it is
the most basic call to this function:

``` r
scuality(dds = aggregated_dds, comparison = "state,young,old",
      out_path = "example/output/path")
```

If you do not wish to generate an output file but rather return the plot
values and wish to speed the process up:

``` r
aggregated_show <- scuality(dds = aggregated_dds, comparison = "state,young,old",
                           do_write = FALSE, do_rlog = FALSE, do_return = TRUE)
```

<br>

##### Output

A file called quality_metrics.pdf will be placed wherever out_path is
pointing to. `aggregated_show` will contain a list with the plots
related to the quality of the genes.

------------------------------------------------------------------------

#### Additional function: scupple

The function `scupple` is a wrapper around `scagg`, `sciff`, `scout` and
`scuality`. It runs through the default steps and simplifies the process
for the user. However, the current form of `scupple` has less parameters
to fine-tune. The following parameters are changed: 1. The parameter
writeorret is a binary choice that either writes the various files
**or** returns the output to an assigned variable. 2. The parameter
`do_rlog` is hardcoded to TRUE. Regularised log is the recommended
option anyway, but that means it will take a bit longer to run through
the steps. 3. You cannot turn off `do_labels` for visualisation.

<br>

##### Input

1.  `sc_counts`: See Input point 1 from the `scagg` function for more
    details.
2.  `sc_table`: See Input point 2 from the `scagg` function for more
    details.
3.  `sc_group`: See Input point 3 from the `scagg` function for more
    details.
4.  `design`: See Input point 3 from the `sciff` function for more
    details.
5.  `comparison`: See Input point 4 from the `sciff` function for more
    details.
6.  `do_prefilter`: See Input point 5 from the `sciff` function for more
    details.
7.  `threads`: See Input point 6 from the `sciff` function for more
    details.
8.  `out_path`: Defaults to work directory. Indicates to which directory
    the output file needs to be written.
9.  `writeorret`: Defaults to TRUE, writing the files to
    quality_metrics.pdf and \[file_name\].DESeq2.xlsx. If set to FALSE,
    will return the data instead.

<br>

##### Running the code

Running the following call will generate the various files in dir
`/example/output/path` quality_metrics file; it is the most basic call
to this function:

``` r
scupple(sc_counts = raw_counts, sc_table = raw_table, sc_group = "sample_name",
        design = "state", comparison = "state,young,old", out_path = "/example/output/path",
        threads = 10, do_prefilter = TRUE, writeorret = TRUE)
```

If you do not wish to generate an output file but rather return the plot
values and wish to speed the process up:

``` r
aggregated_scupple <- scupple(sc_counts = raw_counts, sc_table = raw_table, sc_group = "sample_name",
                              design = "state", comparison = "state,young,old", threads = 10, do_prefilter = TRUE,
                              writeorret = FALSE)
```

<br>

##### Output

Files called quality_metrics.pdf and \[file_name\].DESeq2.xlsx will be
placed wherever out_path is pointing to. `aggregated_scupple` will
contain a list with the plots related to the quality of the genes.

## To Do

-   This is a development version of the pseuDE package, functions are
    subject to change.

## Bibliography

1.  Love MI, Huber W, Anders S (2014). “Moderated estimation of fold
    change and dispersion for RNA-seq data with DESeq2.” Genome Biology,
    15, 550. doi: 10.1186/s13059-014-0550-8.
2.  Squair, J.W., Gautier, M., Kathe, C. et al. Confronting false
    discoveries in single-cell differential expression. Nat Commun 12,
    5692 (2021). <https://doi.org/10.1038/s41467-021-25960-2>
3.  Hao Y, Hao S, Andersen-Nissen E, III WMM, Zheng S, Butler A, Lee MJ,
    Wilk AJ, Darby C, Zagar M, Hoffman P, Stoeckius M, Papalexi E,
    Mimitou EP, Jain J, Srivastava A, Stuart T, Fleming LB, Yeung B,
    Rogers AJ, McElrath JM, Blish CA, Gottardo R, Smibert P, Satija R
    (2021). “Integrated analysis of multimodal single-cell data.” Cell.
    <doi:10.1016/j.cell.2021.04.048>,
    <https://doi.org/10.1016/j.cell.2021.04.048>.
4.  Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022).
    BiocParallel: Bioconductor facilities for parallel evaluation. R
    package version 1.30.3,
    <https://github.com/Bioconductor/BiocParallel>.
