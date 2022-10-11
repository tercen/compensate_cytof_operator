# Compensate CyTOF

##### Description

The `Compensate CyTOF operator` is an operator to perform compensation on mass
cytometry data.

##### Usage

Input projection|.
---|---
`y-axis`        | measurement values
`column`           | observations (row IDs)
`row`        | channels 

Input parameters|.
---|---
`channels`        | specify mass channels stained for & debarcode

Output relations|.
---|---
`value`        | compensated values
`compensation matrix plot`        | Computed tables include a graph of the estimated compensation matrix.

##### Details

This operator uses the compensation approach described in the [CATALYST R package](https://www.bioconductor.org/packages/devel/bioc/vignettes/CATALYST/inst/doc/preprocessing.html#compensation).

##### See Also

[normalise_cytof_operator](https://github.com/tercen/normalise_cytof_operator)
, [debarcode_cytof_operator](https://github.com/tercen/debarcode_cytof_operator)

