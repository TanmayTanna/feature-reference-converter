# Feature Reference Converter

A small, standard-library Python utility for harmonizing `target_gene_id` and
`target_gene_name` in a feature-reference CSV against a supplied annotation GTF.
Python 3.8 or newer is required. No packages, network access, or resource downloads
are required to run it.

## Usage

Use the exact GTF from the destination transcriptome reference:

```sh
python3 update_feature_csv.py destination_annotation.gtf input_features.csv converted_features.csv --strict
```

Gzip-compressed GTF files are also supported. The output must be a new path; the
script refuses to overwrite either the input or an existing output.

## Mapping Rules

1. Recognized control labels are standardized to `Non-Targeting` in both target
   annotation columns. Matching ignores punctuation and case and uses complete
   labels, not substring matches.
2. Rows labelled `Ignore` in either target annotation column retain that sentinel.
3. If a target gene ID exists in the GTF, retain the ID and use its GTF gene name.
4. Otherwise, use the existing gene name only if it maps to exactly one GTF ID.
5. Unresolved or ambiguous targets stop a `--strict` run without writing output.
   Without `--strict`, both target annotation fields become `NA`; these rows must
   be reviewed before downstream use.

All other CSV fields, row order, and column order are preserved. In particular,
the utility does not modify guide sequences, feature IDs, read identifiers,
capture patterns, or feature types. CSV quoting and line endings may differ.

Gene identifiers are matched exactly: version suffixes are not silently stripped,
and symbols are not uppercased. No orthology, genomic alignment, coordinate
liftover, alias database, or guide redesign is performed. A reused ID takes
precedence over a conflicting input symbol, so annotation changes should still
be reviewed. GTF gene records must contain both `gene_id` and `gene_name`.

## Validation

The command prints aggregate counts of unchanged, updated, control, ignored,
unresolved, and ambiguous rows. It does not print target identities.

Run the synthetic tests with:

```sh
python3 -m unittest discover -s tests -v
```

Inspect the output annotations and run the destination software's reference
validation before processing data. Successful conversion is not a substitute for
that validation. This repository contains only generic code, documentation, and
synthetic tests; it does not include reference annotations or experimental data.
