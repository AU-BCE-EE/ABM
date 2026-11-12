# test

Regression and argument tests for `abm()`. Run after every change to the model code:

```
Rscript test/test_abm.R                     # test the source code (compiles src/ via pkgload)
Rscript test/test_abm.R --filter=slurry     # only tests with 'slurry' in the name
Rscript test/test_abm.R --installed         # test the installed package instead
Rscript test/test_abm.R --verbose           # also print warnings from each test
Rscript test/test_abm.R --update-reference  # accept current results as new reference
```

Each test changes one argument (or argument element) from the defaults and checks that
the run succeeds, output is complete and finite, states are non-negative, cumulative values never
decrease, ash/iNDF and slurry mass are conserved (where applicable), and that the argument has the intended effect.

Key results are compared with `reference_values.rds`. A `CHANGED` status means results differ by
more than 0.01 % from the reference. If the change is intended (e.g., a model change), re-run with
`--update-reference` and commit the new reference file.

The script exits with status 1 if any test fails or changes, so it can be used in scripts or CI.
