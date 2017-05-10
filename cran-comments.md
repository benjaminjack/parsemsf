## Test environments
* local macOS install x86_64-apple-darwin16.5.0, R 3.4.0
* R-devel and R 3.4.0 on win-builder (via `devtools::build_win()`)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:
```
Possibly mis-spelled words in DESCRIPTION:
  MSF (2:27, 5:58)
  Proteome (5:80)
  ThermoFisher (2:14, 5:45, 5:474)
  peptide (5:207)
  peptides (5:261)
```
These spellings are correct.
```
Suggests or Enhances not in mainstream repositories:
  dbplyr
```
dbplyr is a dependency for an upcoming release of dplyr. dbplyr will need to be installed when the new version of dplyr is released, but it is not required under the current version of dplyr.

## Downstream dependencies
There are no downstream dependencies.
