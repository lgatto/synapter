# Synapter 2.21

## Changes in version 2.21.1

- Render fragmentmatching.Rmd with rmarkdown::html_document() (see #139)

# Synapter 2.17

## Changes in version 2.17.3
- put commented code chunk back in

## Changes in version 2.17.2
- fix error in vignette

## Changes in version 2.17.1

- Update Laurent's email address

# Synapter 2.5

## Changes in version 2.5.2
- Use `BiocManager::install` [2018-07-16].

## Changes in version 2.5.0
- New version for Bioconductor devel

# Synapter 2.4

## Changes in version 2.4.0
- New version for Bioconductor release

# Synapter 2.3

## Changes in version 2.3.1

- Partly revert 95f4094 because MSnbase:::utils.applyColumnwiseByGroup is gone.
- Use `html_document` instead of `html_document2` in the vignettes [2018-01-17].

# Synapter 2.1

## Changes in version 2.1.1
- Partly revert 95f4094 because MSnbase:::utils.applyColumnwiseByGroup is gone.

# Synapter 2.0
- Version bump for Bioc release version 3.5.

# Synapter 1.99

## Changes in version 1.99.2
- update show,MasterPeptides to check of there's a fragmentlibrary
  slot before trying to access it <2017-04-19 Wed>

## Changes in version 1.99.1

- update NEWS file <2017-04-11 Tue>

## Changes in version 1.99.0

### New features

#### Ion mobility/Grid search.

- Replace 2D grid search (retention time, *m/z*) of synapter1 by 3D grid search
   (retention time, *m/z*, ion mobility); set argument `imdiff = Inf` to get the
   original 2D grid search; closes [#33](https://github.com/lgatto/synapter/issues/33).
- Add `{set,get}ImDiff` methods.
- `getGrid` returns an array instead of a matrix (because of the new 3D grid
    search) [2014-05-16 Fri].
- `plotFeatures(..., what = "all")` gains a new argument: "ionmobilty" to plot
   *m/z* vs ionmobility as well. [2014-05-16 Fri]
- `plotGrid` gets a new argument "maindim" to decided which of the three
   dimension should be used. [2014-05-16 Fri]
- Add `filterNonUniqueIdentMatches` to remove matches of multiple
   identification data to a single quantitation entry (see [#111](https://github.com/lgatto/synapter/issues/111) for details)
   [2016-02-22 Mon].

#### Fragment matching

- Load identification fragments (`final_fragments.csv`) and
   quantitation spectra (`Spectrum.xml`) via `Synapter` constructor.
- New functions: `fragmentMatchingPerformance`, `filterUniqueMatches`,
    `filterNonUniqueMatches`, `filterFragments`,
    `plotCumulativeNumberOfFragments`, `plotFragmentMatchingPerformance`,
    `getIdentificationFragments` and `getQuantitationSpectra`.
- Integrate a fragment library into *master* objects; closes [#63](https://github.com/lgatto/synapter/issues/63) and [#74](https://github.com/lgatto/synapter/issues/74).

#### Misc

- Allow to use an RDS instead of a fasta file as 'Unique Peptides Database',
   adds `createUniquePeptideDbRds`; closes [#55](https://github.com/lgatto/synapter/issues/55) [2014-04-29 Tue].
- Introduce `IisL` argument to `dbUniquePeptideSet` which treats I/L as same
   aminoacid if `IisL == TRUE` (default: `IisL = FALSE`);
   closes [#60](https://github.com/lgatto/synapter/issues/60) [2014-04-30 Wed].
- Add `rescueEMRTs` functions; replaces the argument `mergedEMRTs` in
    `findEMRTs`; closes [#93](https://github.com/lgatto/synapter/issues/93) [2015-07-26 Sun].
- Add `synergise2` which combines the integrates the new 3D grid search, the
    fragment matching; and uses slightly different default arguments than
    `synergise1`; closes [#119](https://github.com/lgatto/synapter/issues/119) [2016-10-25 Di].
- Load isotopic distributions from Pep3D data and also export them to MSnSet, to
    allow the correction of detector saturation; closes [#39](https://github.com/lgatto/synapter/issues/39) [2015-03-29 Sun].
- Add `synapterPlgsAgreement` to find agreement between *synapter* and *PLGS*;
    closes [#73](https://github.com/lgatto/synapter/issues/73).
- Introduce `modelIntensity` to correct systematic intensity shifts
   (similar to `modelRt`); closes [#116](https://github.com/lgatto/synapter/issues/116).

### Improvements

- Extract the ion that was used for identification (`isFid == 1`) from the Pep3D
   file instead of the first instance [2014-05-13 Tue].
- Add `updateObject` and `validObject` method [2014-11-16 Sun].
- Rename `QuantPep3DData$Function` column into `QuantPep3DData$matchedEMRTs`;
   closes [#67](https://github.com/lgatto/synapter/issues/67) [2015-07-26 Sun].
- Use just unique peptides in master creation (see [#107](https://github.com/lgatto/synapter/issues/107)) [2016-01-23 Sat].
- New `rmarkdown` based reports for `synergise1` (synonym to `synergise`) and
   `synergise2`.

### Bugfixes

- Use new loess model in master creation (now based on m-estimator instead of
   least squares, identical to retention time model in classical synergise
   workflow; see [#107](https://github.com/lgatto/synapter/issues/107) for details) [2016-01-23 Sat]
- Fix retention time model calculation in `plotFeatures(..., what="some")`
   [2014-04-28 Mon].

### Internal changes

- Add `testthat` to Suggests [2014-04-25 Fri].
- Add recommended biocView [2014-06-05 Thu].
- Replace `any(is.na(...)` by `anyNA(...)`; *synapter* depends on
   `R >= 3.1.0` now [2014-11-01 Sat].
- Add `ClassVersion` field to `Synapter` class [2014-11-21 Fri].
- Add `Versioned` class as parent class to `MasterPeptides` and
   `MasterFdrResults` [2014-11-22 Sat].
- Adapt `synergise` to new grid search (closes [#81](https://github.com/lgatto/synapter/issues/81)) [2016-10-16 So].
- Replace `hwriter` by `rmarkdown` report in `synergise`; closes [#120](https://github.com/lgatto/synapter/issues/120).
  [2016-10-17 Mon]

### Removed functions/arguments

- Remove `synapterGUI`.
- Remove unused internal functions: `filterCommonSeq`, `filterKeepUniqueSeq`,
   `filterKeepUniqueProt` [2014-11-27 Thu].
- Remove "mergedEMRTs" argument from `findEMRTs`. Now `rescueEMRTs` has to be
   called manually at the end of the processing; close [#93](https://github.com/lgatto/synapter/issues/93) [2015-07-26 Sun]
- Remove "light" version of `writeMergedPeptides` and `writeMachtedPeptides`
   (now always the full `data.frame` is saved; see [#95](https://github.com/lgatto/synapter/issues/95)) [2016-10-16 Sun]
- Update `synapterTiny` and `synapterTinyData` [2016-10-16 So]

# Synapter 1.13

## Changes in version 1.13.1

- Update call to nQuants to accomodate changes in MSnbase (see commit
   #ce637b87e491a2f4c9f93490b2f4c7027cbfce16) <2016-02-29 Mon>
- Defunct synapterGUI <2016-02-29 Mon>

## Changes in version 1.13.0

- Bioc devel 3.2

# Synapter 1.12

## Changes in version 1.12.0

- Bioc release 3.1

# Synapter 1.11

## Changes in version 1.11.2

- fixing bug introduced in 1.11.1 when coercing from
   Synapter to MSnSet object - see this post on the support forum for
   details https://support.bioconductor.org/p/71087/ <2015-09-27 Sun>

## Changes in version 1.11.1

- avoiding error when some fcols are missing when coercing from
   Synapter to MSnSet object - see this post on the support forum for
   details https://support.bioconductor.org/p/71087/ <2015-09-07 Mon>

# Synapter 1.9

## Changes in version 1.9.5

- update cross reference to qvalue::plot.qvalue in Synapter man page;
   closes [#86](https://github.com/lgatto/synapter/issues/86) [2015-03-31 Tue]

## Changes in version 1.9.4

- use biocstyle [2015-01-23 Fri]
- use requireNamespace [2015-01-23 Fri]

## Changes in version 1.9.3

- Filter entries in quantiation final peptides and Pep3D data that don't match
   in their intensity valus; see [#42](https://github.com/lgatto/synapter/issues/42) for details. [2014-11-26 Wed]

## Changes in version 1.9.2

- Only suggesting tcltk and tcltk2, as gui is now deprecated
   [2014-11-21 Fri]

## Changes in version 1.9.1

- Deprecating synatperGUI [2014-11-10 Mon]
- Directing questions to the Bioc support site [2014-11-10 Mon]

## Changes in version 1.9.0

- new devel version for Bioc 3.1

# Synapter 1.7

## Changes in version 1.7.0

- Deprecating synatperGUI [2014-11-10 Mon]
- Directing questions to the Bioc support site [2014-11-10 Mon]

# Synapter 1.5

## Changes in version 1.5.7

- Enable setting number of missed cleavages (closes [#53](https://github.com/lgatto/synapter/issues/53)) [2014-03-24 Mon]

## Changes in version 1.5.6

- change the default value of grid.ppm.from to 2
   (was 5 before) [2014-03-05 Wed]
- change the separator in GridDetails' names to ":"
   (was "." before) [2014-03-07 Fri]
- test for corresponding Pep3D file (closes [#42](https://github.com/lgatto/synapter/issues/42)) [2014-03-19 Wed]

## Changes in version 1.5.5

- fix a bug in the calculation of non unique matches in gridSearch2
   (part of searchGrid); results in a higher number of non unique matches
   (some of the reported -2 will now be reported as 2 in the grid details)
   [2014-03-05 Wed]
- partial rewrite of gridSearch2 (part of searchGrid) for faster grid
   calculation [2014-03-04 Tue]

## Changes in version 1.5.4

- biocViews update

## Changes in version 1.5.3

- modified findMSeEMRTs (part of searchGrid) for faster grid
   calculation [2014-02-28 Fri]
- typos in manual [2014-02-25 Tue]
- replace readFasta by Biostrings::readAAStringSet [2014-02-25 Tue]
- replace digest by cleaver::cleave [2014-02-25 Tue]

## Changes in version 1.5.1

- typo in vignette [2014-02-18 Tue]

## Changes in version 1.5.0

- new devel version for Bioc 2.14

# Synapter 1.3

## Changes in version 1.3.4

- estimateMasterFdr now support list of vector as pepfiles [2013-09-27 Fri]
- import(MSnbase) [2013-09-27 Fri]

## Changes in version 1.3.3

- Updated references. [2013-07-05 Fri]

## Changes in version 1.3.2

- fixed bug when using findEMRTs = "copy" [2013-05-30 Thu]

## Changes in version 1.3.1

- Reporting total number of peptides in dbUniquePeptideSet.
   Fixes issue [#41](https://github.com/lgatto/synapter/issues/41). [2013-05-13 Mon]
- New mergedEMRTs arg in findEMRTs.
   Closes issue [#38](https://github.com/lgatto/synapter/issues/38). [2013-05-13 Mon]
- fixed synapterTiny\$QuantPep3DData, which had the rownames
   as first column synapterTiny\$QuantPep3DData\$X. Detected
   thanks to new mergedEMRTs arg. [2013-05-13 Mon]
- added mergedEMRTs arg to synergise [2013-05-22 Wed]
- Synapter checks that one file per list element is passed [2013-05-28 Tue]
- minor typo/fixes [2013-05-28 Tue]
- new idSource column when matching EMRTs [2013-05-29 Wed]
- new performance2 method that shows identification source
   and NA values contingency table [2013-05-29 Wed]
- new filterPeptideLength method [2013-05-29 Wed]
- added peplen argument to synergise to filter
- n peptide length [2013-05-29 Wed]

## Changes in version 1.3.0

- Bioc 2.13 devel version bump

# Synapter 1.2

## Changes in version 1.2.0

- Bioc 2.12 stable version bump

# Synapter 1.1

## Changes in version 1.1.5

- updated references [2013-03-22 Fri]

## Changes in version 1.1.4

- added citation [2013-03-21 Thu]
- vignette uses knitr engine and scrartcl class [2013-03-21 Thu]

## Changes in version 1.1.3

- knitr 1.0 compatibility, based on Dan's
   updates [2013-01-15 Tue]

## Changes in version 1.1.2

- added note about tcl BWidget package in synapterGUI man [2012-10-04 Thu]

## Changes in version 1.1.1

- fixing vignette [2012-10-02 Tue]

## Changes in version 1.1.0

- new devel version bump [2012-10-01 Mon]

# Synapter 0.99

## Changes in version 0.99.15

- updated ref in package man (2) [2012-09-14 Fri]
- updated Synapter show method to display short
   log [2012-09-14 Fri]

## Changes in version 0.99.14

- updated ref in package man [2012-09-11 Tue]

## Changes in version 0.99.13

- bumping version number to take all previous changes
   into account [2012-08-29 Wed]

## Changes in version 0.99.12

- fixed build problem on windows [2012-08-27 Mon]
- updated refs in vignette [2012-08-28 Tue]

## Changes in version 0.99.11

- using dev='pdf' in vignette and removing
   tikzDevice (not on CRAN anymore) from
   Suggests[2012-08-16 Thu]

## Changes in version 0.99.10

- added tikzDevice suggestion for vignette [2012-08-15 Wed]

## Changes in version 0.99.9

- using knitr instead of pgfSweave [2012-08-13 Mon]
- new verbose arg to loadIdentOnly taken into account
   by estimateMasterFdr and makeMaster [2012-08-13 Mon]

## Changes in version 0.99.8

- vignette update with summary table [2012-07-17 Tue]
- vignette update with additional PLGS slides [2012-07-17 Tue]

## Changes in version 0.99.7

- updated vignette [2012-07-12 Thu]
- removed README.org file [2012-07-12 Thu]
- added github page in DESCRIPTION [2012-07-12 Thu]

## Changes in version 0.99.6

- Masterpeptides has new fdr and method slots  [2012-07-10 Tue]
- added method arg to makeMaster [2012-07-10 Tue]
- fixed non-match precuror.leIDs when using master [2012-07-11 Wed]

## Changes in version 0.99.5

- removed 'precursor.Mobility' columns from light
   merged and matched csv output, as these are not
   available when a MSe master is used [2012-07-09 Mon]

## Changes in version 0.99.4

- use log2 for t.test [2012-07-05 Thu]
- updates to vignette [2012-07-06 Fri]

## Changes in version 0.99.3

- git/svn integration messing up and testing [2012-07-03 Tue]
- late night testing [2012-07-03 Tue]

## Changes in version 0.99.2

- changed merged to total FDR in MasterFdrResults
   plot x label, as in manuscript [2012-06-13 Wed]
- mention multtest in vignette [2012-06-13 Wed]
- qvalue refs [2012-06-14 Thu]
- Synapter.Rd updates, to match manuscript
   nomenclature [2012-06-14 Thu]
- Added 'Getting help' section in vignette and
   package man page [2012-06-14 Thu]
- Changed HDMSe to ident in plotGrid plot
   legend [2012-06-25 Mon]

## Changes in version 0.99.1

- Latex typo in vignette [2012-06-11 Mon]
- Dan's suggestions for [2012-06-13 Wed]

## Changes in version 0.99.0

- bump version for Bioconductor submission [2012-06-08 Fri]

## Changes in version 0.9.1

- more doc completion [2012-06-07 Thu]
- fixed type in dim [2012-06-07 Thu]
- added PLGS ref and screenshots to
   vignette [2012-06-07 Thu]
- finished methods documentation in
   Synapter.Rd [2012-06-07 Thu]
- reduced screenshots sizes to get
   below 2M [2012-06-08 Fri]
- fixed synapterGUI bug introduced
   during naming refactoring [2012-06-08 Fri]

## Changes in version 0.9.0

- moved repo to github [2012-06-06 Wed]
- changed HDMSe/MSe naming convention
   to identification/quantitation [2012-06-06 Wed]
- second round of naming changes (utils) [2012-06-07 Thu]
- Description of package and methodology updated and
   references [2012-06-07 Thu]

## Changes in version 0.8.14

- Updated package title [2012-06-05 Tue]
- removed PLGS vignette [2012-06-06 Wed]

## Changes in version 0.8.13

- more vignette work [2012-06-01 Fri]
- added pgfSweave and xtable to Suggests
   and loaded in vignette [2012-06-01 Fri]

## Changes in version 0.8.12

- updated show MasterPeptides to reflect number of files
   effectively used [2012-05-28 Mon]
- moved MasterFdrResults to S4 framework [2012-05-30 Wed]
- MSnbase is now a dependency [2012-05-30 Wed]
- vignette updates [2012-05-30 Wed]
- synergise now calls filterUniqueDbPeptides even if master,
   as opposed to filterUniqueMSeDbPeptides [2012-05-30 Wed]
- fixes as(, "MSnSet") - peptides sequences are now
   properly made unique [2012-05-30 Wed]

## Changes in version 0.8.11

- Automatic handling of Masterpeptides objects [2012-05-23 Wed]
- Various man updated to describe the above [2012-05-23 Wed]
- Now adding vignettes to windows vignette menu [2012-05-23 Wed]
- fixed/improved writeMasterPeptides and fixed [tk_]choose.dir
   in synergise [2012-05-23 Wed]
- MasterPeptides must be save as 'rds' with saveRDS to be
   reusable in synergise/Synapter [2012-05-23 Wed]

## Changes in version 0.8.10

- updated master-related function names by removing
   explicit reference to HDMSe [2012-05-16 Wed]
- new AllGenerics.R and AllClasses.R files [2012-05-16 Wed]
- removed selection param to makeMaster [2012-05-16 Wed]
- new makeMaster returns an instance of class
   'MasterPeptides' [2012-05-16 Wed]

## Changes in version 0.8.9

- new synergise verbose parameter [2012-05-11 Fri]
- filterUnique[HDMSe|MSe]DbPeptides new verbose param [2012-05-11 Fri]
- started method documentation in Synapter.Rd  [2012-05-11 Fri]
- import Vennerable::Venn [2012-05-15 Tue]
- updated makeMasterHDMSe file - now has a 'selection' parameter,
   fixes issue [#27](https://github.com/lgatto/synapter/issues/27). [2012-05-15 Tue]

## Changes in version 0.8.8

- man/vignette updates [2012-05-02 Wed]
- added PLGS data input to vignette, contributed by Pavel [2012-05-02 Wed]
- added Pavel's PLGS_Data_Processing.pdf slides as vignette [2012-05-03 Thu]
- added Biobase import to benefit of addVigs2WinMenu [2012-05-03 Thu]
- new plgsGuide() function to access the PLGS_Data_Processing vignette [2012-05-03 Thu]
- fixed gui ppm.by label (issue [#25](https://github.com/lgatto/synapter/issues/25)) [2012-05-03 Thu]
- License is GPL-2 [2012-05-03 Thu]
- vignette improvements [2012-05-03 Thu]

## Changes in version 0.8.7

- added Nick's email to vignette [2012-04-30 Mon]
- worked on gui - tree widget now expands when resizing window [2012-04-30 Mon]
- added makeMasterHDMSe and estimateMasterFdr cross
   references [2012-05-01 Tue]
- continued vignette [2012-05-01 Tue]
- added synapterTiny test data [2012-05-01 Tue]
- synapterTinyData to load and initialise the test data [2012-05-01 Tue]
- synapterGuide sortcut to open vignette [2012-05-01 Tue]
- synergise example (notrun, though) [2012-05-01 Tue]

## Changes in version 0.8.6

- Testing if model exists before findEMRTs/searchGrid [2012-04-27 Fri]
- added fdr method choice in synergise [2012-04-27 Fri]
- synergise code updates and documentation [2012-04-27 Fri]
- work on vignette and new bib file [2012-04-27 Fri]
- added package documentation and updated start-up message [2012-04-27 Fri]
- added utils in imports [2012-04-27 Fri]

## Changes in version 0.8.5

- merged filterUniqueDbPeptides and filterUniquePeptides into
   filterUniqueDbPeptides in the Synapter class interface
   (same for MSe/HDMSe specific methods). Both are still distinct
   at the RefClass level though. [2012-04-27 Fri]
- TODO - maybe add fdr method choice in synergise - DONE v 0.8.6

## Changes in version 0.8.4

- new inspectPeptideScores function [2012-04-26 Thu]
- finished gui code - synapterGUI exported [2012-04-27 Fri]
- tested synapterGUI until report generation [2012-04-27 Fri]
- TODO - avoid loading fasta db twice in filterUniqueDbPeptides

## Changes in version 0.8.3

- fixes getIdStats - now stops if any NA are found [2012-04-25 Wed]

## Changes in version 0.8.2

- fixed typo in data.frame column name [2012-04-25 Wed]

## Changes in version 0.8.1

- new masterFdr function that computes merged FDR - still
   under construction [2012-04-06 Fri]
- updated Synapter.Rd [2012-04-23 Mon]
- Only considered *unique peptides* (sequences) to calculate the p-values.
   The .[HD]MSePeptideScores dataframes are also updated so that
   plotPepScores and getPepNumbers also use unique() peptides to visualise
   the scores. [2012-04-23 Mon]
- estimateMasterFdr (was masterFdr) is functional [2012-04-24 Tue]
- estimateMasterFdr output stored in MasterFdrRes S3 class [2012-04-24 Tue]
- TODO Investigate NA p-values - DONE - bug fixed in v 0.8.3

## Changes in version 0.8.0

- added separate add[HD]MSeIdStats method [2012-04-06 Fri]
- added separate filterUnique[HD]MSeDbPeptides method [2012-04-06 Fri]
- added separate filterUnique[HD]MSeSeq method [2012-04-06 Fri]
- new makeMasterHDMSe function [2012-04-06 Fri]
- new Master field in .Synapter class [2012-04-06 Fri]
- new loadMasterData function [2012-04-06 Fri]
- modifed plotFdr and plotPepScores if master [2012-04-06 Fri]
- updated Synapter constructor for master [2012-04-06 Fri]
- new setMaster method [2012-04-06 Fri]
- renames filterUniqueDBpeptides to filterUniqueDbPeptides
   and added filterUnique[HD]MSeDbPeptides methods in interface [2012-04-06 Fri]
- fixes getPepNumbers for master HDMSe [2012-04-06 Fri]

## Changes in version 0.7.8

- renamed plotQvalues to plotFdr [2012-04-05 Thu]
- new getFdrStats method [2012-04-05 Thu]

## Changes in version 0.7.7

- non-unique matched features now get NAs except
   for the first column (Function) [2012-04-04 Wed]
- changed FDR adjustment; it is now possible to choose between
   'BH' (default), 'Bonferroni' and 'qval' - see issue [#22](https://github.com/lgatto/synapter/issues/22) [2012-04-04 Wed]
- updated plotQvalues accordingly [2012-04-04 Wed]
- fix in setBestGridParams (== to %in%) to avoid warnings
   when different nrow for prctModel (x) and details (y)
   best grid params [2012-04-04 Wed]

## Changes in version 0.7.6

- default span set to 0.05 [2012-03-26 Mon]
- fix in setBestGridParams [2012-03-27 Tue]

## Changes in version 0.7.5

- fixed bug in getBestGridParams when what != "auto"
   and mutliple solutions [2012-03-26 Mon]
- type in Synapter.Rd [2012-03-26 Mon]
- changed ppm and nsd defaults in searchGrid to reflect
   defaults in synergise and gui [2012-03-26 Mon]
- new 'object' parameter in synergise [2012-03-26 Mon]
- updated inputFiles accessor (names and order) [2012-03-26 Mon]
- fixed bug in searchGrid when using 'n' [2012-03-26 Mon]

## Changes in version 0.7.4

- using "model" as default parameter selection method [2012-03-23 Fri]
- new "auto" parameter selection model - uses "model" and
   if multiple pairs, uses best "details" if any of ppm/nsd
   match [2012-03-23 Fri]
- remove parameter selection method in gui - using auto [2012-03-24 Sat]
- update html report - provide ppm and nsd values selected
   from grid search and auto method [2012-03-24 Sat]
- updates to Synapter.Rd [2012-03-25 Sun]

## Changes in version 0.7.3

- dropped synergize spelling; too annoying when autocompleting [2012-03-22 Thu]
- when computing grid details, checking if x["1"] or x["-1"] are NA and
   replacing by 0 [2012-03-22 Thu]
- fixed bug in detail grid [2012-03-22 Thu]
- now setting best params in synergize [2012-03-23 Fri]

## Changes in version 0.7.2

- synergi[s/z]e continued [2012-03-21 Wed]
- also create a grid based on grid details.
   This actually fixes a bug in getBestGridParams
   'details' search [2012-03-21 Wed]
- Removed filterUniquePeptides from loadData, forgot
   to do that in version 0.7.1 [2012-03-21 Wed]
- finished synergi[s|z]e [2012-03-22 Thu]

Changes IN VERSION 0.7.1

- synergi[s/z]e definition [2012-03-08 Thu]
- GUI implementation [2012-03-10 Sat]
- cleaned performance method [2012-03-20 Tue]
- new filterUniquePeptides method that filter peptides	that
   appear multiple times in the respective HDMSe/MSe input
   peptide files. Previously, this was done when data was
   loaded. Also updated synergise, example and doc. [2012-03-21 Wed]

## Changes in version 0.7.0

- using filtered HDMSePeptideData (was mfHDMSePeptideData)
   to do the matching; mf[HD]MSePeptideData have been
   dropped for Synapter class  [2012-03-07 Wed]
- new dims method [2012-03-07 Wed]
- best grid parameter selection based on +1/(+1 + -1) proportion
-f correct unique identifications (issue [#18](https://github.com/lgatto/synapter/issues/18)) [2012-03-07 Wed]

## Changes in version 0.6.0

- fixed bug when providing filenames to Synapter [2012-03-06 Tue]
- added 'n' parameter to searchGrid to specify a number of
   features to be used for the grid search [2012-03-06 Tue]
- Test if n ] nrow(HDMSe peptides), warn and use all HDMSe peptides [2012-03-06 Tue]
- moved some parameter initialisation from searchGrid method to
   searchGrid interface [2012-03-06 Tue]
- new filterUniqueDbPeptides method; this is NOT done anymore when
   the data is loaded [2012-03-06 Tue]

## Changes in version 0.5.3

- new write[HD]MSePeptides methods [2012-02-17 Fri]

## Changes in version 0.5.2

- added require(MSnbase) in setAs [2012-02-14 Tue]
- removed pre/post-fix peptides - is not an issue [2012-02-14 Tue]

## Changes in version 0.5.1

- new setAs method to coerce a Synapter instance into
   an MSnSet object [2012-02-14 Tue]
- added MSe prot FPR filter when calculating performance
   [2012-02-14 Tue]

## Changes in version 0.5.0

- filtering HDMSe unique peptides [2012-02-12 Sun]
- added explicit values to be set for setPpmError
   and setRtNsd in Synapter man example [2012-02-13 Mon]
- filtering pep score (q-val) AND prot FPR (NEW) before
   matching EMRTs and calculating performance [2012-02-13 Mon]
- filtering HDMSe unique peptides [2012-02-13 Mon]
- fasta file is also asked for when creating instance [2012-02-13 Mon]
- Filtering unique (i.e. non duplicated peptide sequences) is now done
   when data is loaded, after db unique peptide filtering, to get rid
-f pre/post-fix cases (like ABCD and ABCDXXXX) that are not removed
   by filterUniqueDbPeptides/dbUniquePeptideSet. [2012-02-13 Mon]
- keeping common HDMSe/MSe peptides filter is not included in
   merge Peptides (filterPepSeq removed) [2012-02-13 Mon]
- removed (private) filter and idStats parameters [2012-02-13 Mon]
- updated documentation [2012-02-13 Mon]
- updated performace; removed duplicates count as there are no
   duplicated sequences/precursor.leIDs anymore [2012-02-13 Mon]
- removing duplicated assigned EMRTs; these rows in the matched
   results dataframes are identified with a -1 [2012-02-13 Mon]
- new getLog method [2012-02-14 Tue]
- updaed show method [2012-02-14 Tue]

Changes IN VERSION 0.4.7

- fixed 'performance' typo in Synapter man [2012-02-09 Thu]
- fixed bug in performance/getBestGridValue (issues [#5](https://github.com/lgatto/synapter/issues/5) and [#6](https://github.com/lgatto/synapter/issues/6)) [2012-02-10 Fri]
- Updated performance method [2012-02-10 Fri]
- Added p/q-values to 'light' merged output [2012-02-10 Fri]

## Changes in version 0.4.6

- HDMSePeptideData is updated (predictRt and sdRt columns) when
   rt is modelled. [2012-02-09 Thu]
- Updated doHDMSePredictions to either predict rt using the model
   (as previously) or retrieve that data from the HDMSePeptideData
   data.frame (new). [2012-02-09 Thu]
- doHDMSePredictions now returns a list. [2012-02-09 Thu]
- added a subet parameter to searchGrid [2012-02-09 Thu]
- added a performace method (see issue [#1](https://github.com/lgatto/synapter/issues/1)) [2012-02-09 Thu]

## Changes in version 0.4.5

- removed extra plotRtDiff un Synapter example [2012-02-07 Tue]
- findMSeEMRTs tidy-up, added a matched.mse.spectrumIDs column to
   the MatchedEMRTs data.frame; this stores the MSePep3D spectrumIDs
   that were macthed within an ErrorPpm by RtNsd window. [2012-02-07 Tue]
- added writeMergedPeptides/writeMatchedEMRTs file argument
   examples. [2012-02-07 Tue]
- PpmError and RtNsd checked/set only for plotFeatures(, what = "some")
   and not what = "all" anymore. [2012-02-07 Tue]
- Added the recursor.leID.mse column to the MatchedEMRTs data.frame
   with merged MSe IDs corresponding to the HDMSe IDs used for rt
   modelling [2012-02-07 Tue]
- grd2 is now computed using new MatchedEMRTs columns [2012-02-07 Tue]
- function to get assignment details [2012-02-08 Wed]
- code clean-up: xx\$HDMSePeptideData is now filtered for peptide score fdr
   explicitely by the top-level caller (in synapter-class), rather than
   (redundantely) inside findMSeEMRTs and in gridSearch2 (which also calls
   the latter). Removed fdr parameter in these functions as unused.
   This is also documented now. [2012-02-09 Thu]
- code clean-up: remove old gridSearch call and definition [2012-02-09 Thu]

## Changes in version 0.4.4

- fixed writeMatchedEMRTs/writeMergedPeptides(..., what = "full") [2012-02-07 Tue]
- Synapter constructor has filenames argument [2012-02-07 Tue]

## Changes in version 0.4.3

- continued documentation [2012-01-31 Tue] [2012-02-01 Wed]
- added plotRtDiffs in example [2012-02-01 Wed]
- added warning and nsd[1:3] selection in plotRt model [2012-02-01 Wed]
- fixed bug in filterPepSeq; wring MSe dimensions
   where writton to log [2012-02-03 Fri]
- new getEMRTtable method [2012-02-03 Fri]
- added ... to plotRtDiffs to be passed to hist [2012-02-03 Fri]

## Changes in version 0.4.2

- Matched EMRTs dataframe cols class is set explicitly based
-n HDMSePeptideData and MSePep3DData cols classes [2012-01-30 Mon]
- Grid search on percentage of correct assignments [2012-01-30 Mon]
- updated grid getters/setters accordingly [2012-01-31 Tue]

## Changes in version 0.4.1

- added legend to plotFeature(..., what="some") [2012-01-27 Fri]
- added protein false positive rate functionality [2012-01-27 Fri]

## Changes in version 0.4.0

- introducing ReferenceClasses implementation [2012-01-13 Fri]
- filterPepScore takes medianRndScore and stdRndScore
   as separate parametres instead as a list (params) [2012-01-14 Sat]
- fixed a typo in makeSynapterReport [2012-01-16 Mon]
- changed several function parameters from input$... into ... [2012-01-17 Tue]
- further changes to gridSearch, doHDMSePredictions and findMSeEMRTs
   to integrate with ref classes implementation [2012-01-17 Tue]
- rationalised findMSeEMRTs and gridSearch 'sapply' loops [2012-01-17 Tue]
- new keepUniqueSpectrumIds function, now called explicitely
   before gridSearch and findMSeEMRTs [2012-01-18 Wed]
- Testing that changes produce same makeSynapterReport results
   than version 0.3.2 [2012-01-17 Tue]
- added identification p/q-values calculation and qvalue and multtest
   in imports [2012-01-23 Mon]
- global filtering is now done when data is loaded; can also be done
   manually [2012-01-23 Mon]
- removing PLGS pepscore filtering and HDMSe log input [2012-01-24 Tue]
- added qvalue and peptide score plots [2012-01-24 Tue]
- Synapter ref class interface using S4 methods [2012-01-25 Wed]
- Doing mass error filtering on the individual peptide files
   instead of merged data. [2012-01-26 Thu]
- mf[HD]MSePeptideData is not initialised when files are loaded [2012-01-26 Thu]
- Documentating new interface [2012-01-26 Thu]
- added getRtQs and plotRtDiffs [2012-01-27 Fri]

## Changes in version 0.3.2

- save intermediate filtering tables [2011-12-23 Fri]
- added ... to plotLowess2 [2012-01-09 Mon]

## Changes in version 0.3.1

- changed figures url to relative [2011-12-22 Thu]
- added synapter version and date to report output [2011-12-22 Thu]
- updated filter.error.ppm to specify column name
   to use as filter [2011-12-22 Thu]
- new mass error plot [2011-12-22 Thu]
- added orignal nrows in filtering summary table [2011-12-22 Thu]
- when provided, htmldir is now created automatically [2011-12-22 Thu]

## Changes in version 0.3.0

- modified rt range calculation [2011-12-20 Tue]
- modified feature plot to incorporate rt and mass ranges [2011-12-20 Tue]
- adjusted lowess plots ylims to rt modelling range [2011-12-20 Tue]
- saving light result tables [2011-12-21 Wed]
- added option to perform grid search over ppm and nsd
   to maximise the number of unique matches [2011-12-21 Wed]
- added full range feature plots [2011-12-21 Wed]

## Changes in version 0.2.0

- fixed bug when input was not called 'input' [2011-12-16 Fri]
- can now set output directory in makeSynapterReport [2011-12-16 Fri]
- using CairoSVG [2011-12-16 Fri]

## Changes in version 0.1.0

- initial release [2011-12-15 Thu]
