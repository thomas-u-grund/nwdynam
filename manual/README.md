# The nwdynam manual

`nwdynam_manual.pdf` is the user manual and reference of nwdynam. The LaTeX source is `nwdynam_manual.tex` with one file per chapter (`ch*.tex`, `app_translation.tex`).

## Bumping the version (one-line change)

The manual's title and citation carry no version number; the title page states the series it describes (`\manualseries`, valid for 1.8.x). The exact version and date of the release that produced the listings are two more macros in `nwdynam_manual.tex`:

    \newcommand{\manualseries}{1.8}
    \newcommand{\pkgversion}{1.8.0}
    \newcommand{\pkgdate}{16 November 2026}

`tools/set_version.py` of the package does not write these macros; change them by hand. Change `\manualseries` only for a new minor or major release whose documented behaviour changed. For a patch release change the last two lines, rerun the examples (`make examples`) if their output changes, rebuild (`make pdf`) and check. The tag in the installation example of `ch01_intro.tex` (`nwdynam/v1.8.0/`, inside a listing, so it cannot use the macro) is the only other place with the number.

## The DOI and the draft marks

`\manualdoi` in `nwdynam_manual.tex` is a placeholder (`10.5281/zenodo.TO-BE-REPLACED`): the Zenodo record does not exist yet. Replace it by the all-versions DOI of the record (the owner's rule: the all-versions DOI everywhere) and add it to the citation texts listed in `ZENODO.md`.

The red notes `TODO(owner)` in the PDF mark facts that could not be verified and decisions that are the owner's. They are printed while `\drafttrue` is set in `nwdynam_manual.tex`; set `\draftfalse` for the version that is archived, after every note has been dealt with (`grep -n todoowner *.tex` lists them).

## Building

    make -C manual pdf        # regenerates the comparison table, then pdflatex (3 passes)
    make -C manual check      # undefined references, overfull boxes, glued words
    NWCORE=/path/to/nwcommands make -C manual examples   # reruns the do-files of ex/ in Stata 17 (a few minutes; the data are loaded from the web)
    NWCORE=/path/to/nwcommands make -C manual compare    # the numerical comparison with goldfish (needs R with goldfish 1.6.12)

`ex/mklistings.py` cuts the Stata logs of `ex/*.do` at the marker lines `*# name` into the listings `ex/out/name.txt` that the chapters include (`\showlog{name}`); the do-files start with `ex/setup.doh` (empty PLUS folder, the core, the package and its plugin on the adopath). The examples run with Stata 17 (`run_examples.sh`; `STATA=...` selects another binary), at most one at a time here. `build_compare.py` reads `ex/r/goldfish_results.csv` (written by `ex/r/compare.R`) and `ex/r/nwdynam_results.csv` (written by `ex/ex91_compare_stata.do`) and writes `gen_compare.tex` and `gen_compare_nums.tex`; do not edit those two by hand. `ex/data/socialevolution_calls.csv` is the calls of the data set `Social_Evolution` of goldfish.

The PDF is committed with `git add -f` (LaTeX build files and Stata logs are ignored). The Stata logs `ex/*.log` are deliberately not committed: they carry the licence header (name and serial number) of the Stata that ran them. The committed listings `ex/out/*.txt` are cut from them.

## Archiving

`ZENODO.md` holds the metadata for the Zenodo record and the steps to do.
