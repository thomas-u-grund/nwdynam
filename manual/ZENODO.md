# Zenodo record for the nwdynam manual

Metadata ready to paste into a new Zenodo upload. Nothing has been uploaded and the record does not exist yet.

**DOI.** To be created by the owner. `\manualdoi` in `nwdynam_manual.tex` holds the placeholder `10.5281/zenodo.TO-BE-REPLACED`; after the record is reserved, replace it by the all-versions ("Cite all versions") DOI and tell it to the agent. The citation text always uses the all-versions DOI (master RELEASING.md B.10).

**Labelling.** The title and the citation of the manual contain no version number. The Zenodo version label is the series "X.Y" (here "1.8"; internal to Zenodo, not part of the citation), and the description says that the manual describes nwdynam X.Y.x. The PDF title page has one line saying which series it describes (`\manualseries`); the exact package version of the examples is `\pkgversion`. The software version is cited separately, as printed by `nwdynam version`. A patch release of the package (1.8.1, ...) does not get a new Zenodo version.

**Patch releases.** The manual of the series 1.8 holds for 1.8.x; a patch release changes no documented option or number. Title and citation stay unchanged; no new Zenodo version is needed for a patch.

## Metadata

**Title.** nwdynam: Dynamic network actor models in Stata. Manual

**Version.** 1.8 (the series label, internal to Zenodo: the manual describes nwdynam 1.8.x; not part of the citation)

**Publication date.** The day of the upload.

**Authors.** Grund, Thomas U. (RWTH Aachen University), ORCID 0000-0001-6751-5429 (as on the nwergm record proposal; remove it in the form if the author does not want it).

**Description** (plain text, one paragraph):

> This is the user manual and reference of nwdynam, series 1.8 (examples produced with version 1.8.0 and the nwcommands core 3.26.0), a Stata package for dynamic network actor models (DyNAMs) of relational event data. nwdynam builds on the nwcommands core package and follows the R package goldfish (version 1.6.12) in its model, the definitions of its effects, its likelihoods and its postestimation tools. The manual explains the model (the rate sub-model with the ordinal partial likelihood and with an intercept as a continuous-time hazard, the choice sub-model, and the coordination sub-model for undirected events), the declaration of the event data (event values, events with identical times, two-mode networks, covariates that change over time, actors who join or leave, restricted candidate sets), the syntax of every option, the definition of every effect (structural, covariate, windowed, weighted, and effects of other networks including mixed closures), the postestimation commands (information criteria, log likelihoods and probabilities of events, outliers, change points), worked examples with output, the comparison with goldfish, the limits of the package, the technical background (a native plugin that evaluates the likelihood, the gradient and the Hessian) and a translation table from goldfish to nwdynam. The do-files that produce every listing are in the repository. The package and the manual were prepared with the help of an AI coding assistant working under the author's direction, as stated in the manual.

**License.** Creative Commons Attribution 4.0 International (CC-BY-4.0), as on the nwsaom manual record (the record there). The code stays under GPL-3.0-only.

**Resource type.** Publication, sub-type Software documentation (as the nwsaom record).

**Creators.** As on the nwsaom record, with the ORCID (see Authors).

**Keywords.** Stata; nwcommands; nwdynam; dynamic network actor model; DyNAM; relational events; social network analysis; statistical network models; event history; goldfish

**Language.** English

**Related identifiers.** (relation, identifier, type)

| Relation | Identifier | Type |
|---|---|---|
| Is supplement to | https://github.com/thomas-u-grund/nwdynam/releases/tag/v1.8.0 | URL |
| Is supplement to | https://github.com/thomas-u-grund/nwdynam | URL |
| Requires | https://github.com/thomas-u-grund/nwcommands | URL |
| Is related to (manual of a sibling package) | 10.5281/zenodo.23623807 | DOI |
| Is related to (manual of a sibling package) | 10.5281/zenodo.23661331 | DOI |

Confirm that the GitHub release v1.8.0 is published before the upload, because the first related identifier points to it. The sibling-manual links are optional; they are sensible because the manuals describe packages of the same family (10.5281/zenodo.23623807 is the nwsaom manual, 10.5281/zenodo.23661331 the nwergm manual).

**Other fields.** Funding: none. Communities: none. Access: open.

## Steps for the owner

1. Read the red `TODO(owner)` notes in the draft PDF (`grep -n todoowner manual/*.tex`) and settle each: the Windows/Linux sentence, the origin of the Social Evolution data, the event-network second-network fits (a possible defect, see the report that came with this draft), `estat changepoints, plot`, the timing sentence, size limits. Then set `\draftfalse`.
2. Create the record on zenodo.org (reserve the DOI), tell the all-versions DOI, replace the placeholder in `\manualdoi`, rebuild (`make -C manual pdf check`).
3. Upload `manual/nwdynam_manual.pdf`, fill the fields above, publish.
4. Check that the all-versions DOI of the record is the one in `\manualdoi`.
5. Commit the rebuilt PDF (`git add -f manual/nwdynam_manual.pdf`) on the branch that carries the manual.
6. For a later patch release (1.8.1, ...) nothing is uploaded: only `\pkgversion` and `\pkgdate` change, the examples are rerun if their output changes, and the PDF is committed. A new Zenodo version only for 1.9 or 2.0 with changed documented behaviour: then change `\manualseries` too.

## Suggested edits after the DOI exists (not made)

1. `README.md` (repository), section Manual (new): "Archived at Zenodo: https://doi.org/<DOI>".
2. `CITATION.cff`: add a `preferred-citation`: type manual, title "nwdynam: Dynamic network actor models in Stata. Manual", author Grund, Thomas U., year 2026, publisher Zenodo, doi <DOI>.
3. `NOTICE.md`, section How to cite: add "Grund, T. U. (2026). nwdynam: Dynamic network actor models in Stata. Manual. Zenodo. https://doi.org/<DOI>".
4. Website `nwdynam/citation.md`: the sentence "The package has no DOI" gets the manual's citation and DOI; `citing.md`: nwdynam leaves the list of packages without a DOI and gets its own entry; `index.md` and `nwdynam/index.md`: the "Manual" menu entry for nwdynam (a copy of `website/nwergm/manual.md` with the nwdynam DOI, own `nav_order`, permalink `/nwdynam/manual/`); `news.md`: a dated entry "nwdynam manual archived at Zenodo".
5. RELEASING.md B.10 says that nwdynam and nwrem have no manual; update that sentence.
