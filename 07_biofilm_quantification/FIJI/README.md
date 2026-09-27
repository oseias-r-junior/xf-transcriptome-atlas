# FIJI/masks — user-provided input, not included in this repository

The biofilm-ring quantification scripts in `07_biofilm_quantification/`
(`recompute_from_embedded_roi.py`, `diagnose_threshold_coverage.py`,
`normalize_by_flask_width.py`) read their input photographs from
`FIJI/masks/*.tif` relative to this directory.

**These images are the authors' original manual measurements and are not
distributed with this repository.** To reproduce the biofilm quantification
pipeline, place your own flask-wall photographs here as 8-bit grayscale
TIFFs (`Image > Type > 8-bit` in Fiji/ImageJ), one per sample, with the
ROI used for measurement embedded in each TIFF's metadata (the standard
result of drawing an ROI in Fiji's ROI Manager and saving the image with
"Save" rather than "Export").

Expected layout:

    FIJI/
        masks/
            <label>.tif   # one 8-bit grayscale image per flask photo,
                           # ROI embedded (ImageJ "Iout" IJMetadata tag)

Filenames become the `Label` used throughout
`data/biofilm_ring_area_results*.csv` and `results/biofilm_quantification/`.

See `07_biofilm_quantification/recompute_from_embedded_roi.py`'s own
docstring for the full explanation of the embedded-ROI format and why it
is used instead of a manually-specified `--roi` rectangle.
