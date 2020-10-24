# ImageJ-plugins

A set of ImageJ Jython plugins for molecular biology and bioimaging applications.

*western_processor.py* rotates and crops a set of images to specifications given based on a reference image. Designed for use on blots, gels, or any other image type where side-by-side comparison of multiple exposures is common.

*subcell_loc.py* measures the area and intensity of the cytosol and nucleus of cells in the given multichannel czi image. Will process one multiframe czi image at a time. Designed for quantification of the subcellular location of labelled protein(s).

*F2H_processing.py* measures the area and intensity of the LacO array (or similar relevant tethering method) and nucleoplasm of cells in a fluorescent two-hybrid assay. Will run on a folder of multichannel czi images.
