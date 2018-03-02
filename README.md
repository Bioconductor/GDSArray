# GDSArray
Representing GDS files as array-like objects

GDS files are widely used to represent genotyping or
sequence data. The GDSArray package implements the `GDSArray`
class to represent nodes in GDS files in a matrix-like
representation that allows easy manipulation (e.g., subsetting,
mathematical transformation) in _R_. The data remains on disk
until needed, so that very large files can be processed.
