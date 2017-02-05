#!/usr/bin/env python
# -*- coding: utf-8 -*-

about = {
    "name":             "Multipool",
    "version":          "0.10.2",
    "description":      "Efficient multi-locus genetic mapping with pooled sequencing.",
    "author":           "Matt Edwards",
    "author_email":     "matted@mit.edu",
    "license":          "MIT",
    "url":              "https://github.com/matted/multipool",
    "packages":         [],
    "scripts":          ["mp_inference.py", "mp_prep.py"],
    "zip_safe":         True,
    "install_requires": ["scipy", "numpy"] # pylab is optional; leaving it out for now
}
