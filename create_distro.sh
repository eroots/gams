#!/bin/bash

pyinstaller -y --clean --onedir \
--add-data "gams/GUI/gams_gui.ui;./gams/GUI" \
gams.spec
# gams.py