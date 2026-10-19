#!/usr/bin/env bash
# Builds the browser-only edition of the Python application (Shinylive: Python compiled to
# WebAssembly, no server) into site/. Serve site/ from any static host, e.g. GitHub Pages.
set -euo pipefail
cd "$(dirname "$0")"
rm -rf _site_src site
mkdir -p _site_src/templates
cp edna_py/app.py edna_py/original.py _site_src/
cp app_R/template_sampling_perSize.csv app_R/template_source.csv _site_src/templates/
printf 'numpy\npandas\nscipy\nmpmath\nmatplotlib\n' > _site_src/requirements.txt
shinylive export _site_src site
rm -rf _site_src
echo "site/ built; test with: python -m http.server --directory site 8008"
