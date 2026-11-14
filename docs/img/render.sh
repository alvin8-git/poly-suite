#!/usr/bin/env bash
# Render every Mermaid source here to SVG. The .mmd files are the source of truth.
cd "$(dirname "$0")" && for f in *.mmd; do PUPPETEER_SKIP_DOWNLOAD=1 npx -y @mermaid-js/mermaid-cli@11 --input "$f" --output "${f%.mmd}.svg" --puppeteerConfigFile pconf.json --configFile mmdconf.json --backgroundColor transparent; done
