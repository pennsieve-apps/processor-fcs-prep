# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following
third-party components, installed from PyPI at the versions pinned in `requirements.txt`.
All are permissively licensed; nothing GPL is linked in-process.

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| FlowKit | 1.3.2 | Reads FCS events (gain/log decoding), applies compensation (`Sample`, `Matrix`) | BSD-3-Clause | https://github.com/whitews/FlowKit |
| FlowIO | 1.4.0 | Reads FCS headers, writes FCS 3.1 (`create_fcs`) | BSD-3-Clause | https://github.com/whitews/FlowIO |
| FlowUtils | 1.2.2 | Compensation maths used by FlowKit | BSD-3-Clause | https://github.com/whitews/FlowUtils |
| NumPy | 2.4.6 | Arrays, arcsinh | BSD-3-Clause (bundled parts 0BSD, MIT, Zlib, CC0-1.0) | https://numpy.org |
| pandas | 2.3.3 | FlowKit dependency | BSD-3-Clause | https://pandas.pydata.org |
| SciPy | 1.17.1 | FlowKit dependency | BSD-3-Clause | https://scipy.org |
| Bokeh | 3.9.2 | FlowKit dependency (plotting, unused here) | BSD-3-Clause | https://bokeh.org |
| contourpy | 1.3.3 | FlowKit dependency | BSD-3-Clause | https://github.com/contourpy/contourpy |
| lxml | 6.1.3 | FlowKit dependency (GatingML, unused here) | BSD-3-Clause | https://lxml.de |
| networkx | 3.6.1 | FlowKit dependency | BSD-3-Clause | https://networkx.org |
| anytree | 2.13.0 | FlowKit dependency | Apache-2.0 | https://github.com/c0fec0de/anytree |
| psutil | 7.2.2 | FlowKit dependency | BSD-3-Clause | https://github.com/giampaolo/psutil |
| pyarrow | 26.0.0 | pandas/FlowKit dependency | Apache-2.0 | https://arrow.apache.org |
| Jinja2, MarkupSafe | 3.1.6, 3.0.4 | Bokeh dependencies | BSD-3-Clause | https://palletsprojects.com |
| narwhals | 2.28.0 | Bokeh dependency | MIT | https://github.com/narwhals-dev/narwhals |
| pillow | 13.0.0 | Bokeh dependency | MIT-CMU | https://python-pillow.org |
| tornado | 6.6 | Bokeh dependency | Apache-2.0 | https://www.tornadoweb.org |
| xyzservices | 2026.9.1 | Bokeh dependency | BSD-3-Clause | https://github.com/geopandas/xyzservices |
| PyYAML | 6.0.3 | Bokeh dependency | MIT | https://pyyaml.org |
| packaging | 26.3 | dependency | Apache-2.0 OR BSD-2-Clause | https://github.com/pypa/packaging |
| python-dateutil, pytz, tzdata, six | 3.0.0, 2026.5, 2026.6, 1.17.0 | pandas dependencies | Apache-2.0/BSD dual; MIT; Apache-2.0; MIT | PyPI |
| python:3.11.16-slim-bookworm | — | Base image (Debian 12 + CPython 3.11, pip, setuptools, wheel) | PSF-2.0 (Python); Debian packages under their own free licences; pip/wheel MIT; setuptools MIT | https://hub.docker.com/_/python |

Demo data used for testing (not in the image or the repository): Koenig Z, Andorf S,
Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow cytometry dataset. Zenodo 2025,
doi:10.5281/zenodo.15388817, CC BY 4.0.
