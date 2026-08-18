# A container definition for the container used for the window
# based pruning method method in PopGLen >=0.5.0. A minimal R
# environment with dplyr, data.table, and R.utils (to support
# gzip input for fread()).

FROM ghcr.io/r-hub/r-minimal/r-minimal:4.6.1

RUN installr -a "libgomp zlib" -t "openmp-dev zlib-dev" -d data.table dplyr R.utils
