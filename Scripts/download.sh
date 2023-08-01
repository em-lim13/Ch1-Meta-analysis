#!/bin/bash
here=$(dirname $(realpath $0 ))
destination="${here}/../Data/Rasters"
mkdir -p $destination
echo "Downloading data to $destination"
cat ${here}/urls.txt | xargs -n 1 -P 8 wget -q --directory-prefix=$destination
ls -l $destination
echo "Download complete!"
