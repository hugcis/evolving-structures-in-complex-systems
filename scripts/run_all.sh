#!/bin/bash

echo -n "Retrieving rule files...      "
curl -L https://drive.google.com/uc\?export\=download\&id\=1fymRRN-Yeig560CkXrLTfpl879YLP_UF > maps.zip
echo "Done."

echo -n "Unzipping...     "
unzip maps > /dev/null
echo "Done."

echo "Processing rules"
./scripts/generate_results_from_maps.sh
echo "Done."
