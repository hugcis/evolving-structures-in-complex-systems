#!/bin/bash

echo -n "Retrieving rule files...      "
curl https://doc-10-c4-docs.googleusercontent.com/docs/securesc/ha0ro937gcuc7l7deffksulhg5h7mbp1/d21mvss4rt7sl74ise8v8jm5i9kou7kp/1573120800000/15582562944099467373/\*/1fymRRN-Yeig560CkXrLTfpl879YLP_UF\?e\=download > maps.zip
echo "Done."

echo -n "Unzipping...     "
unzip maps > /dev/null
echo "Done."

echo "Processing rules"
./scripts/generate_results_from_maps.sh
echo "Done."
