#!/bin/bash
# Arguments <from-no>
set -o nounset
set -o errexit

#dart -c ./grid_test.dart 5 $1  # needs fine grid


form=${1:-1}
repeats=${2:-10}
distance=${3:-25.0}


#for i in $(seq 0 ${repeats}) ; do
#    echo "============================================================"
#    echo "$i"
#    echo "============================================================"
#    dart -c ./grid_test.dart $i ${form} ${distance}
#done

echo "Test Parameters"
echo "Form:     ${form}"
echo "Repeats:  ${repeats}"
echo "Distance: ${distance}"
dart -c ./grid_test.dart ${form} ${repeats} ${distance}
