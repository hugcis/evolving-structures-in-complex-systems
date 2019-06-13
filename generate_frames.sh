# Usage: `./generate_frames.sh [rule] [n_states (default 3)]`
#
# This script will run a simulation of a cellular automaton specified by the
# rule hash and number of states and generate a gif file of its evolution.
# TODO: add support for options in both the bin/automaton executable and the
#       python script.

if [ -z $2 ];
then STATES=3;
else STATES=$2;
fi


./bin/automaton 2d -n $STATES -m -f "data_2d_$STATES/map/$1.map" ||
    { echo 'Map file not found'; exit 1; }
i=0;
for fname in `ls rule_gif/*.step | sort -V`; do
    ./test $fname 256 $i;
    pamtogif rule_gif/tmp_$(printf "%05d" $i).ppm > rule_gif/tmp_$(printf "%05d" $i).gif &&
    i=$((i+1));
done;

gifsicle `ls -v rule_gif/tmp*.gif` > rule_gif/temp.gif

rm rule_gif/tmp_*.step
rm rule_gif/tmp_*.ppm
rm rule_gif/tmp_*.gif

open rule_gif/temp.gif
