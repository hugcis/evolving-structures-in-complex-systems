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
./bin/automaton 2d -n $STATES -f "data_2d_$STATES/map/$1.map" || { echo 'Map file not found'; exit 1; }
python plot_2d.py $1 $STATES
open rule_gif/tmp.gif
