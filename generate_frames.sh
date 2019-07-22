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


if [ -z $3 ];
then TIME=1000;
else TIME=$3;
fi

if [ -z $4 ];
then SIZE=256;
else SIZE=$4;
fi

if [ -z $5 ];
then GRAIN=50;
else GRAIN=$5;
fi

if [ -z $6 ];
then DELAY=30;
else DELAY=$6;
fi

tmpdir=$(mktemp -d)

./bin/automaton 2d -n $STATES -m -f "data_2d_$STATES/map/$1.map" -t $TIME \
                -s $SIZE -w $GRAIN -e -o $tmpdir ||
    { echo 'Map file not found'; exit 1; }

i=0;
for fname in `ls $tmpdir/*.step | sort -V`; do
    printf "Processing frame: $((i+1)) / $((TIME / GRAIN)) \r";
    ./step_to_ppm $fname $SIZE $STATES \
        | pamtogif > "$tmpdir/tmp_$(printf "%05d" $i).gif" 2>/dev/null \
        && i=$((i+1));
done;
echo '\nDone.'
gifsicle -d $DELAY --loop `ls -v $tmpdir/tmp*.gif` --scale 3 \
         > rule_gif/temp.gif

rm $tmpdir/tmp_*.step
rm -R $tmpdir

# open rule_gif/temp.gif
