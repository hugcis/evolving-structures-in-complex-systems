# Usage: `./generate_frames.sh [rule] [n_states (default 3)]`
#
# This script will run a simulation of a cellular automaton specified by the
# rule hash and number of states and generate a gif file of its evolution.

STATES=3;
TIME=1000;
SIZE=256;
GRAIN=50;
DELAY=30;
OUTPUT="";
Q="";
PAT="";

OPTIND=1  # Reset in case getopts has been used previously in the shell.

show_help() {
    echo "\
Visualization tool for generating cellular automata GIF evolution.\n\n\
Usage:\n\
    tools/viz/generate_frames.sh [options] [rule_number|file_location]\n\
\n\
Options:\n\
    -h                 Print this message.\n\
    -t <time>          Number of time steps of the simulation \n\
                       (default: 1000).\n\
    -s <size>          Size of the automaton being simulated (default: 256).\n\
    -g <grain>         Grain of the final GIF image - one image every\n\
                       grain steps (default: 50).\n\
    -d <delay>         GIF delay parameter (default: 30).\n\
    -n <states>        Number of states of the automaton (default: 3).\n\
    -o <out_file>      Output GIF file path (default: ./rule_gif/temp.gif).\n\
    -j <pattern_file>  Pattern to load for simulation.\n"
}

while getopts "h?qd:g:s:t:n:j:o:" opt; do
    case "$opt" in
        h|\?)
            show_help
            exit 0
            ;;
        t)  TIME=$OPTARG
            ;;
        s)  SIZE=$OPTARG
            ;;
        g)  GRAIN=$OPTARG
            ;;
        d)  DELAY=$OPTARG
            ;;
        q)  Q=--masking
            ;;
        n)  STATES=$OPTARG
            ;;
        o)  OUTPUT=$OPTARG
            ;;
        j)  PAT="--pattern $OPTARG"
            ;;
    esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

tmpdir=$(mktemp -d)

if echo $@ |grep --quiet -e '[0-9]\+\b'; then
    INPUT="$PWD/data_2d_$STATES/map/$@.map"
else
    INPUT=$@
fi

./bin/automaton 2d --n_states $STATES\
                --temp_output\
                --input_file $INPUT\
                --timesteps $TIME\
                --size $SIZE --write_grain $GRAIN\
                $Q --early_stopping --output $tmpdir $PAT ||
    { echo 'Failure during automaton simulation'; exit 1; }

tools/viz/make_frames.sh $tmpdir $DELAY $SIZE $STATES $TIME $GRAIN $OUTPUT

rm $tmpdir/tmp_*.step
rm -R $tmpdir
