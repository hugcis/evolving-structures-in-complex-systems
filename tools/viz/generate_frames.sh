# Usage: `./generate_frames.sh [rule] [n_states (default 3)]`
#
# This script will run a simulation of a cellular automaton specified by the
# rule hash and number of states and generate a gif file of its evolution.
# TODO: add support for options in both the bin/automaton executable and the
#       python script.

STATES=3;
TIME=1000;
SIZE=256;
GRAIN=50;
DELAY=30;
Q="";
PAT="";

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?qd:g:s:t:n:j:" opt; do
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

tools/viz/make_frames.sh $tmpdir $DELAY $SIZE $STATES $TIME $GRAIN

rm $tmpdir/tmp_*.step
rm -R $tmpdir
