DIR=$1
DELAY=$2
SIZE=$3
STATES=$4
TIME=$5
GRAIN=$6

if [[ -n $7 ]]; then
    OUTPUT=$7
else
    mkdir -p rule_gif/
    OUTPUT="rule_gif/temp.gif"
fi


i=0;
for fname in `ls $DIR/*.step | sort -V`; do
    printf "Processing frame: $((i+1)) / $((TIME / GRAIN)) \r";
    tools/viz/step_to_ppm $fname $SIZE $STATES \
        | pamtogif > "$DIR/tmp_$(printf "%05d" $i).gif" 2>/dev/null \
        && i=$((i+1));
done;
echo '\nDone.'

gifsicle -d $DELAY --loop `ls -v $DIR/tmp*.gif` --scale 3 \
         > $OUTPUT
