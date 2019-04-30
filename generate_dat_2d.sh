list="254960 254961 254962 254963 254964 254965 254966 254967 254968 254969 254970 254971 254972 254973 254974 254975"

for fname in $( ls data_2d/*.dat) ; do
    if [[ !($fname == *"_paq"*)  ]]; then
        RULE=$(echo $fname |sed -E 's/data_2d\/out([[:digit:]]+)\.dat/\1/g');
        echo $RULE;
        if [[ $list =~ (^|[[:space:]])$RULE($|[[:space:]]) ]]; then
            echo "" > data_2d/out$RULE\_paq.log;
            for stepname in $( ls step_2d/out$RULE\_* ); do
                ./paq8l archout $stepname >> data_2d/out$RULE\_paq.log;
            done
        fi
    fi
done
