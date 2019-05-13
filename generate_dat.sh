for rule in {1..255}
do
    echo $rule;
    echo "" > data/out$rule\_pag.log
    for stepname in $( ls steps/out$rule\_* ); do
        ./paq8l archout $stepname >> data/out$rule\_paq.log;
    done
done

rm archout.paq8l
python log_to_dat.py
