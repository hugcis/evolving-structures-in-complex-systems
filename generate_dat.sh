for rule in {1..255}
do
    echo $rule;
    ./paq8l archout steps/out$rule\_* > data/out$rule\_paq.log;
done

rm archout.paq8l
python log_to_dat.py
