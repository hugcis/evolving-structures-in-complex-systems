#!/bin/bash

for i in maps/train/*.map; do
    bin/automaton 2d -n 3 -f $i;
done;

for i in maps/test/*.map; do
    bin/automaton 2d -n 3 -f $i;
done;
