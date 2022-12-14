#!/bin/bash

echo "sample,prop_susp" > prop_susp.csv && 
    for f in .tmp/*.csv; do tail -n 1 $f; done >> prop_susp.csv &&
    rm -r .tmp
