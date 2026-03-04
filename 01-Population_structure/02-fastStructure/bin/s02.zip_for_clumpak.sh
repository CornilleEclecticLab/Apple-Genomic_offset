#!/bin/bash

cd ./out/dom_admixed_rm/s01.faststructure
for i in $(seq 2 10)
do
    zip -q K${i}.zip *.${i}.meanQ
done

zip -q ../input_clumpak.zip K*.zip
rm -rf K*.zip

