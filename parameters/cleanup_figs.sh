#!/bin/bash

for image_file in $(ls figs/)
do
if grep $image_file *.log -c > 1
then
    echo "File $image_file is in use."
    cp figs/$image_file ../submit_arxiv/figs
fi
done
