#!/bin/sh

# on normal bash `seq` has the -f param
# for i in $(seq -f "%02g" 1 21)

# on alpine version, seq uses -w
for i in $(seq -w 1 21)
do
if test -f "pdfs/BP_$i.pdf"; then
    echo "BP_$i.pdf exists"
    cp zero.png "BP_$i.png"
else
    echo "BP_$i.pdf does not exist"
    cp failed.png "BP_$i.png"
fi
done