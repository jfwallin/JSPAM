#!/bin/bash
while read GALARG
do
  ./runone.sh $GALARG
done < data.txt
