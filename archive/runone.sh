#!/bin/sh
mkdir -p output/$4

java -classpath build/jar/jspamarchive.jar edu.gmu.cds.sim.trim.Driver $1 $2 $3 output/$4

echo "set title '$4'" > output/$4/plot
echo "set xrange [-$3:$3]" >> output/$4/plot
echo "set yrange [-$3:$3]" >> output/$4/plot
echo "set size square" >> output/$4/plot
echo "set pointsize 0.2" >> output/$4/plot

cd output/$4
ls *.dat > files.txt
cd ../../


while read STP
do
  echo "plot '$STP' using 2:4 notitle" >> output/$4/plot
done < output/$4/files.txt

echo "pause 10" >> output/$4/plot
