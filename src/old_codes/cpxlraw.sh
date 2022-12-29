#! /bin/sh
for log in *txt
do
	logbak="$log"".bak"
	cp "$log" "$logbak"
done

for name in 20*fit # new filename standards in 2021
do
    # newname="$(echo "$name" | cut -c 11-)"
    newid="$(echo "$name" | cut -c 9-12)" 
    newname="${newid}.fit"
    cp "$name" "$newname"
done
if [ ! -d "./raw" ]; then
  mkdir ./raw
fi
mv 20*fit ./raw/
