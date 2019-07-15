#! /bin/sh
for log in *txt
do
	logbak="$log"".bak"
	cp "$log" "$logbak"
done

for name in 20*fit
do
    newname="$(echo "$name" | cut -c 10-)"
    cp "$name" "$newname"
done
if [ ! -d "./raw" ]; then
  mkdir ./raw
fi
mv 20*fit ./raw/
