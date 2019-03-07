#!/bin/bash

find ../run3 -name "*.root" | grep -v subjets_ca_sjR10 | cut -d / -f 3- | tee exclude.txt
rsync -av --exclude-from=exclude.txt ../run3/* .
