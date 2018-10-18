#!/bin/bash

opts="-avhun"
info="DRY RUN! use doit as first arg to actually sync"
[ "$1" == "doit" ] && opts="-avhu" && info="not a dry run..."
echo ${info}
rsync ${opts} --progress --exclude-from=rsync_ignore.txt /Volumes/one/run4 .
