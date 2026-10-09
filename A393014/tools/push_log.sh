#!/bin/bash
# Copy this machine's search log to another machine (results/remote_<name>/search_log.jsonl there),
# so its campaigns skip lengths already done here.   usage: push_log.sh user@host [name] [path]
host=$1; name=${2:-$(hostname -s | tr '.-' '__')}; path=${3:-dbpal}
cd "$(dirname "$0")/.."
ssh -o BatchMode=yes "$host" "mkdir -p $path/results/remote_$name"
rsync -az results/search_log.jsonl "$host:$path/results/remote_$name/search_log.jsonl"
echo "pushed $(wc -l < results/search_log.jsonl) log lines to $host:$path/results/remote_$name/"
