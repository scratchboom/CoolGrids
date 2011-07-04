#!/bin/sh -e
for h in cluster17 cluster11 cluster05; do
    ssh $h 'ps aux | grep rf-plasma.release | grep -v grep' | sed "s/^/$h: /"
done
