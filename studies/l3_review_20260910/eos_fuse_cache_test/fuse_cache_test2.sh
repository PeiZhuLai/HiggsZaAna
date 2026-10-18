#!/bin/bash
# Follow-up of fuse_cache_test.sh (a new name stayed invisible through the mount for the whole 6 min of stat polling):
# two new names, both looked up through the mount before the publish; then every 10 s
#   n1.bin  stat only (how long the cached "no such file" lasts)
#   n2.bin  ls of the directory first, then stat (does a directory listing refresh it?)
# Stops when both are visible or after 20 min. Test files only, removed at the end.
set -u
S=eosuser.cern.ch
X=/eos/user/p/pelai/HZa/sandbox_varsearch_v1/fuse_cache_test2
F=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/fuse_cache_test2
W=$(cd "$(dirname "$0")" && pwd)/work2; mkdir -p $W
head -c 500000 /dev/urandom > $W/B
pub() { xrdcp -f -s --cksum adler32:source $1 root://$S/$2.partial && xrdfs $S mv $2.partial $2; }
xrdfs $S mkdir -p $X; pub $W/B $X/keep.bin; sleep 20
ls $F > /dev/null; stat $F/n1.bin $F/n2.bin > /dev/null 2>&1; echo "primed: n1/n2 looked up, not there"
T0=$(date +%s); pub $W/B $X/n1.bin; pub $W/B $X/n2.bin
echo "published at +$(( $(date +%s) - T0 )) s; server: $(xrdfs $S ls $X | tr '\n' ' ')"
s1=""; s2=""
while [ $(( $(date +%s) - T0 )) -lt 1200 ]; do
  t=$(( $(date +%s) - T0 ))
  [ -z "$s1" ] && stat $F/n1.bin > /dev/null 2>&1 && s1=$t
  [ -z "$s2" ] && { ls $F > /dev/null 2>&1; stat $F/n2.bin > /dev/null 2>&1 && s2=$t; }
  echo "t=+${t}s n1(stat only) ${s1:-missing}  n2(ls dir + stat) ${s2:-missing}"
  [ -n "$s1" ] && [ -n "$s2" ] && break
  sleep 10
done
echo "RESULT n1 visible at ${s1:-never (>1200)} s, n2 visible at ${s2:-never (>1200)} s"
for f in keep n1 n2; do xrdfs $S rm $X/$f.bin; done; xrdfs $S rmdir $X; rm -rf $W
echo DONE
