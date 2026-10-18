#!/bin/bash
# How long the EOS FUSE mount shows the old state after a publish done like eos_safe_write.publish
# (xrdcp to <final>.partial, then xrdfs mv): (a) overwrite of a file the mount has already read,
# (b) a new name the mount has already looked up and found missing. Polls the mount every 5 s for 6 min.
# Test files only, under sandbox_varsearch_v1/fuse_cache_test/ (removed at the end).
set -u
S=eosuser.cern.ch
X=/eos/user/p/pelai/HZa/sandbox_varsearch_v1/fuse_cache_test
F=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/fuse_cache_test
W=$(cd "$(dirname "$0")" && pwd)/work; mkdir -p $W
head -c 1000000 /dev/urandom > $W/A; head -c 1300000 /dev/urandom > $W/B
mA=$(md5sum < $W/A | cut -c1-32); mB=$(md5sum < $W/B | cut -c1-32)
pub() { xrdcp -f -s --cksum adler32:source $1 root://$S/$2.partial && xrdfs $S mv $2.partial $2; }
xrdfs $S mkdir -p $X
pub $W/A $X/over.bin
sleep 20
echo "prime: over.bin via mount md5 $(md5sum < $F/over.bin | cut -c1-32) (A=$mA)  new.bin: $(ls $F/new.bin 2>&1 | tail -1)"
md5sum < $F/over.bin > /dev/null; stat -c %s $F/new.bin 2>/dev/null
T0=$(date +%s)
pub $W/B $X/over.bin; pub $W/B $X/new.bin
echo "published B (md5 $mB, 1300000 B) at +$(( $(date +%s) - T0 )) s; server: $(xrdfs $S stat $X/over.bin | grep Size) / $(xrdfs $S stat $X/new.bin | grep Size)"
for i in $(seq 0 72); do
  t=$(( $(date +%s) - T0 ))
  o_s=$(stat -c %s $F/over.bin 2>/dev/null || echo NA); o_m=$(md5sum < $F/over.bin 2>/dev/null | cut -c1-32)
  n_s=$(stat -c %s $F/new.bin 2>/dev/null || echo NA); n_m=$(md5sum < $F/new.bin 2>/dev/null | cut -c1-32)
  lab() { [ "$1" = "$mB" ] && echo new || { [ "$1" = "$mA" ] && echo OLD || echo "${1:-none}"; }; }
  echo "t=+${t}s over: size $o_s content $(lab $o_m) | new: size $n_s content $(lab $n_m)"
  sleep 5
done
xrdfs $S rm $X/over.bin; xrdfs $S rm $X/new.bin; xrdfs $S rmdir $X; rm -rf $W
echo DONE
