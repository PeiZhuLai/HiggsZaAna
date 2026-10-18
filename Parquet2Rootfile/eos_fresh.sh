# Existence checks for files published with eos_safe_write.publish (xrdcp to <final>.partial, then xrdfs mv).
# Source it in a chain script; it replaces the usual fresh() and adds eos_mount_sync().
#
# The EOS FUSE mount caches "no such file": a name it looked up before the publish kept answering "No such file"
# after the file was on the server, for 301 s in one test and more than 362 s in another; listing the directory does
# not refresh it, only the cache running out does (studies/l3_review_20260910/eos_fuse_cache_test/). An overwrite
# showed the new size and content within a second. So the server answers whether a file exists, and since the
# downstream steps read through /eos, the checks then wait for the mount to show the same size.
#
#   fresh PATH T        PATH is non-empty and modified at or after T (epoch s). /eos paths: xrdfs stat on the server,
#                       then wait up to EOS_FUSE_WAIT s (default 1200) for the mount. Other paths: stat as before.
#   eos_mount_sync DIR  every file under DIR on the server is visible through the mount with the server size
#                       (wait up to EOS_FUSE_WAIT s). Run it before a step that reads DIR through /eos and skips
#                       missing files (apply_bdt_data.py only warns about a missing era).
# Reasons for a failure go to stderr.

_eos_xrd() {   # /eos path -> "server path"; /eos/home-<l>/<user>/ is /eos/user/<l>/<user>/ on eosuser
  case $1 in
    /eos/home-?/*) local r=${1#/eos/home-}; echo "eosuser.cern.ch /eos/user/${r:0:1}/${r:2}" ;;
    /eos/user/*) echo "eosuser.cern.ch $1" ;;
    /eos/cms/*) echo "eoscms.cern.ch $1" ;;
    /eos/project/?/*) echo "eosproject-${1:13:1}.cern.ch $1" ;;
    *) return 1 ;;
  esac
}

_eos_mount_wait() {   # PATH SIZE: wait until the mount shows PATH with SIZE bytes
  local w=0
  until [ "$(command stat -c %s "$1" 2>/dev/null)" = "$2" ]; do
    if [ $w -ge "${EOS_FUSE_WAIT:-1200}" ]; then
      echo "  $1: on the server ($2 B) but the /eos mount does not show it after $w s" >&2; return 1
    fi
    sleep 15; w=$((w + 15))
  done
}

fresh() {
  local f=$1 t=$2 sp st size mt
  if ! sp=$(_eos_xrd "$f"); then
    [ -s "$f" ] && [ "$(command stat -c %Y "$f")" -ge "$t" ]; return
  fi
  st=$(xrdfs ${sp%% *} stat "${sp#* }" 2>/dev/null) || { echo "  $f: not on the server" >&2; return 1; }
  size=$(awk '/^Size:/ {print $2}' <<< "$st")
  mt=$(date -u -d "$(awk '/^MTime:/ {print $2, $3}' <<< "$st")" +%s)   # xrdfs prints UTC
  [ "${size:-0}" -gt 0 ] || { echo "  $f: empty on the server" >&2; return 1; }
  [ "$mt" -ge "$t" ] || { echo "  $f: on the server but older than the start of this step" >&2; return 1; }
  _eos_mount_wait "$f" "$size"
}

eos_mount_sync() {
  local d=${1%/} sp srv base n=0 bad=0 size p
  sp=$(_eos_xrd "$d") || { echo "  eos_mount_sync: $d is not an /eos path" >&2; return 1; }
  srv=${sp%% *}; base=${sp#* }
  while read -r size p; do
    n=$((n + 1))
    _eos_mount_wait "$d${p#$base}" "$size" || bad=$((bad + 1))
  done < <(xrdfs $srv ls -l -R "$base" | awk '$1 ~ /^-/ {print $4, $NF}')   # mode owner group size date time path
  [ $n -gt 0 ] || { echo "  eos_mount_sync: no files under $d on the server" >&2; return 1; }
  [ $bad -eq 0 ] || { echo "  eos_mount_sync: $bad of $n files under $d not visible through the mount" >&2; return 1; }
}
