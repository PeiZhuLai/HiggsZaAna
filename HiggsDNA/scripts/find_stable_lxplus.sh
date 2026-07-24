#!/bin/bash

results=$(for i in $(seq 900 999); do
  # 排除 901–905，避免掃到 GUP/GPU/特殊用途節點
  if [ "$i" -ge 901 ] && [ "$i" -le 905 ]; then
    continue
  fi

  host=lxplus${i}.cern.ch

  timeout 5s ssh \
    -o BatchMode=yes \
    -o ConnectTimeout=3 \
    -o ConnectionAttempts=1 \
    -o StrictHostKeyChecking=accept-new \
    "$host" \
    "awk '{printf \"%d %.1f\n\", int(\$1), \$1/86400}' /proc/uptime | xargs -I{} echo {} $host" \
    2>/dev/null
done | grep -E '^[0-9]' | sort -n)

# 印出完整清單（uptime 由低到高：越上面代表越剛重開機）
echo "=== uptime seconds | days | host（由低到高）==="
echo "$results" | head -20

# 建議節點：挑 uptime 介於 1h(3600s)–2 天(172800s) 之間、最低的那台
# （避開 <1h 可能還在維護窗、也避開太舊快被重開機的）
echo
echo "=== 建議 tmux 節點 ==="
echo "$results" | awk '$1 >= 3600 && $1 <= 172800 {print; found=1; exit} END {if (!found) print "（清單中無 1h–2 天的節點，改用最上方幾台）"}'
