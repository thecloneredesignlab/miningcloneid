#!/usr/bin/env bash
set -euo pipefail

SSH_USER=4482173@red.moffitt.org
SLURM_BIN=/cm/shared/apps/slurm/current/bin
export SSH_USER SLURM_BIN

remote() {
  ssh "$SSH_USER" "bash -lc 'export PATH=$SLURM_BIN:\$PATH; $1'"
}

jobs=(17434072 17434073 17434074 17434075 17434076 17434077 17434078 17434079 17434080 17434081 17434082)

echo "== QUEUE =="
remote "squeue -h -u 4482173 -o '%i|%T|%R' | egrep '^(17434072|17434073|17434074|17434075|17434076|17434077|17434078|17434079|17434080|17434081|17434082)(_|\\|)' || true"

echo
echo "== JOB SNAPSHOT 17434072 =="
remote "squeue -r -h -j 17434072 -o '%i|%T|%R|%M' | awk -F'|' 'BEGIN{r=0;p=0} \$2==\"RUNNING\"{r++} \$2==\"PENDING\"{p++} END{print \"RUNNING=\" r \", PENDING=\" p}'"
remote "sacct -n -P -j 17434072 --format=JobID,State | egrep '^17434072_[0-9]+\\|' | cut -d'|' -f2 | sort | uniq -c || true"
remote "sacct -n -P -j 17434072 --format=JobID,State | egrep '^17434072_[0-9]+\\|COMPLETED$' | cut -d'|' -f1 | wc -l"

echo
echo "== COMPLETED TASKS AT ITERMAX =="
remote "ids=\$(sacct -n -P -j 17434072 --format=JobID,State | egrep '^17434072_[0-9]+\\|COMPLETED$' | cut -d'|' -f1 | sed 's/^17434072_//'); count=0; total=0; for id in \$ids; do total=\$((total+1)); f=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/log/o2mw_v01_i01_17434072_\${id}.out; line=\$(tail -n 200 \"\$f\" | grep 'Iteration:' | tail -1 || true); n=\$(printf '%s\\n' \"\$line\" | sed -n 's/.*Iteration:[[:space:]]*\\([0-9][0-9]*\\).*/\\1/p'); if [ -n \"\$n\" ] && [ \"\$n\" -ge 500 ]; then count=\$((count+1)); fi; done; printf 'completed=%s itermax=%s reached=%s\\n' \"\$total\" 500 \"\$count\""

echo
for jid in "${jobs[@]:1}"; do
  echo "== JOB ${jid} =="
  remote "squeue -h -j ${jid} -o '%i|%T|%R|%M' || true"
  remote "sacct -n -P -j ${jid} --format=JobID,State | head -20 || true"
  echo
done
