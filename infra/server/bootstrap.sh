#!/bin/bash
# First-boot setup, run by cloud-init: mount the persistent data volume, start
# the stack and create the per-account LLM keys. Safe to re-run by hand.
# Log: /var/log/course-bootstrap.log
set -euo pipefail
exec > >(tee -a /var/log/course-bootstrap.log) 2>&1

cd /opt/course
set -a; source .env; set +a

# --- persistent data volume --------------------------------------------------
# Terraform attaches it only after the instance exists, so wait for it.
device="/dev/disk/by-id/nvme-Amazon_Elastic_Block_Store_${DATA_VOLUME_ID//-/}"
for _ in $(seq 1 60); do
    [ -e "$device" ] && break
    echo "waiting for data volume $DATA_VOLUME_ID ..."
    sleep 5
done
[ -e "$device" ] || { echo "data volume $DATA_VOLUME_ID never appeared" >&2; exit 1; }

blkid "$device" >/dev/null || mkfs.ext4 -L course-data "$device"
mkdir -p /srv/course
grep -q 'LABEL=course-data' /etc/fstab ||
    echo 'LABEL=course-data /srv/course ext4 defaults,nofail 0 2' >>/etc/fstab
mountpoint -q /srv/course || mount /srv/course
mkdir -p /srv/course/home /srv/course/group_work /srv/course/postgres /srv/course/caddy

# --- stack -------------------------------------------------------------------
# Leave ~12 GiB for the host, proxy and database; when RStudio sessions exceed
# the rest, the kernel kills the largest R session instead of the server.
mem_gib=$(awk '/MemTotal/ { print int($2 / 1024 / 1024) }' /proc/meminfo)
sed -i '/^RSTUDIO_MEM_LIMIT=/d' .env
echo "RSTUDIO_MEM_LIMIT=$((mem_gib - 12))g" >>.env

docker compose up -d --build --wait

docker compose exec -T llm-proxy python3 /run/course/create_llm_keys.py

echo "bootstrap finished: https://${SITE_ADDRESS}"
