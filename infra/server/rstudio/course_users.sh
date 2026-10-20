#!/usr/bin/with-contenv bash
# shellcheck shell=bash
# Container start-up (s6 cont-init): create the course accounts, the shared
# group work directories, and give each account access to goose through the
# LLM proxy.
#
# /run/course/users.tsv (root-only) has one account per line:
#   username <TAB> uid <TAB> password <TAB> llm-proxy key
#
# /home and /group_work live on a persistent volume while /etc/passwd does
# not, so this runs on every start and must leave existing files intact.
#
# Unlike rocker's batch_user_creation.sh, accounts are NOT added to the staff
# group: staff can write to the shared R library, which would let one
# participant break (or plant code in) packages everyone else loads.
set -euo pipefail

users_file=/run/course/users.tsv
if [ ! -r "$users_file" ]; then
    echo "course_users: no $users_file, skipping"
    exit 0
fi

# rocker's default account is in staff and has a (random) password; nobody
# should log in with it on the shared server
usermod --lock "${DEFAULT_USER:-rstudio}"

# every course account is in this group, which owns /group_work
getent group course >/dev/null || groupadd --gid 3000 course

# Shared goose settings. The per-account key goes in each home (below).
mkdir -p /etc/goose
cat >/etc/goose/config.yaml <<EOF
GOOSE_PROVIDER: litellm
GOOSE_MODEL: ${COURSE_LLM_MODEL}
LITELLM_HOST: ${COURSE_LLM_HOST}
GOOSE_DISABLE_KEYRING: true
GOOSE_MODE: smart_approve
GOOSE_TELEMETRY_ENABLED: false
GOOSE_CLI_THEME: light
extensions:
  developer:
    enabled: true
    type: builtin
    name: developer
    bundled: true
    timeout: 300
EOF
chmod 644 /etc/goose/config.yaml

while IFS=$'\t' read -r username uid password llm_key; do
    [ -n "$username" ] || continue
    home="/home/$username"

    if ! id -u "$username" >/dev/null 2>&1; then
        useradd --uid "$uid" --user-group --groups course --shell /bin/bash \
            --home-dir "$home" --no-create-home "$username"
    fi
    echo "$username:$password" | chpasswd

    if [ ! -d "$home" ]; then
        cp -rT /etc/skel "$home"
        mkdir -p "$home/.config/rstudio"
        cat >"$home/.config/rstudio/rstudio-prefs.json" <<EOF
{
    "save_workspace": "never",
    "load_workspace": false,
    "always_save_history": false
}
EOF
        chown -R "$username:$username" "$home"
    elif [ "$(stat -c %u "$home")" != "$uid" ]; then
        # the account list changed between runs and the uid moved
        chown -R "$username:$username" "$home"
    fi
    chmod 700 "$home"

    # shortcut to the shared directories, visible in the Files pane
    if [ ! -e "$home/group_work" ] && [ ! -L "$home/group_work" ]; then
        ln -s /group_work "$home/group_work"
        chown -h "$username:$username" "$home/group_work"
    fi

    install -d -o "$username" -g "$username" -m 700 "$home/.config" "$home/.config/goose"
    printf 'LITELLM_API_KEY: %s\n' "$llm_key" >"$home/.config/goose/secrets.yaml"
    chown "$username:$username" "$home/.config/goose/secrets.yaml"
    chmod 600 "$home/.config/goose/secrets.yaml"

    # default instructions for goose; written once, so participants can edit them
    hints="$home/.config/goose/.goosehints"
    if [ ! -e "$hints" ]; then
        cat >"$hints" <<'EOF'
You are helping a bioinformatics course participant working on a shared Rstudio Server. 
They are interacting with you through the terminal. 

- If you provide a long answer, write it to a seperate .md file, so the participant
  can read it in the RStudio IDE. Do this in particular if your answer includes tables, 
  plots, or long code blocks.
- Write analysis code to .R or .qmd files in the project, rather than only
  showing it in your reply, so the participant can read and run it.
- All R packages for the course are installed. Do not install or update
  packages unless the participant asks for it.
EOF
        chown "$username:$username" "$hints"
    fi

    echo "course_users: $username ready"
done <"$users_file"

# Shared directories for the project groups, open to every course account.
# setgid keeps new files in the course group; the default ACL makes them
# group-writable whatever the umask of the account that creates them.
mkdir -p /group_work
chgrp course /group_work
chmod 2770 /group_work
setfacl --default --modify "u::rwx,g::rwx,o::---" /group_work
for i in $(seq 1 "${COURSE_GROUP_COUNT:-0}"); do
    # created inside /group_work, so it inherits the group, setgid and ACL
    mkdir -p "/group_work/group$i"
done
echo "course_users: /group_work ready"
