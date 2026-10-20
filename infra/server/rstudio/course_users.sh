#!/usr/bin/with-contenv bash
# shellcheck shell=bash
# Container start-up (s6 cont-init): create the course accounts and give each
# one access to goose through the LLM proxy.
#
# /run/course/users.tsv (root-only) has one account per line:
#   username <TAB> uid <TAB> password <TAB> llm-proxy key
#
# /home lives on a persistent volume while /etc/passwd does not, so this runs
# on every start and must leave existing home directories intact.
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

# Shared goose settings. The per-account key goes in each home (below).
mkdir -p /etc/goose
cat >/etc/goose/config.yaml <<EOF
GOOSE_PROVIDER: litellm
GOOSE_MODEL: ${COURSE_LLM_MODEL}
LITELLM_HOST: ${COURSE_LLM_HOST}
GOOSE_DISABLE_KEYRING: true
GOOSE_MODE: smart_approve
GOOSE_TELEMETRY_ENABLED: false
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
        useradd --uid "$uid" --user-group --shell /bin/bash --home-dir "$home" --no-create-home "$username"
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

    install -d -o "$username" -g "$username" -m 700 "$home/.config" "$home/.config/goose"
    printf 'LITELLM_API_KEY: %s\n' "$llm_key" >"$home/.config/goose/secrets.yaml"
    chown "$username:$username" "$home/.config/goose/secrets.yaml"
    chmod 600 "$home/.config/goose/secrets.yaml"

    echo "course_users: $username ready"
done <"$users_file"
