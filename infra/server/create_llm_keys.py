"""Create (or update) one LiteLLM key per course account.

Runs inside the llm-proxy container (see bootstrap.sh), so the master key never
leaves it. Keys are chosen by Terraform, which lets the RStudio container write
them into each account without talking to the proxy. Idempotent.

With the argument "report", prints the spend per account instead.
"""

import json
import os
import sys
import urllib.error
import urllib.request

BASE = "http://localhost:4000"
HEADERS = {
    "Authorization": "Bearer " + os.environ["LITELLM_MASTER_KEY"],
    "Content-Type": "application/json",
}


def call(method, path, body=None):
    data = json.dumps(body).encode() if body is not None else None
    request = urllib.request.Request(BASE + path, data=data, headers=HEADERS, method=method)
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def key_exists(key):
    try:
        call("GET", "/key/info?key=" + key)
        return True
    except urllib.error.HTTPError as error:
        if error.code in (400, 404):
            return False
        raise


with open("/run/course/llm-keys.json") as f:
    config = json.load(f)

if sys.argv[1:] == ["report"]:
    total = 0.0
    for account in config["users"]:
        info = call("GET", "/key/info?key=" + account["key"])["info"]
        spend = info.get("spend") or 0.0
        total += spend
        print(f"{account['user']:<16} ${spend:8.2f} of ${info.get('max_budget') or 0:.2f}")
    print(f"{'total':<16} ${total:8.2f}")
    sys.exit()

for account in config["users"]:
    settings = {
        "key": account["key"],
        "models": [config["model"]],
        "max_budget": config["max_budget"],
        "rpm_limit": config["rpm_limit"],
        "max_parallel_requests": 4,
    }
    if key_exists(account["key"]):
        call("POST", "/key/update", settings)
        print(f"updated key for {account['user']}")
    else:
        call("POST", "/key/generate", {**settings, "key_alias": account["user"], "user_id": account["user"]})
        print(f"created key for {account['user']}")
