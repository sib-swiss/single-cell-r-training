# Course server on AWS

One EC2 instance running a shared RStudio Server for all participants, each with
their own account, plus the [goose](https://github.com/aaif-goose/goose) coding
agent wired to an LLM through a proxy that holds the real API key.

```
internet ──443── Caddy (HTTPS) ── RStudio Server (course image + goose, one Linux account per participant)
                                        │
                                        └── LiteLLM proxy ── OpenRouter
                                              └── Postgres (keys, budgets, spend)
```

All four run as Docker containers on the instance (`server/compose.yaml`). The
RStudio container is the public course image (`Docker/Dockerfile`) plus a thin
layer (`server/rstudio/`) that adds goose and creates the accounts.

## How the API key stays hidden

- The OpenRouter key exists only in the LiteLLM container. The RStudio
  container never sees it.
- Each account gets its own LiteLLM key in `~/.config/goose/secrets.yaml`. Participants
  can see that key, but it only works against the proxy inside the server, only for
  the course model, and only up to `llm_budget_per_user_usd`.
- The secrets reach the instance through EC2 user data. IMDSv2 with a hop limit
  of 1 makes the metadata service, and therefore the user data and the instance
  role, unreachable from inside containers.
- Accounts are not in the `staff` group, so they cannot write to the shared R
  library.

All secrets are also in the Terraform state. Keep the state local (or in an
encrypted backend) and `terraform destroy` after the course.

## Before the first run

- Terraform ≥ 1.6, the AWS CLI with credentials for the target account, and the
  [Session Manager plugin](https://docs.aws.amazon.com/systems-manager/latest/userguide/session-manager-working-with-install-plugin.html)
  for admin access to the server (there is no SSH port).
- **vCPU quota**: `r7i.16xlarge` needs 64 vCPUs under "Running On-Demand Standard
  (A, C, D, H, I, M, R, T, Z) instances" in eu-central-1. Check Service Quotas
  well ahead; increases can take days.
- **OpenRouter key**: on openrouter.ai, create a key for this course only, with a
  credit limit of roughly `participants × llm_budget_per_user_usd` plus a margin.

## Running a course

```bash
cd infra
cp terraform.tfvars.example terraform.tfvars   # fill in
terraform init
terraform apply
```

`terraform apply` returns after a minute or two. The server then needs
roughly 10-15 minutes to install Docker, pull the course image and start. Follow
along with:

```bash
aws ssm start-session --target "$(terraform output -raw instance_id)"
sudo tail -f /var/log/course-bootstrap.log    # ends with "bootstrap finished: https://..."
```

The URL is `terraform output url` (`https://<ip-with-dashes>.sslip.io`).
Accounts and passwords are in `infra/credentials.csv`.

### Smoke test before participants arrive

1. Log in as `user01`. Open the Terminal tab and run `goose`, then ask it something
   that needs a tool, e.g. *"list the files in my home directory"*.
2. Log in as `user02` in a private window. Check that `ls /home/user01` is
   denied, and that `cat /run/course/users.tsv` is denied.
3. Check that `curl -s -m 5 http://169.254.169.254/` from the RStudio terminal
   times out or fails.
4. As `user01`, save a file in `~/group_work/group1`; as `user02`, edit and save
   the same file.
5. Run part of the Seurat workflow in a few accounts at once and watch `htop`
   on the host.

### During the course

On the host (via `aws ssm start-session`), in `/opt/course`:

```bash
sudo docker compose exec llm-proxy python3 /run/course/create_llm_keys.py report   # LLM spend per account
sudo docker compose logs -f rstudio                                                # account creation, RStudio
sudo docker compose restart rstudio                                                # logs everyone out; homes are kept
```

Tokens and requests per account, from the proxy's request log:

```bash
sudo docker compose exec db psql -U litellm -d litellm -c '
SELECT t.key_alias                          AS account,
       count(*)                             AS requests,
       sum(s.prompt_tokens)                 AS input_tokens,
       sum(s.completion_tokens)             AS output_tokens,
       round(sum(s.spend)::numeric, 2)      AS usd
FROM "LiteLLM_SpendLogs" s
JOIN "LiteLLM_VerificationToken" t ON t.token = s.api_key
GROUP BY t.key_alias
ORDER BY usd DESC;'
```

For a breakdown per day, add `date_trunc('day', s."startTime") AS day` to the
`SELECT` and use `GROUP BY t.key_alias, day ORDER BY day, usd DESC`. The log
records tokens, cost, model and time for each request, but not the prompts or
responses.

To raise a budget, change `llm_budget_per_user_usd`. The account list and the
budget are part of the instance's user data, so `terraform apply` will **replace
the instance**. Home directories, the LLM spend so far and the URL survive that
(they live on the separate data volume and the Elastic IP), but everyone is
disconnected for ~15 minutes. That is why `participant_count` includes spares.

### After the course

```bash
terraform destroy
```

This deletes the data volume with all home directories too. Also revoke the
OpenRouter key.

## Things to know

- **Memory is shared.** Open-source RStudio Server has no per-user limits. The
  RStudio container may use all but 12 GiB of the instance's memory. When it
  runs out, the kernel kills the largest process, which is normally the R
  session that caused it. Size the instance at ~12-16 GiB per participant.
- **sslip.io** gives a hostname without owning a domain. If its Let's Encrypt
  rate limit is ever exhausted, Caddy falls back to ZeroSSL automatically.
- **goose** runs in `smart_approve` mode, so it asks before risky actions, and
  with the light CLI theme, which matches RStudio's default.
  The shared settings are in `/etc/goose/config.yaml` inside the container. A
  participant's own `~/.config/goose/config.yaml` overrides them. Each account
  also gets `~/.config/goose/.goosehints` (short replies, no wide tables, code in
  files, no package installs); it is written once, so participants can edit it.
- **Group work** happens in `/group_work`, with `group1` ... `group<group_count>`
  created up front. Every course account can read, write and delete everything
  in it, so trainers have nothing to manage when groups form; participants can
  also create extra directories there. Each home has a `group_work` shortcut to
  it. Like `/home`, it is on the data volume.
