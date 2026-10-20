# Everything secret is generated here and ends up in the Terraform state, so
# keep the state local (or in an encrypted backend) and destroy after the course.

locals {
  usernames = concat(
    [for i in range(1, var.participant_count + 1) : format("user%02d", i)],
    var.extra_users,
  )
}

resource "random_password" "user" {
  for_each = toset(local.usernames)

  length  = 12
  special = false
}

# Per-account key for the LLM proxy. Only works against the proxy inside the
# server, and only up to the account's budget.
resource "random_password" "llm_key" {
  for_each = toset(local.usernames)

  length  = 40
  special = false
}

resource "random_password" "litellm_master_key" {
  length  = 48
  special = false
}

resource "random_password" "litellm_salt_key" {
  length  = 48
  special = false
}

resource "random_password" "postgres" {
  length  = 32
  special = false
}

locals {
  # username <TAB> uid <TAB> password <TAB> llm key; read by the RStudio container
  users_tsv = join("", [
    for i, u in local.usernames :
    "${u}\t${2001 + i}\t${random_password.user[u].result}\tsk-${random_password.llm_key[u].result}\n"
  ])

  # read by the LLM proxy to create the per-account keys
  llm_keys_json = jsonencode({
    model      = var.llm_model_alias
    max_budget = var.llm_budget_per_user_usd
    rpm_limit  = var.llm_rpm_per_user
    users      = [for u in local.usernames : { user = u, key = "sk-${random_password.llm_key[u].result}" }]
  })

  litellm_env = <<-EOT
    LITELLM_MASTER_KEY=sk-${random_password.litellm_master_key.result}
    LITELLM_SALT_KEY=sk-${random_password.litellm_salt_key.result}
    DATABASE_URL=postgresql://litellm:${random_password.postgres.result}@db:5432/litellm
    OPENROUTER_API_KEY=${var.openrouter_api_key}
  EOT

  postgres_env = <<-EOT
    POSTGRES_DB=litellm
    POSTGRES_USER=litellm
    POSTGRES_PASSWORD=${random_password.postgres.result}
  EOT
}
