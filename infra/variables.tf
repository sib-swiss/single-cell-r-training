variable "name" {
  description = "Name prefix for all AWS resources; use one per course run, e.g. sc-2026-11."
  type        = string
  default     = "single-cell-course"
}

variable "region" {
  type    = string
  default = "eu-central-1"
}

variable "availability_zone" {
  description = "The instance and the data volume must live in the same AZ."
  type        = string
  default     = "eu-central-1a"
}

# --- accounts -----------------------------------------------------------------

variable "participant_count" {
  description = "Number of participant accounts (user01, user02, ...). Add a few spares: adding accounts later replaces the instance."
  type        = number
  default     = 33
}

variable "extra_users" {
  description = "Additional named accounts, e.g. for trainers. Same rights as participants."
  type        = list(string)
  default     = []
}

# --- server -------------------------------------------------------------------

variable "instance_type" {
  description = "Rule of thumb: 12-16 GiB RAM per participant. r7i.16xlarge = 64 vCPU / 512 GiB."
  type        = string
  default     = "r7i.16xlarge"
}

variable "root_volume_gb" {
  description = "Holds the OS and Docker images."
  type        = number
  default     = 100
}

variable "data_volume_gb" {
  description = "Holds /home, the LiteLLM database and TLS certificates. Survives instance replacement."
  type        = number
  default     = 600
}

variable "allowed_cidrs" {
  description = "Who may reach the web interface (ports 80/443)."
  type        = list(string)
  default     = ["0.0.0.0/0"]
}

# --- software -----------------------------------------------------------------

variable "rstudio_base_image" {
  description = "Course image the server image is built from (see Docker/Dockerfile)."
  type        = string
  default     = "sibswiss/training-singlecell-rstudio:2026.3"
}

variable "goose_version" {
  type    = string
  default = "v1.53.0"
}

variable "litellm_version" {
  type    = string
  default = "v1.102.4"
}

# --- LLM access ---------------------------------------------------------------

variable "openrouter_api_key" {
  description = "A key made for this course only, with a credit limit set on openrouter.ai. Never reaches the participants."
  type        = string
  sensitive   = true
}

variable "openrouter_model" {
  description = "OpenRouter model id."
  type        = string
  default     = "deepseek/deepseek-v4.1-flash-1017"
}

variable "llm_model_alias" {
  description = "Model name participants (and goose) see."
  type        = string
  default     = "deepseek-v4.1-flash"
}

variable "llm_budget_per_user_usd" {
  description = "Spend limit per account over the whole course."
  type        = number
  default     = 10
}

variable "llm_rpm_per_user" {
  description = "Requests per minute per account."
  type        = number
  default     = 60
}
