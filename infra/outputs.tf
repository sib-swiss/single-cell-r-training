output "url" {
  value = local.url
}

output "instance_id" {
  description = "For admin access: aws ssm start-session --target <instance_id>"
  value       = aws_instance.server.id
}

# Hand-out list for participants; not committed (see .gitignore).
resource "local_sensitive_file" "credentials" {
  filename = "${path.module}/credentials.csv"
  content = join("", concat(
    ["url,username,password\n"],
    [for u in local.usernames : "${local.url},${u},${random_password.user[u].result}\n"],
  ))
}
