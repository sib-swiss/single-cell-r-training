# --- network ------------------------------------------------------------------

data "aws_vpc" "default" {
  default = true
}

data "aws_subnet" "default" {
  vpc_id            = data.aws_vpc.default.id
  availability_zone = var.availability_zone
  default_for_az    = true
}

resource "aws_security_group" "server" {
  name        = "${var.name}-server"
  description = "RStudio Server for ${var.name}"
  vpc_id      = data.aws_vpc.default.id

  # 80 is only used for the Let's Encrypt challenge and the redirect to 443
  dynamic "ingress" {
    for_each = [80, 443]
    content {
      from_port   = ingress.value
      to_port     = ingress.value
      protocol    = "tcp"
      cidr_blocks = var.allowed_cidrs
    }
  }

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

# A fixed address, so the URL (and its certificate) survives instance replacement.
# Without a domain of our own, sslip.io resolves 1-2-3-4.sslip.io to 1.2.3.4.
resource "aws_eip" "server" {
  domain = "vpc"
}

locals {
  hostname = "${replace(aws_eip.server.public_ip, ".", "-")}.sslip.io"
  url      = "https://${local.hostname}"
}

# --- admin access (SSM Session Manager, no SSH port) --------------------------

resource "aws_iam_role" "server" {
  name = "${var.name}-server"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "ec2.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy_attachment" "ssm" {
  role       = aws_iam_role.server.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore"
}

resource "aws_iam_instance_profile" "server" {
  name = "${var.name}-server"
  role = aws_iam_role.server.name
}

# --- server -------------------------------------------------------------------

data "aws_ami" "ubuntu" {
  most_recent = true
  owners      = ["099720109477"] # Canonical

  filter {
    name   = "name"
    values = ["ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-amd64-server-*"]
  }
}

resource "aws_ebs_volume" "data" {
  availability_zone = var.availability_zone
  size              = var.data_volume_gb
  type              = "gp3"
  encrypted         = true

  tags = {
    Name = "${var.name}-data"
  }
}

data "cloudinit_config" "server" {
  gzip          = true
  base64_encode = true

  part {
    content_type = "text/cloud-config"
    content = yamlencode({
      package_update = true
      packages       = ["docker.io", "docker-compose-v2"]

      write_files = [
        {
          path        = "/opt/course/.env"
          permissions = "0644"
          content     = <<-EOT
            DATA_VOLUME_ID=${aws_ebs_volume.data.id}
            SITE_ADDRESS=${local.hostname}
            RSTUDIO_BASE_IMAGE=${var.rstudio_base_image}
            GOOSE_VERSION=${var.goose_version}
            LITELLM_VERSION=${var.litellm_version}
            LLM_MODEL_ALIAS=${var.llm_model_alias}
          EOT
        },
        { path = "/opt/course/bootstrap.sh", permissions = "0755", content = file("${path.module}/server/bootstrap.sh") },
        { path = "/opt/course/compose.yaml", permissions = "0644", content = file("${path.module}/server/compose.yaml") },
        { path = "/opt/course/Caddyfile", permissions = "0644", content = file("${path.module}/server/Caddyfile") },
        { path = "/opt/course/create_llm_keys.py", permissions = "0644", content = file("${path.module}/server/create_llm_keys.py") },
        { path = "/opt/course/rstudio/Dockerfile", permissions = "0644", content = file("${path.module}/server/rstudio/Dockerfile") },
        { path = "/opt/course/rstudio/course_users.sh", permissions = "0755", content = file("${path.module}/server/rstudio/course_users.sh") },
        {
          path        = "/opt/course/litellm-config.yaml"
          permissions = "0644"
          content = templatefile("${path.module}/server/litellm-config.yaml.tftpl", {
            model_alias      = var.llm_model_alias
            openrouter_model = var.openrouter_model
          })
        },
        # root-only; mounted into the one container that needs each file
        { path = "/opt/course/secrets/users.tsv", permissions = "0600", content = local.users_tsv },
        { path = "/opt/course/secrets/llm-keys.json", permissions = "0600", content = local.llm_keys_json },
        { path = "/opt/course/secrets/litellm.env", permissions = "0600", content = local.litellm_env },
        { path = "/opt/course/secrets/postgres.env", permissions = "0600", content = local.postgres_env },
      ]

      runcmd = [["/opt/course/bootstrap.sh"]]
    })
  }
}

resource "aws_instance" "server" {
  ami                    = data.aws_ami.ubuntu.id
  instance_type          = var.instance_type
  subnet_id              = data.aws_subnet.default.id
  vpc_security_group_ids = [aws_security_group.server.id]
  iam_instance_profile   = aws_iam_instance_profile.server.name

  user_data_base64            = data.cloudinit_config.server.rendered
  user_data_replace_on_change = true

  # The secrets above sit in the user data. IMDSv2 with a hop limit of 1 keeps
  # them (and the instance role) out of reach of anything inside a container.
  metadata_options {
    http_endpoint               = "enabled"
    http_tokens                 = "required"
    http_put_response_hop_limit = 1
  }

  root_block_device {
    volume_size = var.root_volume_gb
    volume_type = "gp3"
    encrypted   = true
  }

  tags = {
    Name = "${var.name}-server"
  }

  lifecycle {
    # a new Ubuntu AMI should not replace a running course server
    ignore_changes = [ami]
  }
}

resource "aws_volume_attachment" "data" {
  device_name                    = "/dev/sdf"
  volume_id                      = aws_ebs_volume.data.id
  instance_id                    = aws_instance.server.id
  stop_instance_before_detaching = true
}

resource "aws_eip_association" "server" {
  allocation_id = aws_eip.server.id
  instance_id   = aws_instance.server.id
}
