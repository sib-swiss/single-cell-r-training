# Server-only layer on top of the course image: adds the goose coding agent
# and creates the course accounts and group directories at container start. Built on the server by
# compose.yaml, so the public image participants run locally stays unchanged.
ARG BASE_IMAGE=sibswiss/training-singlecell-rstudio:2026.3
FROM ${BASE_IMAGE}

ARG GOOSE_VERSION=v1.53.0
RUN apt-get update \
    && apt-get install -y --no-install-recommends acl bzip2 ca-certificates curl \
    && rm -rf /var/lib/apt/lists/* \
    && curl -fsSL -o /tmp/goose.tar.bz2 "https://github.com/aaif-goose/goose/releases/download/${GOOSE_VERSION}/goose-x86_64-unknown-linux-gnu.tar.bz2" \
    && tar -xjf /tmp/goose.tar.bz2 -C /usr/local/bin ./goose \
    && rm /tmp/goose.tar.bz2 \
    && goose --version

# runs after rocker's own 01_set_env and 02_userconf
COPY course_users.sh /etc/cont-init.d/03_course_users
