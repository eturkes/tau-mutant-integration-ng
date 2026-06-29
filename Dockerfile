FROM rocker/rstudio:4.5.2

# System CLI tooling for agent workflow scripts:
# jq powers .agent/context.sh (the context-window usage gauge).
RUN apt-get update \
 && apt-get install -y --no-install-recommends jq \
 && rm -rf /var/lib/apt/lists/*
