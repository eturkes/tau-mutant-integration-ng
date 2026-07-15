#!/usr/bin/env bash
# Project-local pinned Quarto CLI. Idempotent. Self-contained (bundles deno/pandoc/typst;
# deno needs only baseline glibc, present on trixie -> no apt). Honours no-system-install.
# Repro spec: pinned version + sha256 -> fresh clone regenerates tools/quarto/ identically.
set -euo pipefail

QUARTO_VERSION=1.9.38   # latest STABLE (1.10.x are all prerelease); trust GitHub prerelease flag
QUARTO_SHA256=ea8c897368791ad9f200010c087ea3111b2e556b12a960487dd4e216902aa102
QUARTO_ASSET="quarto-${QUARTO_VERSION}-linux-amd64.tar.gz"
QUARTO_URL="https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/${QUARTO_ASSET}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEST="${ROOT}/tools/quarto/${QUARTO_VERSION}"
BIN="${ROOT}/tools/quarto/bin"
CACHE="${QUARTO_CACHE:-${TMPDIR:-/tmp}}/${QUARTO_ASSET}"   # override dir via QUARTO_CACHE

if [ -f "${DEST}/share/version" ] && [ "$(cat "${DEST}/share/version")" = "${QUARTO_VERSION}" ]; then
  echo "quarto ${QUARTO_VERSION} already installed: ${DEST}"
else
  mkdir -p "${DEST}" "${BIN}"
  if ! { [ -f "${CACHE}" ] && echo "${QUARTO_SHA256}  ${CACHE}" | sha256sum -c - >/dev/null 2>&1; }; then
    echo "downloading quarto ${QUARTO_VERSION} (~136MB) ..."
    rm -f "${CACHE}"   # drop any pre-placed file/symlink at the predictable cache path
    curl -fL -o "${CACHE}" "${QUARTO_URL}"
  fi
  echo "${QUARTO_SHA256}  ${CACHE}" | sha256sum -c -
  tar xzf "${CACHE}" -C "${DEST}" --strip-components=1   # strip leading quarto-<ver>/
fi

# wrapper: launcher hard-errors unless its resolved location is in a dir named `bin` -> symlink ok
ln -sfn "../${QUARTO_VERSION}/bin/quarto" "${BIN}/quarto"
"${BIN}/quarto" --version
