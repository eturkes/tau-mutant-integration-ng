#!/usr/bin/env bash
# Pinned rv (A2-ai R package manager) on PATH. Idempotent.
# rv's .Rprofile/activate.R discovers rv via Sys.which("rv") + shells out to `rv info`,
# so rv MUST be on PATH. Reproducibility = pinned version + sha256
# below (the committed source of truth), not the install location. cargo-binstall fails
# for rv (crates.io `rv` is an unrelated lib) -> use the release asset.
set -euo pipefail

RV_VERSION=0.22.0
RV_SHA256=3c0bd966193e319863387a74e41e211e6071b4b25222175791fdfc9f016b5f48
RV_ASSET="rv-v${RV_VERSION}-x86_64-unknown-linux-gnu.tar.gz"   # bare `rv` ELF inside
RV_URL="https://github.com/A2-ai/rv/releases/download/v${RV_VERSION}/${RV_ASSET}"

BIN_DIR="${RV_BIN_DIR:-$HOME/.local/bin}"
DEST="${BIN_DIR}/rv"
CACHE="${RV_CACHE:-${TMPDIR:-/tmp}}/${RV_ASSET}"   # download cache (override dir via RV_CACHE)

if [ -x "${DEST}" ] && [ "$("${DEST}" --version 2>/dev/null | awk '{print $NF}')" = "${RV_VERSION}" ]; then
  echo "rv ${RV_VERSION} already installed: ${DEST}"
  exit 0
fi

mkdir -p "${BIN_DIR}"
if ! { [ -f "${CACHE}" ] && echo "${RV_SHA256}  ${CACHE}" | sha256sum -c - >/dev/null 2>&1; }; then
  echo "downloading rv ${RV_VERSION} ..."
  rm -f "${CACHE}"   # drop any pre-placed file/symlink at the predictable cache path
  curl -fL -o "${CACHE}" "${RV_URL}"
fi
echo "${RV_SHA256}  ${CACHE}" | sha256sum -c -
tar xzf "${CACHE}" -C "${BIN_DIR}" rv
chmod +x "${DEST}"

resolved="$(command -v rv 2>/dev/null || true)"
if [ "${resolved}" != "${DEST}" ]; then
  echo "WARNING: 'rv' on PATH = '${resolved:-<none>}', not the pinned ${DEST}." >&2
  echo "         Put ${BIN_DIR} at the front of PATH so rv's R activation uses the pinned binary." >&2
fi
"${DEST}" --version
