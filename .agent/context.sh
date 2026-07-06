#!/bin/sh
# Codex context gauge -> "N% used/window" from newest Codex JSONL session for this cwd.
# Override: CODEX_SESSION_FILE=/path/session.jsonl CODEX_CONTEXT_WINDOW=200000 .agent/context.sh
root="${CODEX_HOME:-$HOME/.codex}/sessions"
cwd=$(pwd -P)
f="${CODEX_SESSION_FILE:-}"

if [ -z "$f" ] && [ -d "$root" ]; then
  f=$(
    find "$root" -type f -name '*.jsonl' -printf '%T@ %p\n' 2>/dev/null |
      sort -rn |
      while IFS= read -r row; do
        p=${row#* }
        first=$(sed -n '1p' "$p" 2>/dev/null)
        pcwd=$(printf '%s\n' "$first" | jq -r '(.payload.cwd // .cwd // empty)' 2>/dev/null)
        [ "$pcwd" = "$cwd" ] && { printf '%s\n' "$p"; break; }
      done
  )
fi

u=""
cw=""
if [ -n "$f" ] && [ -r "$f" ]; then
  u=$(
    jq -s '
      [ .[] |
        select(.type == "event_msg" and .payload.type == "token_count") |
        .payload.info |
        (.last_token_usage.input_tokens // .total_token_usage.input_tokens // empty)
      ] | last // empty
    ' "$f" 2>/dev/null
  )
  cw=$(
    jq -s '
      [ .[] |
        select(.type == "event_msg" and .payload.type == "token_count") |
        (.payload.info.model_context_window // empty)
      ] | last // empty
    ' "$f" 2>/dev/null
  )

  if [ -z "$u" ]; then
    u=$(
      jq -s '
        def usage:
          .message.usage? // .payload.message.usage? // .payload.usage? //
          .response.usage? // .payload.response.usage? // .usage?;
        def total:
          (.input_tokens // .prompt_tokens // 0) +
          (.cache_creation_input_tokens // 0) +
          (.cache_read_input_tokens // 0);
        [ .[] | usage | select(type == "object") | total ] | last // empty
      ' "$f" 2>/dev/null
    )
  fi
fi

w="${CODEX_CONTEXT_WINDOW:-${cw:-200000}}"
awk -v u="$u" -v w="$w" '
function h(n){ if(n>=1000000){s=sprintf("%.1fM",n/1000000);sub(/\.0M$/,"M",s);return s}
              return sprintf("%dK",int(n/1000+0.5)) }
BEGIN{ if(u==""){ print "? ?/" h(w); exit }
       print int(u*100/w+0.5) "% " h(u) "/" h(w) }'
