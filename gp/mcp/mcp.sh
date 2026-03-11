#!/usr/bin/env bash
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
#
# mcp.sh — Unified management script for cc-math MCP server
# Usage: ./mcp.sh {setup|doctor|test|register|start-http|stop}
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SERVER_JS="$SCRIPT_DIR/cc_math_server.js"
LIB_GP="$SCRIPT_DIR/../cc_lib_v10.gp"
PIDFILE="$SCRIPT_DIR/.mcp-http.pid"
PORT="${MCP_PORT:-3000}"

# Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'

ok()   { printf "${GREEN}✓${NC} %s\n" "$1"; }
fail() { printf "${RED}✗${NC} %s\n" "$1"; }
warn() { printf "${YELLOW}!${NC} %s\n" "$1"; }

# ---------------------------------------------------------------------------
# doctor — check prerequisites
# ---------------------------------------------------------------------------
cmd_doctor() {
    local errors=0
    echo "cc-math MCP server — environment check"
    echo "──────────────────────────────────────────"

    # Node.js
    if command -v node &>/dev/null; then
        local nv; nv=$(node --version)
        if [[ "${nv#v}" > "18" ]] || [[ "${nv#v}" == 18* ]]; then
            ok "Node.js $nv"
        else
            fail "Node.js $nv (need ≥18)"; ((errors++))
        fi
    else
        fail "Node.js not found"; ((errors++))
    fi

    # npm
    if command -v npm &>/dev/null; then
        ok "npm $(npm --version)"
    else
        fail "npm not found"; ((errors++))
    fi

    # GP/PARI
    if command -v gp &>/dev/null; then
        local gpv; gpv=$(gp --version-short 2>/dev/null || echo "unknown")
        ok "GP/PARI $gpv"
    else
        fail "GP/PARI (gp) not found — install: apt install pari-gp / brew install pari"; ((errors++))
    fi

    # cc_lib_v10.gp
    if [[ -f "$LIB_GP" ]]; then
        ok "cc_lib_v10.gp found"
    else
        fail "cc_lib_v10.gp not found at $LIB_GP"; ((errors++))
    fi

    # Server file
    if [[ -f "$SERVER_JS" ]]; then
        ok "cc_math_server.js found"
    else
        fail "cc_math_server.js not found"; ((errors++))
    fi

    # node_modules
    if [[ -d "$SCRIPT_DIR/node_modules/@modelcontextprotocol/sdk" ]]; then
        local sdkv; sdkv=$(node --input-type=commonjs -e "console.log(JSON.parse(require('fs').readFileSync('$SCRIPT_DIR/node_modules/@modelcontextprotocol/sdk/package.json')).version)" 2>/dev/null || echo "?")
        ok "MCP SDK v$sdkv installed"
    else
        warn "node_modules not installed — run: ./mcp.sh setup"; ((errors++))
    fi

    echo "──────────────────────────────────────────"
    if (( errors == 0 )); then
        ok "All checks passed"
    else
        fail "$errors issue(s) found"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# setup — install dependencies + verify
# ---------------------------------------------------------------------------
cmd_setup() {
    echo "Installing dependencies..."
    cd "$SCRIPT_DIR"
    npm install --production
    echo ""
    cmd_doctor
}

# ---------------------------------------------------------------------------
# test — smoke test: spawn server, send MCP init, verify response
# ---------------------------------------------------------------------------
cmd_test() {
    echo "Smoke test: spawning server in stdio mode..."

    local init_msg='{"jsonrpc":"2.0","method":"initialize","id":1,"params":{"capabilities":{},"clientInfo":{"name":"mcp-test","version":"1.0.0"},"protocolVersion":"2024-11-05"}}'

    # Cross-platform: use Node.js to spawn with timeout (no coreutils needed)
    local response
    response=$(node -e "
const { spawn } = require('child_process');
const child = spawn('node', ['$SERVER_JS'], { stdio: ['pipe', 'pipe', 'ignore'] });
let out = '';
child.stdout.on('data', d => { out += d.toString(); child.kill(); });
child.stdin.write('$init_msg\n');
child.stdin.end();
const timer = setTimeout(() => { child.kill(); process.exit(1); }, 15000);
child.on('exit', () => { clearTimeout(timer); process.stdout.write(out.split('\n')[0]); });
" 2>/dev/null) || {
        fail "Server failed to respond within 15 seconds"
        return 1
    }

    if echo "$response" | grep -q '"result"'; then
        ok "Server responded with valid MCP init result"
        local info; info=$(echo "$response" | node -e "process.stdin.on('data',d=>{const r=JSON.parse(d);console.log(r.result?.serverInfo?.name+' v'+r.result?.serverInfo?.version)})" 2>/dev/null || echo "?")
        echo "  Server: $info"
    else
        fail "Unexpected response: $response"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# register — auto-generate Claude Code MCP settings
# ---------------------------------------------------------------------------
cmd_register() {
    # Find project root (where CLAUDE.md lives)
    local project_root="$SCRIPT_DIR/../.."
    project_root="$(cd "$project_root" && pwd)"

    local settings_dir="$project_root/.claude"
    local settings_file="$settings_dir/settings.json"

    mkdir -p "$settings_dir"

    local abs_server; abs_server="$(cd "$SCRIPT_DIR" && pwd)/cc_math_server.js"

    if [[ -f "$settings_file" ]]; then
        # Merge into existing settings using Node.js
        node -e "
const fs = require('fs');
const s = JSON.parse(fs.readFileSync('$settings_file', 'utf8'));
if (!s.mcpServers) s.mcpServers = {};
s.mcpServers['cc-math'] = { command: 'node', args: ['$abs_server'] };
fs.writeFileSync('$settings_file', JSON.stringify(s, null, 2) + '\n');
console.log('Updated existing $settings_file');
"
    else
        # Create new settings file
        cat > "$settings_file" <<SETTINGS
{
  "mcpServers": {
    "cc-math": {
      "command": "node",
      "args": ["$abs_server"]
    }
  }
}
SETTINGS
        echo "Created $settings_file"
    fi

    ok "Registered cc-math MCP server"
    echo "  Server: $abs_server"
    echo "  Config: $settings_file"
}

# ---------------------------------------------------------------------------
# start-http — start HTTP mode (for remote / Docker use)
# ---------------------------------------------------------------------------
cmd_start_http() {
    if [[ -f "$PIDFILE" ]]; then
        local old_pid; old_pid=$(cat "$PIDFILE")
        if kill -0 "$old_pid" 2>/dev/null; then
            warn "Server already running (PID $old_pid)"
            return 1
        fi
        rm -f "$PIDFILE"
    fi

    echo "Starting HTTP server on port $PORT..."
    MCP_HTTP=1 MCP_PORT="$PORT" nohup node "$SERVER_JS" --http > "$SCRIPT_DIR/mcp-http.log" 2>&1 &
    local pid=$!
    echo "$pid" > "$PIDFILE"

    # Wait briefly and check it's still alive
    sleep 1
    if kill -0 "$pid" 2>/dev/null; then
        ok "Server started (PID $pid, port $PORT)"
        echo "  Log: $SCRIPT_DIR/mcp-http.log"
        echo "  Stop: ./mcp.sh stop"
    else
        fail "Server exited immediately — check $SCRIPT_DIR/mcp-http.log"
        rm -f "$PIDFILE"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# stop — kill HTTP-mode server
# ---------------------------------------------------------------------------
cmd_stop() {
    if [[ ! -f "$PIDFILE" ]]; then
        warn "No pidfile found — server not running?"
        return 0
    fi

    local pid; pid=$(cat "$PIDFILE")
    if kill -0 "$pid" 2>/dev/null; then
        kill "$pid"
        rm -f "$PIDFILE"
        ok "Stopped server (PID $pid)"
    else
        warn "Process $pid not running — cleaning up pidfile"
        rm -f "$PIDFILE"
    fi
}

# ---------------------------------------------------------------------------
# Main dispatcher
# ---------------------------------------------------------------------------
case "${1:-}" in
    setup)      cmd_setup ;;
    doctor)     cmd_doctor ;;
    test)       cmd_test ;;
    register)   cmd_register ;;
    start-http) cmd_start_http ;;
    stop)       cmd_stop ;;
    *)
        echo "Usage: $0 {setup|doctor|test|register|start-http|stop}"
        echo ""
        echo "Commands:"
        echo "  setup      Install npm dependencies + verify prerequisites"
        echo "  doctor     Check environment (node, gp, lib, deps)"
        echo "  test       Smoke-test: spawn server, send MCP init, verify"
        echo "  register   Auto-configure Claude Code settings.json"
        echo "  start-http Start HTTP mode on port ${PORT} (for remote/Docker)"
        echo "  stop       Stop HTTP-mode server"
        exit 1
        ;;
esac
