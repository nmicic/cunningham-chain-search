# cc-math MCP Server — Setup Guide

GP/PARI math engine exposed as an MCP server for Claude Code.
Provides Cunningham Chain analysis, search, safe prime generation, and raw GP access.

## Prerequisites

| Dependency | Version | Install |
|---|---|---|
| Node.js | ≥ 18 | `nvm install 20` or system package |
| npm | (bundled) | comes with Node.js |
| GP/PARI | any | `apt install pari-gp` / `brew install pari` |

## Quick Start (Local / stdio)

```bash
cd gp/mcp

# 1. Check environment
./mcp.sh doctor

# 2. Install dependencies
./mcp.sh setup

# 3. Register with Claude Code
./mcp.sh register

# 4. Verify it works
./mcp.sh test
```

After `register`, restart Claude Code — the cc-math tools will appear automatically.

## Manual Setup (if mcp.sh doesn't suit your needs)

```bash
cd gp/mcp
npm install --production
```

Add to `.claude/settings.json` (use absolute paths):

```json
{
  "mcpServers": {
    "cc-math": {
      "command": "node",
      "args": ["/full/path/to/gp/mcp/cc_math_server.js"]
    }
  }
}
```

## Docker (Remote / Headless)

```bash
cd gp/mcp

# Build and start
docker compose up --build -d

# Verify
curl -s http://localhost:3000/health
# → {"status":"ok","server":"cc-math"}

# Stop
docker compose down
```

### Configure Claude Code for remote HTTP

```json
{
  "mcpServers": {
    "cc-math": {
      "url": "http://hostname:3000/mcp"
    }
  }
}
```

Replace `hostname` with the server's IP or hostname.

## HTTP Mode (without Docker)

```bash
# Start (background, port 3000)
./mcp.sh start-http

# Custom port
MCP_PORT=8080 ./mcp.sh start-http

# Stop
./mcp.sh stop
```

Or run directly:

```bash
node cc_math_server.js --http
# or
MCP_HTTP=1 node cc_math_server.js
```

## Environment Variables

| Variable | Default | Description |
|---|---|---|
| `MCP_HTTP` | unset | Set to `1` to enable HTTP mode |
| `MCP_PORT` | `3000` | HTTP listen port |
| `CC_LIB_PATH` | `../cc_lib_v10.gp` | Override path to GP library |

## Troubleshooting

### `./mcp.sh doctor` reports failures

- **Node.js not found**: Install via `nvm`, Homebrew, or system package manager
- **gp not found**: `sudo apt install pari-gp` (Linux) or `brew install pari` (macOS)
- **cc_lib_v10.gp not found**: Ensure you're running from the correct directory; the library should be at `../cc_lib_v10.gp` relative to the `mcp/` folder
- **node_modules not installed**: Run `./mcp.sh setup`

### Server starts but tools fail

- Check GP is working: `echo "isprime(7)" | gp -q` should print `1`
- Check library loads: `echo '\r cc_lib_v10.gp' | gp -q` should show no errors
- Check logs (HTTP mode): `cat mcp-http.log`

### Docker build fails on `pari-gp`

On ARM hosts (Apple Silicon), the `pari-gp` Debian package should work. If not, build with `--platform linux/amd64`.

### Claude Code doesn't see the tools

1. Verify settings: `cat .claude/settings.json`
2. Paths must be absolute (not `~/...`)
3. Restart Claude Code after changing settings
4. Test manually: `./mcp.sh test`

## Available Tools

| Tool | Description |
|---|---|
| `gp_eval` | Run arbitrary GP/PARI code |
| `cc_analyze` | Full chain analysis (root, walk, shadow, classify) |
| `cc_check` | Quick triage — is a number part of a CC? |
| `cc_search` | Search for CC roots at given bit size |
| `cc_safe_prime` | Generate safe primes (delta-sieve algorithm) |
| `cc_immunization` | Shadow spiral residue analysis |
