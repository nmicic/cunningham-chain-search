#!/usr/bin/env node
/**
 * Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
 * SPDX-License-Identifier: Apache-2.0
 *
 * cc-math MCP Server
 *
 * Wraps GP/PARI with cc_lib_v10.gp loaded, exposing Cunningham Chain
 * analysis tools as MCP resources and tools.
 *
 * Tools:
 *   gp_eval(code)          — run GP code, return output
 *   cc_analyze(number)      — full chain analysis of a number
 *   cc_search(bits, target) — search for CC roots
 *   cc_check(number)        — quick "is this a CC root?" check
 *   cc_explain(topic)       — explain a CC concept with live examples
 */

import { McpServer } from "@modelcontextprotocol/sdk/server/mcp.js";
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js";
import { StreamableHTTPServerTransport } from "@modelcontextprotocol/sdk/server/streamableHttp.js";
import { z } from "zod";
import { spawn } from "node:child_process";
import { createServer } from "node:http";
import { resolve, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import { randomUUID } from "node:crypto";

const __dirname = dirname(fileURLToPath(import.meta.url));
const LIB_PATH = process.env.CC_LIB_PATH || resolve(__dirname, "..", "cc_lib_v10.gp");

// ---------------------------------------------------------------------------
// GP process management — persistent session, load library once
// ---------------------------------------------------------------------------

let gpProcess = null;
let gpReady = false;

const PROMPT_MARKER = "___MCP_READY___";
const MARKER_CMD = `printf("${PROMPT_MARKER}\\n");\n`;

function startGP() {
  if (gpProcess && !gpProcess.killed) return;

  gpProcess = spawn("gp", ["-q", "--default", "parisizemax=512000000"], {
    stdio: ["pipe", "pipe", "pipe"],
  });

  gpProcess.on("exit", () => {
    gpProcess = null;
    gpReady = false;
  });

  gpProcess.on("error", (err) => {
    console.error("GP process error:", err.message);
    gpProcess = null;
    gpReady = false;
  });
}

/**
 * Send code to GP and collect output until we see the marker.
 * Returns the captured output (without the marker line).
 */
function gpEval(code, timeoutMs = 60000) {
  return new Promise((resolve, reject) => {
    if (!gpProcess || gpProcess.killed) {
      startGP();
    }

    let output = "";
    let errOutput = "";
    let settled = false;

    const timer = setTimeout(() => {
      if (!settled) {
        settled = true;
        reject(new Error(`GP timed out after ${timeoutMs / 1000}s`));
      }
    }, timeoutMs);

    const onData = (chunk) => {
      output += chunk.toString();
      if (output.includes(PROMPT_MARKER)) {
        settled = true;
        clearTimeout(timer);
        gpProcess.stdout.removeListener("data", onData);
        gpProcess.stderr.removeListener("data", onErr);
        // Strip everything after (and including) the marker
        const idx = output.indexOf(PROMPT_MARKER);
        const result = output.slice(0, idx).trimEnd();
        resolve(result);
      }
    };

    const onErr = (chunk) => {
      errOutput += chunk.toString();
    };

    gpProcess.stdout.on("data", onData);
    gpProcess.stderr.on("data", onErr);

    // Send code followed by marker
    const safeCode = code.replace(/\n*$/, "\n");
    gpProcess.stdin.write(safeCode + MARKER_CMD);
  });
}

/**
 * Ensure GP is running and library is loaded.
 */
async function ensureReady() {
  if (gpReady) return;
  startGP();
  // Load library — capture and discard self-test output
  await gpEval(`\\r ${LIB_PATH}`);
  gpReady = true;
}

// ---------------------------------------------------------------------------
// Helper: run GP and return clean output
// ---------------------------------------------------------------------------

async function runGP(code, timeoutMs = 60000) {
  await ensureReady();
  return gpEval(code, timeoutMs);
}

// ---------------------------------------------------------------------------
// MCP Server factory — registers all tools and resources on a new McpServer.
// The GP backend (runGP) is shared across all instances.
// ---------------------------------------------------------------------------

function createMcpServer() {
  const srv = new McpServer({
    name: "cc-math",
    version: "1.0.0",
  });

  // --- Tool: gp_eval ---
  srv.tool(
    "gp_eval",
    "Run arbitrary GP/PARI code with cc_lib_v10 loaded. " +
      "Use for any Cunningham Chain computation. " +
      "The library has 37 functions across 10 layers: " +
      "chain analysis, CRT search, constructor search, safe primes, " +
      "delta-sieve, bitmask sieve, BiTwin, primorial scan.",
    { code: z.string().describe("GP/PARI code to execute") },
    async ({ code }) => {
      try {
        const result = await runGP(code);
        return { content: [{ type: "text", text: result || "(no output)" }] };
      } catch (e) {
        return {
          content: [{ type: "text", text: `Error: ${e.message}` }],
          isError: true,
        };
      }
    }
  );

  // --- Tool: cc_analyze ---
  srv.tool(
    "cc_analyze",
    "Full analysis of a number as a Cunningham Chain member. " +
      "Finds root, walks chain, shows shadow spiral, and classifies. " +
      "Input can be decimal or hex (prefix 0x).",
    {
      number: z.string().describe("Prime number to analyze (decimal or 0x hex)"),
      kind: z
        .number()
        .optional()
        .default(1)
        .describe("Chain kind: 1=first (2p+1), 2=second (2p-1)"),
    },
    async ({ number, kind }) => {
      try {
        const code =
          `{my(p = ${number}, R = cc_root(p, ${kind}));` +
          `printf("Input: %s\\n", Str(p));` +
          `printf("Root:  %s (back %d steps)\\n\\n", Str(R[1]), R[2]);` +
          `cc_walk(R[1], ${kind});` +
          `printf("\\n");` +
          `cc_shadow(R[1], ${kind});` +
          `printf("\\n");` +
          `cc_classify(p)}`;
        const result = await runGP(code);
        return { content: [{ type: "text", text: result }] };
      } catch (e) {
        return {
          content: [{ type: "text", text: `Error: ${e.message}` }],
          isError: true,
        };
      }
    }
  );

  // --- Tool: cc_check ---
  srv.tool(
    "cc_check",
    "Quick check: is this number part of a Cunningham Chain? " +
      "Returns chain length, whether it's a root, and classification.",
    {
      number: z.string().describe("Number to check (decimal or 0x hex)"),
    },
    async ({ number }) => {
      try {
        const code =
          `{my(p = ${number}, len = cc_chain_test(p), is_root = cc_is_root(p));` +
          `printf("Number: %s\\n", Str(p));` +
          `printf("Bits:   %d\\n", #binary(p));` +
          `printf("Prime:  %s\\n", if(isprime(p), "yes", "no"));` +
          `printf("CC1a length: %d\\n", len);` +
          `printf("Is CC1a root: %s\\n", if(is_root, "yes", "no"));` +
          `if(!is_root && len > 0, my(R = cc_root(p)); printf("Root: %s (back %d steps)\\n", Str(R[1]), R[2]));` +
          `my(len2 = cc_chain_test(p, 2));` +
          `printf("CC1b length: %d\\n", len2);` +
          `cc_classify(p)}`;
        const result = await runGP(code);
        return { content: [{ type: "text", text: result }] };
      } catch (e) {
        return {
          content: [{ type: "text", text: `Error: ${e.message}` }],
          isError: true,
        };
      }
    }
  );

  // --- Tool: cc_search ---
  srv.tool(
    "cc_search",
    "Search for Cunningham Chain roots at a given bit size. " +
      "Uses CRT + sieve + PRP pipeline.",
    {
      bits: z.number().describe("Target bit size (e.g. 30, 40, 60, 89)"),
      target: z
        .number()
        .describe("Minimum chain length to find (e.g. 4, 5, 7)"),
      max_candidates: z
        .number()
        .optional()
        .default(50000)
        .describe("Maximum candidates to test"),
    },
    async ({ bits, target, max_candidates }) => {
      try {
        const code =
          `{my(hits = cc_search(${bits}, ${target}, 1, 0, 0, ${max_candidates}, 1));` +
          `printf("\\nFound %d chains of CC%d+ at %d bits\\n", #hits, ${target}, ${bits});` +
          `for(i = 1, min(5, #hits), printf("  [%d] CC%d root: %s (%d bits)\\n", i, hits[i][2], Str(hits[i][1]), #binary(hits[i][1])));` +
          `if(#hits > 5, printf("  ... and %d more\\n", #hits - 5))}`;
        const result = await runGP(code, 120000);
        return { content: [{ type: "text", text: result }] };
      } catch (e) {
        return {
          content: [{ type: "text", text: `Error: ${e.message}` }],
          isError: true,
        };
      }
    }
  );

  // --- Tool: cc_safe_prime ---
  srv.tool(
    "cc_safe_prime",
    "Generate a safe prime (p = 2q+1 where both p and q are prime). " +
      "Uses delta-sieve algorithm that beats OpenSSL.",
    {
      bits: z
        .number()
        .describe("Bit size of safe prime p (e.g. 64, 256, 512, 1024)"),
    },
    async ({ bits }) => {
      try {
        const timeout = bits <= 256 ? 30000 : bits <= 1024 ? 120000 : 300000;
        const code =
          `{my(r = cc_sp_search(${bits}));` +
          `if(r, printf("Safe prime found:\\n  q = %s\\n  p = 2q+1 = %s\\n  bits(p) = %d\\n", Str(r[1]), Str(r[2]), r[3]), printf("Not found within candidate limit\\n"))}`;
        const result = await runGP(code, timeout);
        return { content: [{ type: "text", text: result }] };
      } catch (e) {
        return {
          content: [{ type: "text", text: `Error: ${e.message}` }],
          isError: true,
        };
      }
    }
  );

  // --- Tool: cc_immunization ---
  srv.tool(
    "cc_immunization",
    "Analyze immunization of a Cunningham Chain root — which small primes " +
      "can never kill the chain and why (residue analysis).",
    {
      number: z.string().describe("Chain root to analyze (decimal or 0x hex)"),
      kind: z.number().optional().default(1).describe("Chain kind: 1 or 2"),
    },
    async ({ number, kind }) => {
      try {
        const code =
          `{my(p = ${number});` +
          `printf("Immunization analysis for %s\\n\\n", Str(p));` +
          `cc_shadow(p, ${kind});` +
          `printf("\\n");` +
          `my(v = cc_walk_quiet(p, ${kind}), len = v[2]);` +
          `printf("Chain length: %d\\n", len);` +
          `printf("\\nImmune primes (residue = q-1 for all analysis primes q):\\n");` +
          `foreach(cc_analysis_primes, q, if(p % q == q - 1, printf("  prime %d: residue %d = %d-1 (IMMUNE)\\n", q, p % q, q)));` +
          `printf("\\nVulnerable primes:\\n");` +
          `foreach(cc_analysis_primes, q, if(p % q != q - 1, my(kp = cc_kill_pos(p % q, q, ${kind})); printf("  prime %d: residue %d, kills at position %d\\n", q, p % q, kp)))}`;
        const result = await runGP(code);
        return { content: [{ type: "text", text: result }] };
      } catch (e) {
        return {
          content: [{ type: "text", text: `Error: ${e.message}` }],
          isError: true,
        };
      }
    }
  );

  // --- Resource: function reference ---
  srv.resource(
    "cc-functions",
    "cc://functions",
    async () => {
      try {
        const result = await runGP("cc_help()");
        return {
          contents: [
            {
              uri: "cc://functions",
              mimeType: "text/plain",
              text: result,
            },
          ],
        };
      } catch (e) {
        return {
          contents: [
            {
              uri: "cc://functions",
              mimeType: "text/plain",
              text: `Error loading help: ${e.message}`,
            },
          ],
        };
      }
    }
  );

  return srv;
}

// ---------------------------------------------------------------------------
// Start — stdio (default) or HTTP (--http flag / MCP_HTTP env)
// ---------------------------------------------------------------------------

const useHttp = process.argv.includes("--http") || process.env.MCP_HTTP === "1";

if (useHttp) {
  const port = parseInt(process.env.MCP_PORT || "3000", 10);
  const sessions = new Map();  // sessionId → { server, transport }

  const httpServer = createServer(async (req, res) => {
    // Health check
    if (req.method === "GET" && req.url === "/health") {
      res.writeHead(200, { "Content-Type": "application/json" });
      res.end(JSON.stringify({ status: "ok", server: "cc-math" }));
      return;
    }

    // MCP endpoint
    if (req.url === "/mcp") {
      // Check for existing session
      const sessionId = req.headers["mcp-session-id"];
      if (sessionId && sessions.has(sessionId)) {
        const session = sessions.get(sessionId);
        await session.transport.handleRequest(req, res);
        return;
      }

      // New session: fresh McpServer + transport
      const transport = new StreamableHTTPServerTransport({
        sessionIdGenerator: () => randomUUID(),
      });
      const server = createMcpServer();
      transport.onclose = () => {
        if (transport.sessionId) sessions.delete(transport.sessionId);
      };
      await server.connect(transport);
      await transport.handleRequest(req, res);
      // Session ID is assigned after handleRequest processes the init
      if (transport.sessionId) {
        sessions.set(transport.sessionId, { server, transport });
      }
      return;
    }

    res.writeHead(404);
    res.end("Not found");
  });

  httpServer.listen(port, () => {
    console.log(`cc-math MCP server (HTTP) listening on port ${port}`);
    console.log(`  MCP endpoint: http://localhost:${port}/mcp`);
    console.log(`  Health check: http://localhost:${port}/health`);
  });
} else {
  const server = createMcpServer();
  const transport = new StdioServerTransport();
  await server.connect(transport);
}
