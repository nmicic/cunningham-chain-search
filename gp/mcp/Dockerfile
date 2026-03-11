# Build context must be gp/ (set via docker-compose)
FROM node:20-slim

# Install GP/PARI
RUN apt-get update && \
    apt-get install -y --no-install-recommends pari-gp && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Node.js dependencies
COPY mcp/package.json mcp/package-lock.json ./
RUN npm ci --production

# Copy GP library
COPY cc_lib_v10.gp ./lib/cc_lib_v10.gp

# Copy server + knowledge base
COPY mcp/cc_math_server.js ./
COPY mcp/CC_KNOWLEDGE.md ./

EXPOSE 3000

# LIB_PATH override: in container, lib is at /app/lib/cc_lib_v10.gp
ENV CC_LIB_PATH=/app/lib/cc_lib_v10.gp
ENV MCP_HTTP=1
ENV MCP_PORT=3000

CMD ["node", "cc_math_server.js", "--http"]
