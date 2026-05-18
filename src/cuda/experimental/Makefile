# cc20 first-kind immune (engine version v16) — standalone build.
#
# Defaults target RTX 5090 / Blackwell (sm_120). Override for other GPUs:
#     make ARCH=sm_89                       # Ada / RTX 4090
#     make ARCH=sm_86                       # Ampere / RTX 3090
#     make NVCC=/usr/local/cuda/bin/nvcc    # non-default CUDA install
#     make OPT="-O3 -g"                     # build with line info

NVCC      ?= nvcc
ARCH      ?= sm_120
OPT       ?= -O3
GMP_LDFLAGS ?= -lgmp

SRC       := cc20_first_kind_immune_v16.cu
BIN_DIR   := bin
BIN       := $(BIN_DIR)/cc20_first_kind_immune_v16

NVCC_FLAGS := $(OPT) -arch=$(ARCH) -std=c++17 \
              -Xcompiler -Wall -Xcompiler -Wno-deprecated-gpu-targets

.PHONY: all clean dirs

all: $(BIN)

dirs:
	@mkdir -p $(BIN_DIR)

$(BIN): $(SRC) | dirs
	$(NVCC) $(NVCC_FLAGS) $< -o $@ $(GMP_LDFLAGS) -lpthread

clean:
	rm -f $(BIN)
	rm -rf $(BIN_DIR)
