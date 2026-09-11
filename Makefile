# LexiPair-BioRange | Makefile (MinGW / MSYS2 / WSL fallback)
# Slide 17: C++ compiled with GCC -O3, single-threaded.
CXX ?= g++
CXXFLAGS = -O3 -std=c++17 -Iinclude -DNDEBUG
SRC = src/lexipair_index.cpp src/lsm_index.cpp src/fractional_cascade.cpp
APPS = build_index static_query dynamic_update fc_benchmark scalability ablation vgp_compare validate
BIN = bin

all: $(APPS)

$(BIN):
	mkdir -p $(BIN)

%: apps/%.cpp $(SRC) | $(BIN)
	$(CXX) $(CXXFLAGS) apps/$@.cpp $(SRC) -o $(BIN)/$@

clean:
	rm -rf $(BIN)
