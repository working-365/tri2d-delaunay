#######################################################################################
# Help
#######################################################################################
.PHONY: help
help:
	@echo "Makefile Usage:"
	@echo "  make all"
	@echo "      Run all targets: csim, csynth, cosim, and export."
	@echo ""
	@echo "  make csim"
	@echo "      Run C simulation."
	@echo ""
	@echo "  make csynth"
	@echo "      Run HLS synthesis."
	@echo ""
	@echo "  make cosim"
	@echo "      Run C/RTL co-simulation."
	@echo ""
	@echo "  make export"
	@echo "      Export IP core."
	@echo ""
	@echo "  make clean"
	@echo "      Remove generated files."
	@echo ""

#######################################################################################
# Configuration
#######################################################################################
SHELL := bash

# === 修改成你自己的文件 ===
TOP_FUNCTION := tri2d

SOURCES := tri2d.cpp tri2d.hpp
SOURCES_TB := ${SOURCES} tri2d_tb.cpp

CONFIG := hls_config.cfg
WORK_DIR := hls_component/tri2d

COMMON_OPTIONS := --mode hls --config ${CONFIG} --work_dir ${WORK_DIR}

#######################################################################################
# Targets
#######################################################################################

# Default target
all: csim csynth cosim export

# C Simulation
csim: clean
csim: ${SOURCES_TB} ${CONFIG}
	vitis-run --csim ${COMMON_OPTIONS}

# HLS Synthesis
csynth: ${SOURCES} ${CONFIG}
	v++ --compile ${COMMON_OPTIONS}

# Co-simulation
cosim: ${SOURCES_TB} ${CONFIG}
	vitis-run --cosim ${COMMON_OPTIONS}

# Export IP
export: ${SOURCES} ${CONFIG}
	vitis-run --package ${COMMON_OPTIONS}

#######################################################################################
# Clean
#######################################################################################
.PHONY: clean
clean:
	rm -rf ${WORK_DIR} vitis_c++/ .Xil/ xcd.log
