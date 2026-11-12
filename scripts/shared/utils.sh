#!/bin/bash
# Shared bash utilities for analysis scripts

# Logging functions
log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2
}

# Exit with error message
die() {
    log_error "$@"
    exit 1
}

# Check if command exists
require_command() {
    command -v "$1" > /dev/null 2>&1 || die "Required command not found: $1"
}

# Check if file exists
require_file() {
    [[ -f "$1" ]] || die "Required file not found: $1"
}

# Check if directory exists
require_dir() {
    [[ -d "$1" ]] || die "Required directory not found: $1"
}

# Validate VCF path (check for .vcf.gz and .vcf.gz.tbi)
require_vcf() {
    local vcf="$1"
    require_file "$vcf"
    require_file "${vcf}.tbi" || die "VCF index not found: ${vcf}.tbi"
}
