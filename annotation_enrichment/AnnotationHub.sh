#!/bin/bash
# AnnotationHub metadata sqlite3 download script
# Purpose: Resume interrupted download, unlimited automatic retry on network failure
# Database: annotationhub.sqlite3 for offline AnnotationHub cache
# wget arguments explanation
# -c                   Continue partially downloaded files (resume breakpoint)
# --tries=0            Unlimited retry attempts (0 = no limit)
# --timeout=600        Timeout after 600 seconds (10 min), trigger retry
# --retry-connrefused  Retry even if connection is refused

# Modify output directory here
OUTPUT_DIR="/share/home/yuqun/AnnotationHub_download"
OUTPUT_FILE="${OUTPUT_DIR}/annotationhub.sqlite3"
URL="https://annotationhub.bioconductor.org/metadata/annotationhub.sqlite3"

# Create directory if it does not exist
mkdir -p ${OUTPUT_DIR}

# wget download command
wget -c \
  --tries=0 \
  --timeout=600 \
  --retry-connrefused \
  "${URL}" \
  -O "${OUTPUT_FILE}"
