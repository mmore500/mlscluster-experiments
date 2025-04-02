#!/bin/bash

set -euo pipefail

which air >/dev/null 2>&1 || \
    curl -LsSf "https://github.com/posit-dev/air/releases/tag/0.4.1/download/air-installer.sh" | sh

any_changed=0

for f in R/app.R; do
    # Compute checksum before formatting
    before=$(md5sum "$f")

    # Run formatting
    air format "$f"

    # Compute checksum after formatting
    after=$(md5sum "$f")

    # If the file changed, mark it
    if [ "$before" != "$after" ]; then
        echo "File changed: $f"
        any_changed=1
    fi
done

if [ "$any_changed" -ne 0 ]; then
    echo "Formatting changes detected. Exiting with failure."
    exit 1
else
    echo "No formatting changes detected. Exiting successfully."
    exit 0
fi
