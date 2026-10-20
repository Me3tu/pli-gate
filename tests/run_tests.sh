#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# Run the whole PLI test suite. No arguments needed.
set -euo pipefail
exec python3 "$(dirname "${BASH_SOURCE[0]}")/test_pli.py" "$@"
