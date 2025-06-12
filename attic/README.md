# Attic - Dead Code Archive

This directory contains code that was identified as dead (never called) during the Phase 0 dependency analysis conducted on 2025-06-09.

## Purpose

Rather than deleting potentially valuable code, we've moved it here to:
1. Clean up the active R/ directory
2. Maintain a record of past implementation attempts
3. Allow for code archaeology if needed

## Contents

The files here were moved because:
- They contained only functions that were never called in production code
- They represented duplicate implementations of functionality
- They were part of abandoned architectural approaches

## Important Notes

- This code is NOT loaded or sourced by the package
- It is NOT included in package builds
- It exists purely for historical reference
- The code remains in git history regardless

## Files Moved

See the git commit history for details on when each file was moved and why.

---

*Generated during fmriparametry refactoring project*