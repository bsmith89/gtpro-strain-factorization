#!/usr/bin/env bash

set -euxo pipefail

git fetch bueno master
git fetch wynton master
git rebase bueno/master
git rebase wynton/master
git push wynton master:origin
git push bueno master:origin
