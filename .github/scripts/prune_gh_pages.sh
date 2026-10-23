#!/usr/bin/env bash
# Usage: prune_gh_pages.sh <gh-pages clone> <file with closed PR numbers>
#
# Delete the previews of closed PRs, and replace each release directory that
# is not the latest patch of its minor version with a symlink to that patch.
# Like Documenter with `forcepush = true`, push one commit without history.
set -euo pipefail

repo=$1
closed_prs=$(cat "$2")
index=$(mktemp -d)/index

for attempt in $(seq 10); do
    git -C "$repo" fetch --quiet --depth=1 origin gh-pages
    base=$(git -C "$repo" rev-parse FETCH_HEAD)
    rm -f "$index"
    export GIT_INDEX_FILE=$index
    git -C "$repo" read-tree "$base"

    for pr in $closed_prs; do
        git -C "$repo" rm -r --cached --sparse --quiet --ignore-unmatch "previews/PR$pr"
    done

    superseded=$(git -C "$repo" ls-tree "$base" |
        awk '$2 == "tree" && $4 ~ /^v[0-9]+\.[0-9]+\.[0-9]+$/ { print $4 }' |
        sort -V |
        awk -F. '{ release[NR] = $0; minor[NR] = $1 "." $2; latest[$1 "." $2] = $0 }
            END { for (i = 1; i <= NR; i++) if (release[i] != latest[minor[i]]) print release[i], latest[minor[i]] }')
    while read -r release latest; do
        [ -n "$release" ] || continue
        git -C "$repo" rm -r --cached --sparse --quiet "$release"
        link=$(printf '%s' "$latest" | git -C "$repo" hash-object -w --stdin)
        git -C "$repo" update-index --add --cacheinfo "120000,$link,$release"
    done <<<"$superseded"

    tree=$(git -C "$repo" write-tree)
    unset GIT_INDEX_FILE
    if [ "$tree" = "$(git -C "$repo" rev-parse "$base^{tree}")" ]; then
        echo "Nothing to prune"
        exit 0
    fi
    commit=$(git -C "$repo" commit-tree "$tree" -m "Prune previews and superseded releases")
    if git -C "$repo" push --force-with-lease="gh-pages:$base" origin "$commit:refs/heads/gh-pages"; then
        exit 0
    fi
    echo "Push attempt $attempt was rejected, retrying"
    sleep $((10 + RANDOM % 20))
done
exit 1
