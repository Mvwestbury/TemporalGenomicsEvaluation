# Race-safe mkdir for cluster scripts: source this right after `set -euo pipefail`.
# The nodes' default mkdir is uutils coreutils 0.10.0, whose `mkdir -p` fails with "mkdir: Already exists" when
# another process creates a directory on the path at the same moment (parallel array tasks on NFS).
# This replaces mkdir with GNU mkdir (/usr/bin/gnumkdir, race-safe -p). Without it: with -p, retry and accept
# directories that already exist; without -p the normal "fail if it exists" behaviour (used for lock dirs) is kept.
# The function is exported, so bash child scripts (e.g. PlainMap) use it too.
mkdir() {
    if [[ -x /usr/bin/gnumkdir ]]; then /usr/bin/gnumkdir "$@"; return; fi
    local a p=0 try
    for a in "$@"; do [[ $a == -p || $a == --parents || ($a == -[!-]* && $a == *p*) ]] && p=1; done
    (( p )) || { command mkdir "$@"; return; }
    for try in 1 2 3 4 5; do
        command mkdir "$@" 2>/dev/null && return 0
        local ok=1
        for a in "$@"; do [[ $a == -* ]] && continue; [[ -d $a ]] || ok=0; done
        (( ok )) && return 0
        sleep 0.$((RANDOM % 9 + 1))
    done
    command mkdir "$@"
}
export -f mkdir
