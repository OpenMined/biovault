#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

print_submodule() {
    local path="$1"
    local indent="$2"
    local prefix="$3"

    if [[ ! -d "$path/.git" && ! -f "$path/.git" ]]; then
        echo -e "${indent}${prefix}${CYAN}${path##*/}/${NC} ${RED}[uninitialized]${NC}"
        return
    fi

    local branch=$(git -C "$path" rev-parse --abbrev-ref HEAD 2>/dev/null)
    local dirty=""
    local dirty_color=""

    if [[ -n $(git -C "$path" status --porcelain -uno 2>/dev/null) ]]; then
        dirty=" [dirty]"
        dirty_color="${RED}"
    else
        dirty_color="${GREEN}"
    fi

    local branch_display=""
    if [[ "$branch" == "HEAD" ]]; then
        local tag=$(git -C "$path" describe --tags --exact-match 2>/dev/null)
        if [[ -n "$tag" ]]; then
            branch_display="${YELLOW}($tag)${NC}"
        else
            local short_sha=$(git -C "$path" rev-parse --short HEAD 2>/dev/null)
            branch_display="${YELLOW}(detached: $short_sha)${NC}"
        fi
    else
        branch_display="${BLUE}[$branch]${NC}"
    fi

    echo -e "${indent}${prefix}${CYAN}${path##*/}/${NC} ${branch_display}${dirty_color}${dirty}${NC}"
}

traverse_submodules() {
    local base_path="$1"
    local indent="$2"
    local is_last="$3"

    local submodules=$(git -C "$base_path" config --file .gitmodules --get-regexp path 2>/dev/null | awk '{print $2}' | sort)

    if [[ -z "$submodules" ]]; then
        return
    fi

    local count=$(echo "$submodules" | wc -l | tr -d ' ')
    local i=0

    while IFS= read -r submodule; do
        ((i++))
        local full_path="$base_path/$submodule"
        local name=$(basename "$submodule")

        local current_prefix="├── "
        local next_indent="${indent}│   "
        if [[ $i -eq $count ]]; then
            current_prefix="└── "
            next_indent="${indent}    "
        fi

        print_submodule "$full_path" "$indent" "$current_prefix"
        traverse_submodules "$full_path" "$next_indent" $([[ $i -eq $count ]] && echo 1 || echo 0)
    done <<< "$submodules"
}

show_tree() {
    root_branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null)
    root_dirty=""
    if [[ -n $(git status --porcelain -uno 2>/dev/null) ]]; then
        root_dirty=" ${RED}[dirty]${NC}"
    fi
    echo -e "${CYAN}$(basename "$(pwd)")/${NC} ${BLUE}[$root_branch]${NC}${root_dirty}"

    traverse_submodules "." ""

    echo ""
    echo -e "${GREEN}Legend:${NC}"
    echo -e "  ${BLUE}[branch]${NC}     - on branch"
    echo -e "  ${YELLOW}(tag)${NC}        - detached at tag"
    echo -e "  ${YELLOW}(detached)${NC}   - detached HEAD"
    echo -e "  ${RED}[dirty]${NC}      - uncommitted changes"
    echo -e "  ${RED}[uninitialized]${NC} - submodule not checked out"
}

collect_dirty_submodules() {
    local base_path="$1"
    local submodules=$(git -C "$base_path" config --file .gitmodules --get-regexp path 2>/dev/null | awk '{print $2}')

    for submodule in $submodules; do
        local full_path="$base_path/$submodule"
        if [[ -d "$full_path/.git" || -f "$full_path/.git" ]]; then
            if [[ -n $(git -C "$full_path" status --porcelain -uno 2>/dev/null) ]]; then
                echo "$full_path"
            fi
            collect_dirty_submodules "$full_path"
        fi
    done
}

do_branch() {
    local branch_name="$1"

    echo -e "${YELLOW}Current state:${NC}"
    echo ""
    show_tree
    echo ""

    local dirty_modules=$(collect_dirty_submodules ".")

    if [[ -z "$dirty_modules" ]]; then
        echo -e "${GREEN}No dirty submodules found.${NC}"
        exit 0
    fi

    if [[ -z "$branch_name" ]]; then
        local root_branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null)
        echo -ne "${YELLOW}Enter branch name${NC} [${BLUE}$root_branch${NC}]: "
        read -r branch_name
        if [[ -z "$branch_name" ]]; then
            branch_name="$root_branch"
        fi
    fi

    echo -e "${CYAN}Dirty submodules that will be branched:${NC}"
    echo "$dirty_modules" | while read -r mod; do
        echo -e "  → ${mod#./}"
    done
    echo ""
    echo -e "${YELLOW}This will create/checkout branch '${branch_name}' in each dirty submodule.${NC}"
    echo -ne "${YELLOW}Continue? [y/N]: ${NC}"
    read -r confirm

    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
        echo -e "${BLUE}Aborted.${NC}"
        exit 0
    fi

    echo ""
    echo "$dirty_modules" | while read -r mod; do
        echo -ne "  → ${mod#./}: "
        if git -C "$mod" show-ref --verify --quiet "refs/heads/$branch_name"; then
            git -C "$mod" checkout "$branch_name" 2>/dev/null
            echo -e "${BLUE}checked out existing branch${NC}"
        else
            git -C "$mod" checkout -b "$branch_name" 2>/dev/null
            echo -e "${GREEN}created new branch${NC}"
        fi
    done

    echo ""
    echo -e "${GREEN}Done! New state:${NC}"
    echo ""
    show_tree
}

collect_uninitialized_submodules() {
    local base_path="$1"
    local submodules=$(git -C "$base_path" config --file .gitmodules --get-regexp path 2>/dev/null | awk '{print $2}')

    for submodule in $submodules; do
        local full_path="$base_path/$submodule"
        if [[ ! -d "$full_path/.git" && ! -f "$full_path/.git" ]]; then
            echo "$full_path"
        else
            collect_uninitialized_submodules "$full_path"
        fi
    done
}

do_init() {
    echo -e "${YELLOW}Current state:${NC}"
    echo ""
    show_tree
    echo ""

    local uninit_modules=$(collect_uninitialized_submodules ".")

    if [[ -z "$uninit_modules" ]]; then
        echo -e "${GREEN}All submodules are already initialized.${NC}"
        exit 0
    fi

    echo -e "${CYAN}Uninitialized submodules that will be initialized:${NC}"
    echo "$uninit_modules" | while read -r mod; do
        echo -e "  → ${mod#./}"
    done
    echo ""
    echo -e "${YELLOW}This will run 'git submodule update --init' recursively.${NC}"
    echo -e "${YELLOW}No existing data will be modified or deleted.${NC}"
    echo -ne "${YELLOW}Continue? [y/N]: ${NC}"
    read -r confirm

    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
        echo -e "${BLUE}Aborted.${NC}"
        exit 0
    fi

    echo ""
    echo -e "${CYAN}Initializing submodules recursively...${NC}"

    local pass=1
    while true; do
        uninit_modules=$(collect_uninitialized_submodules ".")
        if [[ -z "$uninit_modules" ]]; then
            break
        fi

        echo -e "${YELLOW}Pass $pass:${NC}"
        echo "$uninit_modules" | while read -r mod; do
            local parent_dir=$(dirname "$mod")
            local submod_name=$(basename "$mod")
            echo -ne "  → ${mod#./}: "
            if git -C "$parent_dir" submodule update --init "$submod_name" 2>/dev/null; then
                echo -e "${GREEN}initialized${NC}"
            else
                echo -e "${RED}failed${NC}"
            fi
        done

        ((pass++))
        if [[ $pass -gt 10 ]]; then
            echo -e "${RED}Too many passes, possible circular dependency. Aborting.${NC}"
            break
        fi
    done

    echo ""
    echo -e "${GREEN}Done! New state:${NC}"
    echo ""
    show_tree
}

do_reset() {
    echo -e "${YELLOW}Current state:${NC}"
    echo ""
    show_tree
    echo ""
    echo -e "${RED}WARNING: This will:${NC}"
    echo -e "  1. Discard all uncommitted changes in submodules"
    echo -e "  2. Remove untracked files in submodules"
    echo -e "  3. Fetch, checkout 'main' (or 'master'), and pull in all submodules"
    echo -e "  4. Initialize any uninitialized submodules"
    echo ""
    echo -ne "${YELLOW}Are you sure you want to reset all submodules? [y/N]: ${NC}"
    read -r confirm

    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
        echo -e "${BLUE}Aborted.${NC}"
        exit 0
    fi

    echo ""
    echo -e "${CYAN}Resetting submodules...${NC}"
    echo ""

    git submodule foreach --recursive '
        echo "  → $name"
        git checkout -- . 2>/dev/null
        git clean -fd 2>/dev/null
        git fetch origin 2>/dev/null
        if git show-ref --verify --quiet refs/heads/main; then
            git checkout main 2>/dev/null
            git pull origin main 2>/dev/null
        elif git show-ref --verify --quiet refs/heads/master; then
            git checkout master 2>/dev/null
            git pull origin master 2>/dev/null
        elif git show-ref --verify --quiet refs/remotes/origin/main; then
            git checkout -b main origin/main 2>/dev/null
        elif git show-ref --verify --quiet refs/remotes/origin/master; then
            git checkout -b master origin/master 2>/dev/null
        else
            echo "    (no main/master branch found locally or on origin)"
        fi
    '

    echo ""
    echo -e "${CYAN}Updating submodules...${NC}"
    git submodule update --init --recursive

    echo ""
    echo -e "${GREEN}Done! New state:${NC}"
    echo ""
    show_tree
}

case "$1" in
    --reset)
        do_reset
        ;;
    --branch)
        do_branch "$2"
        ;;
    --init)
        do_init
        ;;
    --help|-h)
        echo "Usage: ./mods.sh [OPTIONS]"
        echo ""
        echo "Options:"
        echo "  (none)              Show submodule tree with branch/dirty status"
        echo "  --branch [name]     Checkout/create branch in all dirty submodules"
        echo "                      (prompts for name if not provided, defaults to root branch)"
        echo "  --init              Initialize uninitialized submodules (safe, no data loss)"
        echo "  --reset             Reset all submodules to main/master (DESTRUCTIVE)"
        echo "  --help, -h          Show this help"
        ;;
    *)
        show_tree
        ;;
esac
