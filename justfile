builddir    := "build/debugoptimized"
debugdir    := "build/debug"
releasedir  := "build/release"
testdir     := "build/test"
snapshotdir := "build/snapshot"

bats_run  := "test/bats/bin/bats -p -T -j "

default:
    @just --list

_setup dir *args:
    [ -d "{{dir}}" ] || meson setup "{{dir}}" {{args}}

[group('build')]
build: (_setup builddir)
    meson compile -C "{{builddir}}"

[group('build')]
release: (_setup releasedir "--buildtype=release")
    meson compile -C "{{releasedir}}"

[group('build')]
debug: (_setup debugdir "--buildtype=debug")
    meson compile -C "{{builddir}}"

[group('build')]
install prefix="/usr/local": (_setup releasedir "--buildtype=release")
    meson configure "{{releasedir}}" --prefix="{{prefix}}"
    meson install -C "{{releasedir}}"

[group('build')]
[doc('Re-run meson configure on existing build dirs')]
reconfigure:
    #!/usr/bin/env bash
    set -euo pipefail
    for d in "{{builddir}}" "{{debugdir}}" "{{releasedir}}" "{{testdir}}" "{{snapshotdir}}"; do
        if [ -d "$d" ]; then meson setup --reconfigure "$d"; fi
    done

[group('build')]
clean:
    rm -rf build

[group('dev')]
run bin="eq3nr" *args: build
    "{{builddir}}/{{bin}}" {{args}}

[group('dev')]
[doc('Re-run a recipe on source changes (e.g. `just watch test`)')]
watch cmd="build" *args:
    watchexec -r -c -e f90,F90,build -- just {{cmd}} {{args}}

[group('test')]
[doc('Library unit tests (args forwarded to meson test)')]
test-unit *args: (_setup builddir)
    meson test -C "{{builddir}}" {{args}}

[group('test')]
[doc('CLI tests (bats); optional [suite] [name-regex]')]
test-cli suite="" filter="" jobs="20":
    #!/usr/bin/env bash
    set -euo pipefail
    target="test"
    if [ -n "{{suite}}" ]; then target="test/{{suite}}.bats"; fi
    if [ -n "{{filter}}" ]; then
        {{bats_run}} "{{jobs}}" -f "{{filter}}" "$target"
    else
        {{bats_run}} "{{jobs}}" "$target"
    fi

[group('test')]
[doc('Generate test snapshots')]
snapshot build="release":
    #!/usr/bin/env bash
    set -euo pipefail
    case "$(arch)" in
      "x86_64")
        ARCH="amd64"
        ;;
      "aarch64")
        ARCH="arm64"
        ;;
      *)
        ARCH="$(arch)"
        ;;
    esac

    BUILD_DIR="{{justfile_directory()}}/{{snapshotdir}}"

    TEST_DIR="{{justfile_directory()}}/test"
    DATA_DIR="${TEST_DIR}/data"
    SNAPSHOTS_DIR="${DATA_DIR}/snapshots"
    OUTPUT_DIR="${SNAPSHOTS_DIR}/${ARCH}"

    printf "OUTPUT_DIR: %s\n" "${OUTPUT_DIR}"

    [ -d "${BUILD_DIR}" ] || meson setup "${BUILD_DIR}" "--buildtype={{build}}" >/dev/null
    meson compile -C "${BUILD_DIR}" >/dev/null
    export PATH="${BUILD_DIR}:${PATH}"

    mkdir -p "${OUTPUT_DIR}/eqpt"
    pushd "${OUTPUT_DIR}/eqpt" >/dev/null 2>&1 || exit 2
    for DATA0 in $(fd -e d0 . "${DATA_DIR}/eqpt"); do
        eqpt "${DATA0}" >/dev/null 2>&1 || exit 3
        perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.po
        perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.po
    done
    popd >/dev/null 2>&1 || exit 4

    for INPUT in $(fd -e 3i . "${DATA_DIR}/eq3nr"); do
        THERMO=$(basename "$(dirname "${INPUT}")")
        DATA1="${OUTPUT_DIR}/eqpt/${THERMO}.d1"
        mkdir -p "${OUTPUT_DIR}/eq3nr/${THERMO}"
        pushd "${OUTPUT_DIR}/eq3nr/${THERMO}" >/dev/null 2>&1 || exit 5
        eq3nr "${DATA1}" "${INPUT}" >/dev/null 2>&1 || exit 6
        perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.3o
        perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.3o
        popd >/dev/null 2>&1 || exit 7
    done

    for INPUT in $(fd -e 6i . "${DATA_DIR}/eq6"); do
        THERMO=$(basename "$(dirname "${INPUT}")")
        DATA1="${OUTPUT_DIR}/eqpt/${THERMO}.d1"
        mkdir -p "${OUTPUT_DIR}/eq6/${THERMO}"
        pushd "${OUTPUT_DIR}/eq6/${THERMO}" >/dev/null 2>&1 || exit 8
        eq6 "${DATA1}" "${INPUT}" >/dev/null 2>&1 || exit 9
        perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.6o
        perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.6o
        rm -f *.6ba *.6bb *.6t
        popd >/dev/null 2>&1 || exit 10
    done

[group('test')]
test: test-cli

[group('test')]
ci: build test
