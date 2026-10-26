builddir   := "build/debugoptimized"
debugdir   := "build/debug"
releasedir := "build/release"
testdir    := "build/test"

bats_run  := "./test/bats/bin/bats -p -T -j "

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
    for d in "{{builddir}}" "{{testdir}}" "{{releasedir}}"; do
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
test: test-cli

[group('test')]
ci: build test
