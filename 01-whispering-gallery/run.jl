#!/usr/bin/env julia
# ── Run this Kaimon Slate notebook live on your machine ──────────────────────────────────────
# Auto-generated. Installs Kaimon + KaimonSlate into a dedicated environment, gets this notebook's
# reproducible bundle, and serves it — the notebook's exact environment reconstructs on open.
# Re-runnable (idempotent). Prerequisite: Julia 1.10+ (juliaup / https://julialang.org/downloads).
#
# Steps are separate functions so this is easy to extend or audit.
using Pkg, Sockets, Downloads

# Don't auto-register into the user's Kaimon config — this is a self-contained standalone run.
ENV["KAIMONSLATE_NO_AUTOREGISTER"] = "1"

# A dedicated project env keeps this off your default environment.
const ENVDIR = joinpath(first(DEPOT_PATH), "environments", "kaimonslate-run")

# The reproducible bundle: its filename (shipped beside this script), and — for the published-page
# one-liner, which runs this script from a temp dir with no sibling — a URL to fetch it from.
const BUNDLE_NAME = "01_whispering_gallery.standalone.jl"
const BUNDLE_URL = get(() -> "https://kahliburke.github.io/portfolio/01-whispering-gallery/01_whispering_gallery.standalone.jl", ENV, "SLATE_BUNDLE_URL")

# Pick a free TCP port (the default 8765 may already be taken, e.g. by a running Kaimon).
function free_port()
    s = Sockets.listen(Sockets.localhost, 0); p = Int(Sockets.getsockname(s)[2]); close(s); return p
end
function ensure_julia()
    VERSION >= v"1.10" || error("Julia 1.10+ required (found $VERSION). See https://julialang.org/downloads")
end

# LibGit2 (Pkg's git transport) warns about credential attributes newer git sends that its parser
# doesn't know ("Unknown git credential attribute … capability[]/state[]") — harmless noise. Drop
# just those while installing; everything else logs normally.
struct _QuietGit <: Base.CoreLogging.AbstractLogger; sink::Any; end
Base.CoreLogging.min_enabled_level(l::_QuietGit) = Base.CoreLogging.min_enabled_level(l.sink)
Base.CoreLogging.shouldlog(l::_QuietGit, a...) = Base.CoreLogging.shouldlog(l.sink, a...)
Base.CoreLogging.catch_exceptions(l::_QuietGit) = Base.CoreLogging.catch_exceptions(l.sink)
function Base.CoreLogging.handle_message(l::_QuietGit, lvl, msg, _mod, grp, id, file, line; kw...)
    occursin("git credential attribute", string(msg)) && return
    Base.CoreLogging.handle_message(l.sink, lvl, msg, _mod, grp, id, file, line; kw...)
end

# Install one package, always at its LATEST — never a pinned snapshot. Precedence:
#   1. a LOCAL checkout via `env_path` (dev/testing, or a fork), if set;
#   2. the registered RELEASE once the package is published (latest compatible) — detected by the
#      UUID being reachable in a registry, so this SELF-SWITCHES at release with no edit here;
#   3. until then, track `main` on GitHub. `rev="main"` is deliberate: a plain `Pkg.add(url=…)` on a
#      package already pinned to an OLDER `main` commit is a no-op ("already satisfied") and never
#      advances to the tip — pinning the branch re-resolves to the CURRENT tip each run.
# `force_main` skips step 2 even when a registry has the UUID — for a package whose registered
# release is STALE/incompatible with its sibling (so we must track the branch tip until a
# compatible version is registered). Quiet-git wraps LibGit2 noise.
function _add_pkg(name, uuid, url, env_path; force_main::Bool = false)
    src = strip(get(ENV, env_path, ""))
    registered = !force_main && try
        any(r -> haskey(r.pkgs, Base.UUID(uuid)), Pkg.Registry.reachable_registries())
    catch; false; end
    Base.CoreLogging.with_logger(_QuietGit(Base.CoreLogging.current_logger())) do
        if !isempty(src)
            @info "Using a local $name checkout" path = src
            Pkg.develop(path = expanduser(String(src)))
        elseif registered
            Pkg.add(Pkg.PackageSpec(name = name, uuid = uuid))     # registered → latest release
        else
            Pkg.add(url = url, rev = "main")                        # unregistered → track the branch tip
        end
    end
end

function install_packages()
    @info "Installing Kaimon + KaimonSlate into $ENVDIR (first run compiles — this can take several minutes)…"
    mkpath(ENVDIR); Pkg.activate(ENVDIR)
    # Kaimon FIRST: it provides the compute gate the notebook's env reconstructs through (and the
    # agent). KaimonSlate second — both want HTTP 2, so they co-resolve into one env.
    # NOTE: Kaimon's currently-REGISTERED release predates KaimonSlate's move to HTTP 2 (its compat
    # pins HTTP 1.x), so it can't co-resolve with KaimonSlate. Track Kaimon's `main` tip (HTTP 2)
    # for now via force_main — DROP that once a HTTP-2 Kaimon is registered (then it self-switches).
    _add_pkg("Kaimon", "d3856c55-31fd-4246-b7e8-380411123c01", "https://github.com/kahliburke/Kaimon.jl", "SLATE_KAIMON_PATH"; force_main = true)
    _add_pkg("KaimonSlate", "f7b954f5-0334-4562-ac21-b005218ce1da", "https://github.com/kahliburke/KaimonSlate.jl", "SLATE_KAIMONSLATE_PATH")
end

# Prefer the bundle shipped NEXT TO this script (extracted site tarball, or downloaded alongside
# run.jl); only reach out to BUNDLE_URL when there's no sibling (the published-page one-liner).
function fetch_bundle()
    sib = joinpath(@__DIR__, BUNDLE_NAME)
    isfile(sib) && return sib
    startswith(BUNDLE_URL, "http") ||
        error("Bundle $BUNDLE_NAME not found next to run.jl, and no download URL is set. " *
              "Put $BUNDLE_NAME in this folder (download it from the page) and re-run.")
    dst = joinpath(pwd(), BUNDLE_NAME)
    @info "Downloading the notebook bundle" BUNDLE_URL
    Downloads.download(BUNDLE_URL, dst)
    return dst
end

function setup()
    ensure_julia()
    install_packages()
    nb = fetch_bundle()
    println("""
┌ Kaimon is installed alongside KaimonSlate, so the in-notebook AI agent is available.
└ See https://github.com/kahliburke/Kaimon.jl for standalone Kaimon usage.
""")
    return nb
end

# Where to set this notebook up on your machine. A full-history bundle expands to a real GIT
# CHECKOUT (branch & PR with matching commits); a source-only bundle to a plain project. Setting
# SLATE_INSTALL_DIR skips the prompt (non-interactive). '-' runs it once from a throwaway cache.
function choose_install_dir()
    pre = strip(get(ENV, "SLATE_INSTALL_DIR", ""))
    isempty(pre) || return String(pre)
    stem = endswith(BUNDLE_NAME, ".standalone.jl") ? BUNDLE_NAME[1:end-14] : first(splitext(BUNDLE_NAME))
    default = joinpath(pwd(), stem)
    print("\nWhere should I set up this notebook?\n  [", default, "]  ",
          "(Enter = accept · type a path · '-' = run once, don't keep files): ")
    flush(stdout)
    ans = try; strip(readline()); catch; ""; end
    ans == "-" && return ""
    isempty(ans) ? default : abspath(expanduser(String(ans)))
end

# `using` + serve run at TOP LEVEL (not inside a function): the world age advances between
# top-level statements, so the just-loaded `KaimonSlate` binding is visible for the serve call.
const NB = setup()
let dir = choose_install_dir()
    isempty(dir) || (ENV["SLATE_INSTALL_DIR"] = dir; println("→ Setting up in ", dir))
end
port = free_port()
println("Starting the notebook server (first run compiles the environment)…")
# We only need Kaimon as the compute-gate CLIENT (it spawns/drives the notebook's worker) — pure
# code, no services. Its __init__ would auto-start a gate SERVER if a kaimon.toml [gate] / KAIMON_GATE_*
# env is configured (a non-issue on a fresh machine, but a port clash if the viewer already runs
# Kaimon) — a non-tcp mode makes that auto-serve early-return.
get(ENV, "KAIMON_GATE_MODE", "") == "" && (ENV["KAIMON_GATE_MODE"] = "off")
# And don't fire the agent doc-index background service just to VIEW a notebook (it would harvest the
# whole package doc set). The in-notebook agent can still index on demand.
get(ENV, "KAIMONSLATE_NO_AUTOINDEX", "") == "" && (ENV["KAIMONSLATE_NO_AUTOINDEX"] = "1")
import Kaimon        # defines Main.Kaimon → the compute gate that spawns the worker + reconstructs the env
using KaimonSlate
# serve_notebook blocks; once the hub answers HTTP it prints a framed banner with the live URL.
KaimonSlate.serve_notebook(NB; port = port)
