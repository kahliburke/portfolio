    #!/usr/bin/env julia
    # ── Run this Kaimon Slate notebook live on your machine ──────────────────────────────────────
    # Auto-generated. Installs Kaimon + KaimonSlate into a dedicated environment, gets this notebook's
    # reproducible bundle, and serves it — the notebook's exact environment reconstructs on open.
    # Re-runnable (idempotent). Prerequisite: Julia 1.10+ (juliaup / https://julialang.org/downloads).
    #
    # Steps are separate functions so this is easy to extend or audit.

    # A dedicated project env keeps this off your default environment.
    const ENVDIR = joinpath(first(DEPOT_PATH), "environments", "kaimonslate-run")

    # Make it active BEFORE the first `using`, so every package this script loads resolves in the
    # one environment. Activating after `using Pkg` instead means Pkg and Downloads are already
    # loaded at the DEFAULT environment's versions, and installing here resolves those same packages
    # differently — Julia then reports dozens of "precompiled but different versions are currently
    # loaded", which it can no longer fix, because loading has happened. `set_active_project` is
    # Base, so this costs no package load, and ENVDIR is only path arithmetic.
    mkpath(ENVDIR); Base.set_active_project(ENVDIR)
    using Pkg, Sockets, Downloads

    # Don't auto-register into the user's Kaimon config — this is a self-contained standalone run.
    ENV["KAIMONSLATE_NO_AUTOREGISTER"] = "1"

    # The reproducible bundle: its filename (shipped beside this script), and — for the published-page
    # one-liner, which runs this script from a temp dir with no sibling — a URL to fetch it from.
    const BUNDLE_NAME = "power_balance_1d.standalone.jl"
    const BUNDLE_URL = get(() -> "", ENV, "SLATE_BUNDLE_URL")

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
    # `default_path` is a checkout baked in at EXPORT time — set when the exporting Slate was itself a
    # working checkout, so the app runs against the build that produced it. The env var still wins.
    # Quiet-git wraps LibGit2 noise.
    function _add_pkg(name, uuid, url, env_path; default_path = "")
        src = strip(get(ENV, env_path, default_path))
        # A baked path that has since moved is worse than useless — it would `develop` a dead
        # directory. Fall through to the normal resolution and say so.
        if !isempty(src) && !isdir(expanduser(String(src)))
            @warn "The $name checkout this app was exported from is gone — resolving $name normally" path = src
            src = ""
        end
        registered = try
            any(r -> haskey(r.pkgs, Base.UUID(uuid)), Pkg.Registry.reachable_registries())
        catch; false; end
        Base.CoreLogging.with_logger(_QuietGit(Base.CoreLogging.current_logger())) do
            if !isempty(src)
                @info "Using a local $name checkout" path = src
                Pkg.develop(path = expanduser(String(src)))
            elseif registered
                spec = Pkg.PackageSpec(name = name, uuid = uuid)
                Pkg.add(spec)                                           # registered → latest release
                # …but ENVDIR persists across runs, and `add` on a package the env already holds at an
                # OLDER release is "already satisfied" — the same no-op documented above for `rev`. So a
                # visitor who ran an earlier bundle keeps that first resolve forever, dependencies and
                # all, long after the package has dropped or changed them. `update` advances it.
                try; Pkg.update(spec); catch e
                    @warn "Could not update $name to the latest release — continuing with the installed version" exception = e
                end
            else
                Pkg.add(url = url, rev = "main")                        # unregistered → track the branch tip
            end
        end
    end

    function install_packages()
        @info "Installing Kaimon + KaimonSlate into $ENVDIR (first run compiles — this can take several minutes)…"
        Pkg.activate(ENVDIR)   # already the active project (set before the first `using`); this re-asserts it
        # Kaimon FIRST: it provides the compute gate the notebook's env reconstructs through (and the
        # agent). KaimonSlate second — both want HTTP 2, so they co-resolve into one env.
        _add_pkg("Kaimon", "d3856c55-31fd-4246-b7e8-380411123c01", "https://github.com/kahliburke/Kaimon.jl", "SLATE_KAIMON_PATH")
        _add_pkg("KaimonSlate", "f7b954f5-0334-4562-ac21-b005218ce1da", "https://github.com/kahliburke/KaimonSlate.jl", "SLATE_KAIMONSLATE_PATH";
                 default_path = raw"/Users/kburke/devel/KaimonSlate.jl")
    end

    # Prefer the bundle shipped NEXT TO this script (extracted site tarball, or downloaded alongside
    # run.jl); only reach out to BUNDLE_URL when there's no sibling (the published-page one-liner).
    function fetch_bundle()
        sib = joinpath(@__DIR__, BUNDLE_NAME)
        isfile(sib) && return sib
        startswith(BUNDLE_URL, "http") ||
            error("Bundle $BUNDLE_NAME not found next to run.jl, and no download URL is set. " *
                  "Put $BUNDLE_NAME in this folder (download it from the page) and re-run.")
        # A TEMP dir, never `pwd()`. The download is a transient input that gets expanded into the
        # install directory — leaving it in whatever folder the one-liner was run from drops a second
        # real notebook `.jl` there, and the expanded copy inherits its document id, so the two read as
        # one document in two places and the reader is asked to split a copy they never made.
        dst = joinpath(mktempdir(; prefix = "slate-bundle-"), BUNDLE_NAME)
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

    # An APP does not ask. "Where should I set this up?" is a question for someone who downloaded a
    # notebook to a temp dir; for an app the folder you are standing in IS the distribution, so the
    # answer is already known. Asking anyway would be worse than noise: a prompt blocks on stdin,
    # which makes the app impossible to start unattended — no systemd unit, no container, no
    # `nohup ./run.sh &`. It expands beside the launcher (a dot-dir, so re-exporting over this
    # folder never mixes bundle files with expanded ones). SLATE_INSTALL_DIR still overrides.
    function app_install_dir()
        pre = strip(get(ENV, "SLATE_INSTALL_DIR", ""))
        isempty(pre) || return abspath(expanduser(String(pre)))
        return joinpath(@__DIR__, ".app")
    end

    # `using` + serve run at TOP LEVEL (not inside a function): the world age advances between
    # top-level statements, so the just-loaded `KaimonSlate` binding is visible for the serve call.
    const NB = setup()
    let dir = choose_install_dir()
        isempty(dir) || (ENV["SLATE_INSTALL_DIR"] = dir; println("→ Setting up in ", dir))
        # Keep THIS run's Slate state — prefs/secrets, the publish ledger, local site builds — OFF the
        # machine-global homes (~/.config, ~/.local/share, ~/.cache under kaimonslate), so a standalone
        # never reads or clobbers a Kaimon extension's (or another standalone's) secrets/ledger/cache.
        # A real install keeps it under the notebook's OWN folder (self-contained + portable; delete the
        # folder → gone); a throwaway ('-') run uses a temp home so it leaves no trace. Because every run
        # gets its own home, multiple standalone hubs can run at once. An explicit KAIMONSLATE_HOME wins.
        if get(ENV, "KAIMONSLATE_HOME", "") == "" && get(ENV, "KAIMONSLATE_CONFIG_HOME", "") == ""
            ENV["KAIMONSLATE_HOME"] = isempty(dir) ? mktempdir(; prefix = "kaimonslate-run-") :
                                                     joinpath(dir, ".kaimonslate")
        end
        # …and the same for KAIMON's cache, which is a separate home and the one that actually leaks.
        # A worker's gate announces itself by dropping a session-metadata file into `<cache>/kaimon/sock`,
        # and a Kaimon running on this machine WATCHES that directory and connects to every local gate it
        # finds there. So on a machine that also runs Kaimon, it adopts THIS run's workers: two clients
        # on one gate, and the hub's calls to its own worker start timing out — which reads as workers
        # "failing" and being respawned, while the worker is healthy and idle throughout.
        #
        # Pointing the cache at this run's own folder announces its gates where nothing else is
        # watching. A throwaway run gets a temp dir, so it still leaves no trace. An explicitly-set
        # XDG_CACHE_HOME / LOCALAPPDATA wins, exactly as with the home above.
        #
        # Note WHERE: beside the install dir, never inside it. `dir` is where the bundle is about to
        # be reconstructed, and that refuses to extract into a non-empty directory which isn't already
        # a Slate install (it must not clobber someone's folder). The cache is populated as soon as
        # the hub touches its history — before the notebook is opened, hence before extraction — so a
        # cache under `dir` makes the install dir non-empty and fails the run with a message about
        # SLATE_INSTALL_DIR that points nowhere near the cause. This folder holds the launcher, so
        # the state stays self-contained: delete the folder and it's all gone.
        let cachevar = Sys.iswindows() ? "LOCALAPPDATA" : "XDG_CACHE_HOME"
            if get(ENV, cachevar, "") == ""
                ENV[cachevar] = isempty(dir) ? mktempdir(; prefix = "kaimonslate-cache-") :
                                               joinpath(@__DIR__, ".cache")
            end
        end
        # A notebook can place cells on a REGION — a named compute target. The name travels in the
        # notebook, the host does not, and the isolated home above has no registry, so seed it from
        # the definitions carried beside this launcher. Only when the home has none of its own, so an
        # operator who configures a region there (pointing the app at their own machine) always wins.
        let src = joinpath(@__DIR__, "regions.json")
            cfg = get(ENV, "KAIMONSLATE_CONFIG_HOME", "")
            isempty(cfg) && (h = get(ENV, "KAIMONSLATE_HOME", ""); cfg = isempty(h) ? "" : joinpath(h, "config"))
            if isfile(src) && !isempty(cfg) && !isfile(joinpath(cfg, "regions.json"))
                try; mkpath(cfg); cp(src, joinpath(cfg, "regions.json")); catch e
                    println("! could not seed region definitions: ", e)
                end
            end
        end
    end
    port = something(tryparse(Int, get(ENV, "SLATE_PORT", "")), free_port())
    host = "127.0.0.1"
    println("Starting the notebook server (first run compiles the environment)…")
    # We only need Kaimon as the compute-gate CLIENT (it spawns/drives the notebook's worker) — pure
    # code, no services. Its __init__ would auto-start a gate SERVER if a kaimon.toml [gate] / KAIMON_GATE_*
    # env is configured (a non-issue on a fresh machine, but a port clash if the viewer already runs
    # Kaimon) — a non-tcp mode makes that auto-serve early-return.
    get(ENV, "KAIMON_GATE_MODE", "") == "" && (ENV["KAIMON_GATE_MODE"] = "off")
    # And don't fire the agent doc-index background service just to VIEW a notebook (it would harvest the
    # whole package doc set). The in-notebook agent can still index on demand.
    get(ENV, "KAIMONSLATE_NO_AUTOINDEX", "") == "" && (ENV["KAIMONSLATE_NO_AUTOINDEX"] = "1")
    # Don't auto-open a browser. Land in the terminal with b/p/q keys instead, so you open in YOUR
    # browser — set SLATE_BROWSER (e.g. "Google Chrome") to override the OS default — and choose whether
    # to go live (b) or just preview the stored render (p). Explicit b/p key opens ignore this flag.
    get(ENV, "KAIMONSLATE_NO_OPEN", "") == "" && (ENV["KAIMONSLATE_NO_OPEN"] = "1")
    import Kaimon        # defines Main.Kaimon → the compute gate that spawns the worker + reconstructs the env
    using KaimonSlate
# serve_notebook blocks; once the hub answers HTTP it prints a framed banner with the URL + keys.
# `inactive=true`: open as a static preview (no worker/precompile) until you press `b` to go live.
KaimonSlate.serve_notebook(NB; port = port, inactive = true)
    