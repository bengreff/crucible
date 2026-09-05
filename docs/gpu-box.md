# GPU box `backhouse` — access + usage

The CRUCIBLE GPU tier (plan ruling #9; brought up in session 24 / plan S12). This is the
authoritative how-to for reaching and driving the box; SESSION_LOG session 24 has the story and the
first measured numbers. Read this before any GPU session (S13+) or before pointing another task at
the machine.

## The machine

- **`backhouse`** — a Windows 11 desktop, **RTX 4070 Ti SUPER (16 GB VRAM)**, driver 591.86.
  *(Plan ruling #9 originally assumed an RTX 4080; the real card is a 4070 Ti SUPER — same Ada gen,
  same 16 GB, ~6% less bandwidth. S12 measured the real card.)*
- On Ben's **Tailscale** tailnet: `100.81.254.12`, Windows user `ben`.
- Real work runs inside **WSL2 Ubuntu 24.04** (Linux user `greff`): 28 cores, 31 GB RAM,
  ~950 GB free on `/`, CUDA 13.3 toolkit (`/usr/local/cuda`, `nvcc` arch `sm_89`), Rust 1.93.1,
  Python 3.13. `nvidia-smi` works inside WSL.

## Remote access (works from anywhere, not just the same wifi)

Tailscale is a WireGuard mesh VPN, so the box is reachable over **any** internet connection once
both machines are signed into the tailnet — same LAN, home, hotspot, wherever. It connects directly
when possible and falls back to an encrypted relay otherwise.

Physical caveats (the only things that break reachability): the PC must be **powered on, awake, and
online** — a desktop that **sleeps/hibernates won't answer** (set it to never sleep on AC), and
Tailscale must be running (it autostarts on boot).

- **From Ben's Mac:** `ssh backhouse` just works — the alias + key are in `~/.ssh/config`
  (`HostName 100.81.254.12`, `User ben`, `IdentityFile ~/.ssh/id_ed25519`).
- **From a different machine:** install Tailscale (same account) + add your public key to the PC's
  `C:\ProgramData\ssh\administrators_authorized_keys` (admin account → that system file, not
  `~/.ssh`, and it must be ACL-locked to SYSTEM+Administrators or sshd ignores it).

## The critical mechanic: SSH lands in `cmd.exe`, bounce into WSL

Windows OpenSSH drops you at `cmd.exe`, not bash. **Pipe bash scripts over stdin** — this sidesteps
all cmd.exe quoting/metachar problems (inline `"wsl bash -lc '... | grep ...'"` breaks: cmd eats the
pipe).

```bash
# run a bash script in WSL (default user greff)
ssh backhouse "wsl bash -l" < script.sh

# run as root (for apt; no Linux password needed — WSL root is free to the Windows owner)
ssh backhouse "wsl -u root bash -l" < script.sh

# copy files Mac→box (the --rsync-path bounces rsync into the WSL filesystem)
rsync -az -e ssh --rsync-path="wsl rsync" ./localdir backhouse:/home/greff/DEST/
```

`sudo` in WSL asks for a password nobody has — run root things via `wsl -u root` instead. **Do not**
write a `NOPASSWD` sudoers file (unrequested persistent privilege change; the harness guard blocks
it anyway).

## Long jobs survive only if you keep them alive

**WSL kills background / `nohup`'d jobs when the launching SSH session exits.** Two options:

- **Simple** — keep the job in the **foreground** of an SSH session you hold open (a background task
  on your side). Dies if your connection drops.
- **Robust (survives disconnects) — `tmux` inside WSL (session 28, measured):** the distro keeps
  running after the launching ssh session closes, and a tmux server escapes the job-object kill:
  ```bash
  # stage a script on the box (through the stdin-piped bash), then:
  ssh backhouse "wsl bash -l" <<< 'tmux new-session -d -s run "bash /home/greff/run_gpu.sh"'
  ssh backhouse "wsl bash -l" <<< 'tmux ls; tail -5 ~/gpu_run.log'     # poll
  ```
  A 40 s CUDA job launched this way outlived its ssh session; the engine harness runs this way.
  *(The earlier `powershell Start-Process wsl …` recipe is NOT reliable for scripts: an inline
  `bash -c 'sleep 100; echo …'` survived, but `bash -c /home/greff/script.sh` never executed under
  Start-Process in session 28, through both cmd-quoted and `-File`-launched forms — retired.)*

## CRUCIBLE on the box

- Repo lives at `~/inquiry-project` in WSL (kept in sync by rsync, or `git pull` if you set up
  GitHub auth). Full workspace builds in ~27 s incl. hdf5-from-source.
- The GPU crate `crates/gpu` (`crucible-gpu`) has its own empty `[workspace]` so the laptop
  `check.sh` never builds it. Build/run on the box:
  ```bash
  . ~/.cargo/env
  export PATH=/usr/local/cuda/bin:$PATH
  export LD_LIBRARY_PATH=/usr/local/cuda/lib64
  cd ~/inquiry-project/crates/gpu && cargo build --release && ./target/release/gpu_spike
  ```
- Standalone spike benchmarks: `crates/gpu/cuda/bench/` (`nvcc -O3 -arch=sm_89 <file>.cu && ./a.out`).

## Sharing the box with other tasks

One GPU, 16 GB VRAM — GPU work contends. If a non-CRUCIBLE task uses the box (e.g. NN training):
keep it out of `~/inquiry-project` (use `~/ml`, `~/nn`, …); `nvidia-smi` before a heavy run; and if
two GPU jobs would overlap, confirm with Ben first. Standard PyTorch/TF wheels bundle their own CUDA
runtime and only need the **driver** (present, WSL-exposed) — they do not need the system CUDA 13.3
toolkit (that's for compiling custom kernels).
