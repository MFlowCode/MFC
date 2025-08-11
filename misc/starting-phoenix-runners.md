# Launching Phoenix Runners

The Phoenix runners were repeatedly failing due to a network error.
Spencer managed to fix it via [this PR](https://github.com/MFlowCode/MFC/pull/933) and by running things through a socks5 proxy on each login node that holds a runner.
These are documented for Spencer or his next of kin.

__The runners are started via the following process__

1. Log in to the login node <x> via `ssh login-phoenix-rh9-<x>.pace.gatech.edu`. `<x>` can be `1` through `6` on Phoenix.
  * Detour: Make sure no stray `ssh` daemons are sitting around: `pkill -9 sshd`.
  * You can probably keep your terminal alive via `fuser -k -9 ~/nohup.out`, which kills (signal 9) whatever process is writing to that no-hangup file (the daemon we care about)
2. Log back into the same login node because you may have just nuked your session
  * Detour: Make sure stray runners on that login node are dead (one liner): `pkill -9 -f -E 'run.sh|Runner.listener|Runner.helper'`
  * If cautious, check that no runner processes are left over. `top` followed by `u` and `<type your user name>` and return.
3. Execute from your home directory: `nohup ssh -N -D 1080 -vvv login-phoenix-rh9-<x>.pace.gatech.edu &`, replacing `<x>` with the login node number
  * This starts a proxy to tunnel a new ssh session through
4. Navigate to your runner's directory (or create a runner directory if you need).
  * Right now they are in Spencer's `scratch/mfc-runners/action-runner-<runner#>`
5. Run the alias `start_runner`, which dumps output `~/runner.out`
  * If one doesn't have this alias yet, create and source it in your `.bashrc` or similar:
```bash
alias start_runner=' \
  http_proxy="socks5://localhost:1080" \
  https_proxy="socks5://localhost:1080" \
  no_proxy="localhost,127.0.0.1,github.com,api.github.com,pipelines.actions.githubusercontent.com,alive.github.com,pypi.org,files.pythonhosted.org,fftw.org,www.fftw.org" \
  NO_PROXY="localhost,127.0.0.1,github.com,api.github.com,pipelines.actions.githubusercontent.com,alive.github.com,pypi.org,files.pythonhosted.org,fftw.org,www.fftw.org" \
  RUNNER_DEBUG=1 \
  ACTIONS_STEP_DEBUG=1 \
  GITHUB_ACTIONS_RUNNER_PREFER_IP_FAMILY=ipv4 \
  DOTNET_SYSTEM_NET_SOCKETS_KEEPALIVE_TIME=00:01:00 \
  DOTNET_SYSTEM_NET_SOCKETS_KEEPALIVE_INTERVAL=00:00:20 \
  DOTNET_SYSTEM_NET_SOCKETS_KEEPALIVE_RETRYCOUNT=5 \
  nohup ./run.sh > ~/runner.out 2>&1 &'
```
6. You're done


### For inquisitive minds 

__Why the `start_runner` alias?__

1. `alias start_runner='…'`  
   Defines a new shell alias named `start_runner`. Whenever you run `start_runner`, the shell will execute everything between the single quotes as if you’d typed it at the prompt.

2. `http_proxy="socks5://localhost:1080"`  
   Sets the `http_proxy` environment variable so that any HTTP traffic from the runner is sent through a SOCKS5 proxy listening on `localhost:1080`.

3. `https_proxy="socks5://localhost:1080"`  
   Tells HTTPS-aware tools to use that same local SOCKS5 proxy for HTTPS requests.

4. `no_proxy="localhost,127.0.0.1,github.com,api.github.com,pipelines.actions.githubusercontent.com,alive.github.com,pypi.org,files.pythonhosted.org,fftw.org,www.fftw.org"`  
   Lists hosts and domains that should bypass the proxy entirely. Commonly used for internal or high-volume endpoints where you don’t want proxy overhead.

5. `NO_PROXY="localhost,127.0.0.1,github.com,api.github.com,pipelines.actions.githubusercontent.com,alive.github.com,pypi.org,files.pythonhosted.org,fftw.org,www.fftw.org"`  
   Same list as `no_proxy`—some programs only check the uppercase `NO_PROXY` variable.

6. `RUNNER_DEBUG=1`  
   Enables debug-level logging in the GitHub Actions runner itself, so you’ll see more verbose internal messages in its logs.

7. `ACTIONS_STEP_DEBUG=1`  
   Turns on step-level debug logging for actions you invoke—handy if you need to trace exactly what each action is doing under the hood.

8. `GITHUB_ACTIONS_RUNNER_PREFER_IP_FAMILY=ipv4`  
   Forces the runner to resolve DNS names to IPv4 addresses only. Useful if your proxy or network has spotty IPv6 support.

9. `DOTNET_SYSTEM_NET_SOCKETS_KEEPALIVE_TIME=00:01:00`  
   For .NET–based tasks: sets the initial TCP keepalive timeout to 1 minute (after 1 minute of idle, a keepalive probe is sent).

10. `DOTNET_SYSTEM_NET_SOCKETS_KEEPALIVE_INTERVAL=00:00:20`  
    If the first keepalive probe gets no response, wait 20 seconds between subsequent probes.

11. `DOTNET_SYSTEM_NET_SOCKETS_KEEPALIVE_RETRYCOUNT=5`  
    If probes continue to go unanswered, retry up to 5 times before declaring the connection dead.

12. `nohup ./run.sh > ~/runner.out 2>&1 &`  
    - `nohup … &` runs `./run.sh` in the background and makes it immune to hangups (so it keeps running if you log out).  
    - `> ~/runner.out` redirects **stdout** to the file `runner.out` in your home directory.  
    - `2>&1` redirects **stderr** into the same file, so you get a combined log of everything the script prints.

__Why the extra ssh command?__

1. `http_proxy="socks5://localhost:1080"`  
   Routes all HTTP traffic through a local SOCKS5 proxy on port 1080.

2. `https_proxy="socks5://localhost:1080"`  
   Routes all HTTPS traffic through the same proxy.

3. `no_proxy="localhost,127.0.0.1,github.com,api.github.com,pipelines.actions.githubusercontent.com,alive.github.com,pypi.org,files.pythonhosted.org,fftw.org,www.fftw.org"`  
   Specifies hosts and domains that bypass the proxy entirely. Includes specific things that MFC's CMake will try to `wget` (e.g., `fftw`) or some other non `git` command. Allows `git clone` to work.

4. `NO_PROXY="localhost,127.0.0.1,github.com,api.github.com,pipelines.actions.githubusercontent.com,alive.github.com,pypi.org,files.pythonhosted.org,fftw.org,www.fftw.org"`  
   Same bypass list for applications that only check the uppercase variable.

5. `RUNNER_DEBUG=1`  
   Enables verbose internal logging in the GitHub Actions runner.

6. `GITHUB_ACTIONS_RUNNER_PREFER_IP_FAMILY=ipv4`  
   Forces DNS resolution to IPv4 to avoid IPv6 issues.

7. `DOTNET_SYSTEM_NET_SOCKETS_KEEPALIVE_TIME=00:01:00`  
   (For .NET tasks) sends the first TCP keepalive probe after 1 minute of idle.

8. `DOTNET_SYSTEM_NET_SOCKETS_KEEPALIVE_INTERVAL=00:00:20`  
   Waits 20 seconds between subsequent TCP keepalive probes.

9. `DOTNET_SYSTEM_NET_SOCKETS_KEEPALIVE_RETRYCOUNT=5`  
   Retries keepalive probes up to 5 times before closing the connection.

10. `nohup ./run.sh > ~/runner.out 2>&1 &`  
    Runs `run.sh` in the background, immune to hangups, redirecting both stdout and stderr to `~/runner.out`.
