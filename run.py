#!/usr/bin/env python3
"""
Script to run both the compound compilation and web interface.
"""

import os
import sys
import asyncio
import subprocess
from pathlib import Path
import webbrowser
import time


def ensure_directory(path: str) -> None:
    """Ensure directory exists."""
    Path(path).mkdir(parents=True, exist_ok=True)


async def main():
    # Setup directories
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")
    ensure_directory(data_dir)

    print("\n=== Receptor Ligand Data Pipeline ===\n")

    try:
        # Step 1: Compile compounds
        print("Step 1: Compiling compound data...")
        try:
            from scripts.compile_compounds import main as compile_main

            await compile_main()
            print("\nCompound compilation complete!")
        except Exception as e:
            print(f"\nError during compound compilation: {str(e)}")
            print("Continuing with existing data...")

        # Step 2: Start web server
        print("\nStep 2: Starting web interface...")
        web_dir = os.path.join(base_dir, "web")
        print(f"Starting Flask server in {web_dir}...")

        # Start server with enhanced output
        server_process = subprocess.Popen(
            [sys.executable, "app.py"],
            cwd=web_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=1,
            universal_newlines=True,
        )

        # Wait for server to start and check output
        start_time = time.time()
        server_started = False
        while time.time() - start_time < 10:  # 10 second timeout
            if server_process.poll() is not None:
                out, err = server_process.communicate()
                print("Server failed to start!")
                print("Output:", out)
                print("Error:", err)
                sys.exit(1)

            # Check if server is responding
            try:
                import urllib.request

                urllib.request.urlopen("http://localhost:5001/")
                server_started = True
                break
            except (urllib.error.URLError, ConnectionRefusedError):
                await asyncio.sleep(0.5)

        if not server_started:
            print("Timeout waiting for server to start!")
            server_process.terminate()
            sys.exit(1)
        if server_process.poll() is not None:
            print("Error: Flask server failed to start")
            sys.exit(1)

        # Open browser
        print("\nStarting browser...")
        webbrowser.open("http://localhost:5001")

        print("\nWeb interface running at: http://localhost:5001")
        print("\nPress Ctrl+C to stop the server")

        # Keep running until interrupted
        while True:
            await asyncio.sleep(1)

    except KeyboardInterrupt:
        print("\nShutting down...")
        if "server_process" in locals():
            try:
                server_process.terminate()
                server_process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                server_process.kill()
            except Exception as e:
                print(f"Error during shutdown: {str(e)}")
        print("Done!")


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\nExiting...")
        sys.exit(0)
