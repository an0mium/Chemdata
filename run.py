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
    data_dir = os.path.join(base_dir, 'data')
    ensure_directory(data_dir)

    print("\n=== Receptor Ligand Data Pipeline ===\n")

    try:
        # Step 1: Compile compounds
        print("Step 1: Compiling compound data...")
        from scripts.compile_compounds import main as compile_main
        await compile_main()
        print("\nCompound compilation complete!")

        # Step 2: Start web server
        print("\nStep 2: Starting web interface...")
        web_dir = os.path.join(base_dir, 'web')
        server_process = subprocess.Popen(
            [sys.executable, 'app.py'],
            cwd=web_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        # Wait for server to start
        time.sleep(2)
        if server_process.poll() is not None:
            # Server failed to start
            stdout, stderr = server_process.communicate()
            print("Error starting web server:")
            print(stderr.decode())
            sys.exit(1)

        # Open browser
        print("\nStarting browser...")
        webbrowser.open('http://localhost:5000')

        print("\nWeb interface running at: http://localhost:5000")
        print("\nPress Ctrl+C to stop the server")

        # Keep running until interrupted
        while True:
            await asyncio.sleep(1)

    except KeyboardInterrupt:
        print("\nShutting down...")
        if 'server_process' in locals():
            server_process.terminate()
            server_process.wait()
        print("Done!")


if __name__ == '__main__':
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\nExiting...")
        sys.exit(0)
