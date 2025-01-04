#!/usr/bin/env python3
"""
Main entry point for the receptor ligand data pipeline.
Handles both data processing and web interface.
"""

import os
import sys
import asyncio
import subprocess
from pathlib import Path
import webbrowser
import time
from datetime import datetime
from tqdm import tqdm
from typing import Optional

from binding_data_processor import BindingDataProcessor
from api_client import PubMedClient


def ensure_directory(path: str) -> None:
    """Ensure directory exists."""
    Path(path).mkdir(parents=True, exist_ok=True)


async def process_bindingdb_data(input_file: str, output_file: str) -> Optional[str]:
    """
    Process BindingDB TSV file into web-friendly format.

    Args:
        input_file: Path to input TSV file
        output_file: Path to output TSV file

    Returns:
        Path to output file if successful, None otherwise
    """
    print("\nStep 1: Processing BindingDB data...")

    try:
        processor = BindingDataProcessor(PubMedClient())

        # Count total lines first for progress bar
        print("Counting entries in BindingDB...")
        total_lines = sum(1 for _ in open(input_file, "r", encoding="utf-8"))
        print(f"Found {total_lines:,} entries")

        # Process with progress bar
        with tqdm(total=total_lines, desc="Processing compounds") as pbar:
            processor.process_bindingdb_file(
                input_file, output_file, progress_callback=lambda: pbar.update(1)
            )

        print("\nBindingDB processing complete!")
        print(f"Output saved to: {output_file}")
        return output_file

    except Exception as e:
        print(f"\nError during BindingDB processing: {str(e)}")
        return None


async def start_web_interface(data_file: str) -> subprocess.Popen:
    """
    Start the web interface.

    Args:
        data_file: Path to data file for web app

    Returns:
        Server process handle
    """
    print("\nStep 2: Starting web interface...")

    # Get web directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    web_dir = os.path.join(base_dir, "web")
    print(f"Starting Flask server in {web_dir}...")

    # Set data file environment variable
    env = os.environ.copy()
    env["CHEMDATA_FILE"] = data_file

    # Start server process
    server_process = subprocess.Popen(
        [sys.executable, "app.py"],
        cwd=web_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1,
        universal_newlines=True,
    )

    # Wait for server to start
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

    # Open browser
    print("\nStarting browser...")
    webbrowser.open("http://localhost:5001")

    print("\nWeb interface running at: http://localhost:5001")
    print("\nPress Ctrl+C to stop the server")

    return server_process


async def main():
    """Main entry point."""
    # Setup directories
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")
    ensure_directory(data_dir)

    print("\n=== Receptor Ligand Data Pipeline ===\n")

    try:
        # Process BindingDB data
        input_file = os.path.join(data_dir, "BindingDB_All.tsv")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = os.path.join(data_dir, f"receptor_compounds_{timestamp}.tsv")

        processed_file = await process_bindingdb_data(input_file, output_file)
        if not processed_file:
            # Try to use existing data file
            existing_files = sorted(
                [
                    f
                    for f in os.listdir(data_dir)
                    if f.startswith("receptor_compounds_")
                ],
                reverse=True,
            )
            if existing_files:
                processed_file = os.path.join(data_dir, existing_files[0])
                print(f"\nUsing existing data file: {processed_file}")
            else:
                print("\nNo existing data file found!")
                sys.exit(1)

        # Start web interface
        server_process = await start_web_interface(processed_file)

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
