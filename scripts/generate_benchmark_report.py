#!/usr/bin/env python3
"""Generate HTML report from benchmark results."""

import argparse
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from jinja2 import Environment, FileSystemLoader

# Configure plotting style
plt.style.use("seaborn")
sns.set_palette("husl")


class BenchmarkReport:
    """Generate HTML report from benchmark results."""

    def __init__(self, output_path: str):
        """Initialize report generator.

        Args:
            output_path: Path to output HTML file
        """
        self.output_path = output_path
        self.plots_dir = Path(output_path).parent / "plots"
        self.plots_dir.mkdir(exist_ok=True)

        # Load templates
        self.env = Environment(
            loader=FileSystemLoader("web/templates"),
            trim_blocks=True,
            lstrip_blocks=True,
        )

    def parse_deps_results(self, deps_file: str) -> Dict:
        """Parse dependency benchmark results.

        Args:
            deps_file: Path to dependency results file

        Returns:
            Dictionary containing parsed results
        """
        results = {}
        with open(deps_file) as f:
            lines = f.readlines()

        # Parse pip times
        pip_cold = float(lines[0].split(": ")[1].split()[0])
        pip_warm = float(lines[1].split(": ")[1].split()[0])
        pip_uninstall = float(lines[2].split(": ")[1].split()[0])

        # Parse uv times
        uv_cold = float(lines[4].split(": ")[1].split()[0])
        uv_warm = float(lines[5].split(": ")[1].split()[0])
        uv_uninstall = float(lines[6].split(": ")[1].split()[0])

        # Parse speedups
        cold_speedup = float(lines[9].split(": ")[1].strip("x"))
        warm_speedup = float(lines[10].split(": ")[1].strip("x"))
        uninstall_speedup = float(lines[11].split(": ")[1].strip("x"))

        results["pip"] = {
            "cold_install": pip_cold,
            "warm_install": pip_warm,
            "uninstall": pip_uninstall,
        }
        results["uv"] = {
            "cold_install": uv_cold,
            "warm_install": uv_warm,
            "uninstall": uv_uninstall,
        }
        results["speedups"] = {
            "cold_install": cold_speedup,
            "warm_install": warm_speedup,
            "uninstall": uninstall_speedup,
        }

        # Generate dependency comparison plot
        self.plot_deps_comparison(results)

        return results

    def parse_data_results(self, data_dir: str) -> Dict:
        """Parse data processing benchmark results.

        Args:
            data_dir: Path to data results directory

        Returns:
            Dictionary containing parsed results
        """
        results = {}

        # Parse load time
        with open(os.path.join(data_dir, "load_time.txt")) as f:
            for line in f:
                if "User time" in line:
                    results["load_time"] = float(line.split(": ")[1])
                elif "Maximum resident set size" in line:
                    # Convert KB to MB
                    results["load_memory"] = int(line.split(": ")[1]) / 1024

        # Parse process time
        with open(os.path.join(data_dir, "process_time.txt")) as f:
            for line in f:
                if "User time" in line:
                    results["process_time"] = float(line.split(": ")[1])
                elif "Maximum resident set size" in line:
                    # Convert KB to MB
                    results["process_memory"] = int(line.split(": ")[1]) / 1024

        # Generate data processing plot
        self.plot_data_performance(results)

        return results

    def parse_ml_results(self, ml_dir: str) -> Dict:
        """Parse ML benchmark results.

        Args:
            ml_dir: Path to ML results directory

        Returns:
            Dictionary containing parsed results
        """
        results = {}

        # Parse load time
        with open(os.path.join(ml_dir, "load_time.txt")) as f:
            for line in f:
                if "User time" in line:
                    results["load_time"] = float(line.split(": ")[1])
                elif "Maximum resident set size" in line:
                    # Convert KB to MB
                    results["load_memory"] = int(line.split(": ")[1]) / 1024

        # Parse inference time
        with open(os.path.join(ml_dir, "inference_time.txt")) as f:
            for line in f:
                if "User time" in line:
                    # Calculate per-inference time
                    results["inference_time"] = float(line.split(": ")[1]) / 100
                elif "Maximum resident set size" in line:
                    # Convert KB to MB
                    results["inference_memory"] = int(line.split(": ")[1]) / 1024

        # Parse batch inference time
        with open(os.path.join(ml_dir, "batch_time.txt")) as f:
            for line in f:
                if "User time" in line:
                    # Calculate per-sample time
                    results["batch_time"] = float(line.split(": ")[1]) / 100
                elif "Maximum resident set size" in line:
                    # Convert KB to MB
                    results["batch_memory"] = int(line.split(": ")[1]) / 1024

        # Parse GPU utilization if available
        gpu_util_file = os.path.join(ml_dir, "gpu_util.txt")
        if os.path.exists(gpu_util_file):
            with open(gpu_util_file) as f:
                utils = []
                for line in f:
                    if not line.startswith("#"):
                        utils.append(int(line.split()[1]))
                results["gpu_utilization"] = sum(utils) / len(utils)

        # Generate ML performance plots
        self.plot_ml_performance(results)

        return results

    def parse_web_results(self, web_dir: str) -> Dict:
        """Parse web benchmark results.

        Args:
            web_dir: Path to web results directory

        Returns:
            Dictionary containing parsed results
        """
        results = {}

        # Parse Apache Bench results
        with open(os.path.join(web_dir, "ab_results.txt")) as f:
            for line in f:
                if "Requests per second" in line:
                    rps = line.split(":")[1].split("[")[0]
                    results["requests_per_second"] = float(rps)
                elif "Time per request" in line and "across all" not in line:
                    tpr = line.split(":")[1].split("[")[0]
                    results["time_per_request"] = float(tpr)
                elif "Transfer rate" in line:
                    tr = line.split(":")[1].split("[")[0]
                    results["transfer_rate"] = float(tr)

        # Parse API results
        with open(os.path.join(web_dir, "api_results.txt")) as f:
            for line in f:
                if "Requests per second" in line:
                    rps = line.split(":")[1].split("[")[0]
                    results["api_requests_per_second"] = float(rps)
                elif "Time per request" in line and "across all" not in line:
                    tpr = line.split(":")[1].split("[")[0]
                    results["api_time_per_request"] = float(tpr)

        # Parse WebSocket results if available
        ws_file = os.path.join(web_dir, "ws_results.txt")
        if os.path.exists(ws_file):
            with open(ws_file) as f:
                data = json.load(f)
                results["ws_latency"] = data["average_latency"]
                results["ws_throughput"] = data["messages_per_second"]

        # Generate web performance plots
        self.plot_web_performance(results)

        return results

    def parse_memory_results(self, memory_dir: str) -> Dict:
        """Parse memory benchmark results.

        Args:
            memory_dir: Path to memory results directory

        Returns:
            Dictionary containing parsed results
        """
        results = {}

        # Parse data processing memory
        with open(os.path.join(memory_dir, "data_memory.txt")) as f:
            for line in f:
                if "Maximum resident set size" in line:
                    # Convert KB to MB
                    results["data_memory"] = int(line.split(": ")[1]) / 1024

        # Parse ML inference memory
        with open(os.path.join(memory_dir, "ml_memory.txt")) as f:
            for line in f:
                if "Maximum resident set size" in line:
                    # Convert KB to MB
                    results["ml_memory"] = int(line.split(": ")[1]) / 1024

        # Parse memory leaks
        with open(os.path.join(memory_dir, "leaks.txt")) as f:
            leaks = []
            current_leak = {}
            for line in f:
                if "Memory leak detected" in line:
                    if current_leak:
                        leaks.append(current_leak)
                    current_leak = {"size": int(line.split(": ")[1].split()[0])}
                elif "Location" in line and current_leak:
                    current_leak["location"] = line.split(": ")[1].strip()
            if current_leak:
                leaks.append(current_leak)
            results["memory_leaks"] = leaks

        # Generate memory usage plots
        self.plot_memory_usage(results)

        return results

    def parse_gpu_results(self, gpu_dir: str) -> Optional[Dict]:
        """Parse GPU benchmark results.

        Args:
            gpu_dir: Path to GPU results directory

        Returns:
            Dictionary containing parsed results, or None if no GPU
        """
        if not os.path.exists(gpu_dir):
            return None

        results = {}

        # Parse memory usage
        with open(os.path.join(gpu_dir, "memory.txt")) as f:
            memory_usage = []
            for line in f:
                if not line.startswith("#"):
                    memory_usage.append(int(line.split()[1]))
            results["memory_usage"] = sum(memory_usage) / len(memory_usage)

        # Parse utilization
        with open(os.path.join(gpu_dir, "utilization.txt")) as f:
            utilization = []
            for line in f:
                if not line.startswith("#"):
                    utilization.append(int(line.split()[1]))
            results["utilization"] = sum(utilization) / len(utilization)

        # Parse multi-GPU results if available
        multi_gpu_file = os.path.join(gpu_dir, "multi_gpu.txt")
        if os.path.exists(multi_gpu_file):
            with open(multi_gpu_file) as f:
                data = json.load(f)
                results["multi_gpu_speedup"] = data["speedup"]
                results["gpu_scaling_efficiency"] = data["scaling_efficiency"]

        # Generate GPU performance plots
        self.plot_gpu_performance(results)

        return results

    def plot_deps_comparison(self, results: Dict) -> None:
        """Generate dependency comparison plots.

        Args:
            results: Parsed dependency results
        """
        # Bar plot comparing pip vs uv
        fig, ax = plt.subplots(figsize=(10, 6))

        x = np.arange(3)
        width = 0.35

        pip_times = [
            results["pip"]["cold_install"],
            results["pip"]["warm_install"],
            results["pip"]["uninstall"],
        ]
        uv_times = [
            results["uv"]["cold_install"],
            results["uv"]["warm_install"],
            results["uv"]["uninstall"],
        ]

        ax.bar(x - width / 2, pip_times, width, label="pip")
        ax.bar(x + width / 2, uv_times, width, label="uv")

        ax.set_ylabel("Time (seconds)")
        ax.set_title("Package Manager Performance")
        ax.set_xticks(x)
        ax.set_xticklabels(["Cold Install", "Warm Install", "Uninstall"])
        ax.legend()

        plt.savefig(self.plots_dir / "deps_comparison.png")
        plt.close()

    def plot_data_performance(self, results: Dict) -> None:
        """Generate data processing performance plots.

        Args:
            results: Parsed data processing results
        """
        # Bar plot of processing times
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        times = [results["load_time"], results["process_time"]]
        ax1.bar(["Load", "Process"], times)
        ax1.set_ylabel("Time (seconds)")
        ax1.set_title("Data Processing Time")

        memory = [results["load_memory"], results["process_memory"]]
        ax2.bar(["Load", "Process"], memory)
        ax2.set_ylabel("Memory (MB)")
        ax2.set_title("Memory Usage")

        plt.tight_layout()
        plt.savefig(self.plots_dir / "data_performance.png")
        plt.close()

    def plot_ml_performance(self, results: Dict) -> None:
        """Generate ML performance plots.

        Args:
            results: Parsed ML results
        """
        # Bar plot of inference times
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        times = [results["inference_time"], results["batch_time"]]
        ax1.bar(["Single", "Batch"], times)
        ax1.set_ylabel("Time per sample (seconds)")
        ax1.set_title("Inference Time")

        memory = [results["inference_memory"], results["batch_memory"]]
        ax2.bar(["Single", "Batch"], memory)
        ax2.set_ylabel("Memory (MB)")
        ax2.set_title("Memory Usage")

        plt.tight_layout()
        plt.savefig(self.plots_dir / "ml_performance.png")
        plt.close()

        # GPU utilization if available
        if "gpu_utilization" in results:
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.bar(["GPU"], [results["gpu_utilization"]])
            ax.set_ylabel("Utilization (%)")
            ax.set_title("GPU Utilization")
            plt.savefig(self.plots_dir / "gpu_utilization.png")
            plt.close()

    def plot_web_performance(self, results: Dict) -> None:
        """Generate web performance plots.

        Args:
            results: Parsed web results
        """
        # Bar plot of request rates
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        rates = [results["requests_per_second"], results["api_requests_per_second"]]
        ax1.bar(["Web", "API"], rates)
        ax1.set_ylabel("Requests per second")
        ax1.set_title("Request Rate")

        times = [results["time_per_request"], results["api_time_per_request"]]
        ax2.bar(["Web", "API"], times)
        ax2.set_ylabel("Time (ms)")
        ax2.set_title("Response Time")

        plt.tight_layout()
        plt.savefig(self.plots_dir / "web_performance.png")
        plt.close()

        # WebSocket performance if available
        if "ws_latency" in results:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

            ax1.bar(["WebSocket"], [results["ws_latency"]])
            ax1.set_ylabel("Latency (ms)")
            ax1.set_title("WebSocket Latency")

            ax2.bar(["WebSocket"], [results["ws_throughput"]])
            ax2.set_ylabel("Messages per second")
            ax2.set_title("WebSocket Throughput")

            plt.tight_layout()
            plt.savefig(self.plots_dir / "websocket_performance.png")
            plt.close()

    def plot_memory_usage(self, results: Dict) -> None:
        """Generate memory usage plots.

        Args:
            results: Parsed memory results
        """
        # Bar plot of memory usage
        fig, ax = plt.subplots(figsize=(8, 5))

        memory = [results["data_memory"], results["ml_memory"]]
        ax.bar(["Data Processing", "ML Inference"], memory)
        ax.set_ylabel("Memory (MB)")
        ax.set_title("Peak Memory Usage")

        plt.savefig(self.plots_dir / "memory_usage.png")
        plt.close()

        # Memory leak summary if any found
        if results["memory_leaks"]:
            fig, ax = plt.subplots(figsize=(10, 6))

            leaks = results["memory_leaks"]
            locations = [leak["location"] for leak in leaks]
            sizes = [leak["size"] for leak in leaks]

            ax.barh(locations, sizes)
            ax.set_xlabel("Leak Size (bytes)")
            ax.set_title("Memory Leaks")

            plt.tight_layout()
            plt.savefig(self.plots_dir / "memory_leaks.png")
            plt.close()

    def plot_gpu_performance(self, results: Dict) -> None:
        """Generate GPU performance plots.

        Args:
            results: Parsed GPU results
        """
        if not results:
            return

        # Bar plot of GPU metrics
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        ax1.bar(["Memory"], [results["memory_usage"]])
        ax1.set_ylabel("Memory Usage (%)")
        ax1.set_title("GPU Memory Usage")

        ax2.bar(["Utilization"], [results["utilization"]])
        ax2.set_ylabel("Utilization (%)")
        ax2.set_title("GPU Utilization")

        plt.tight_layout()
        plt.savefig(self.plots_dir / "gpu_performance.png")
        plt.close()

        # Multi-GPU scaling if available
        if "multi_gpu_speedup" in results:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

            ax1.bar(["Speedup"], [results["multi_gpu_speedup"]])
            ax1.set_ylabel("Speedup Factor")
            ax1.set_title("Multi-GPU Speedup")

            ax2.bar(["Efficiency"], [results["gpu_scaling_efficiency"]])
            ax2.set_ylabel("Efficiency (%)")
            ax2.set_title("GPU Scaling Efficiency")

            plt.tight_layout()
            plt.savefig(self.plots_dir / "multi_gpu_performance.png")
            plt.close()

    def generate_report(
        self,
        deps_file: str,
        data_dir: str,
        ml_dir: str,
        web_dir: str,
        memory_dir: str,
        gpu_dir: str,
    ) -> None:
        """Generate HTML report from benchmark results.

        Args:
            deps_file: Path to dependency results file
            data_dir: Path to data results directory
            ml_dir: Path to ML results directory
            web_dir: Path to web results directory
            memory_dir: Path to memory results directory
            gpu_dir: Path to GPU results directory
        """
        # Parse results
        deps_results = self.parse_deps_results(deps_file)
        data_results = self.parse_data_results(data_dir)
        ml_results = self.parse_ml_results(ml_dir)
        web_results = self.parse_web_results(web_dir)
        memory_results = self.parse_memory_results(memory_dir)
        gpu_results = self.parse_gpu_results(gpu_dir)

        # Load template
        template = self.env.get_template("benchmark_report.html")

        # Get relative path to plots directory
        plots_path = os.path.relpath(self.plots_dir, os.path.dirname(self.output_path))

        # Render template
        html = template.render(
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            deps_results=deps_results,
            data_results=data_results,
            ml_results=ml_results,
            web_results=web_results,
            memory_results=memory_results,
            gpu_results=gpu_results,
            plots_dir=plots_path,
        )

        # Write report
        with open(self.output_path, "w") as f:
            f.write(html)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Generate benchmark report")
    parser.add_argument("--deps", required=True, help="Path to dependency results")
    parser.add_argument("--data", required=True, help="Path to data results dir")
    parser.add_argument("--ml", required=True, help="Path to ML results dir")
    parser.add_argument("--web", required=True, help="Path to web results dir")
    parser.add_argument("--memory", required=True, help="Path to memory results dir")
    parser.add_argument("--gpu", required=True, help="Path to GPU results dir")
    parser.add_argument("--output", required=True, help="Path to output HTML file")

    args = parser.parse_args()

    report = BenchmarkReport(args.output)
    report.generate_report(
        args.deps,
        args.data,
        args.ml,
        args.web,
        args.memory,
        args.gpu,
    )


if __name__ == "__main__":
    main()
