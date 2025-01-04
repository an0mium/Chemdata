"""Test runner for structure processor tests.

Provides:
1. Test discovery and execution
2. Test result reporting
3. Coverage reporting
4. Performance profiling
5. HTML test report generation
"""

import os
import sys
import time
import unittest
import coverage
import cProfile
import pstats
from pathlib import Path
from typing import List, Optional, Set, Tuple
from unittest.runner import TextTestResult
from unittest.suite import TestSuite

from .test_config import test_config, logger


class StructureTestResult(TextTestResult):
    """Enhanced test result with timing and profiling."""

    def __init__(self, *args, **kwargs):
        """Initialize test result."""
        super().__init__(*args, **kwargs)
        self.test_timings = {}
        self.test_profiles = {}
        self.skipped_tests = set()
        self.failed_tests = set()
        self.error_tests = set()

    def startTest(self, test):
        """Start timing test."""
        self._started_at = time.time()
        super().startTest(test)

    def addSuccess(self, test):
        """Record successful test timing."""
        elapsed = time.time() - self._started_at
        name = self.getDescription(test)
        self.test_timings[name] = elapsed
        super().addSuccess(test)

    def addError(self, test, err):
        """Record test error."""
        name = self.getDescription(test)
        self.error_tests.add(name)
        super().addError(test, err)

    def addFailure(self, test, err):
        """Record test failure."""
        name = self.getDescription(test)
        self.failed_tests.add(name)
        super().addFailure(test, err)

    def addSkip(self, test, reason):
        """Record skipped test."""
        name = self.getDescription(test)
        self.skipped_tests.add(name)
        super().addSkip(test, reason)

    def add_profile(self, test_name: str, profile_stats: pstats.Stats):
        """Add profiling data for test."""
        self.test_profiles[test_name] = profile_stats


class StructureTestRunner:
    """Test runner with enhanced reporting."""

    def __init__(
        self,
        verbosity: int = 2,
        failfast: bool = False,
        profile: bool = False,
        coverage: bool = True,
        html_report: bool = True,
    ):
        """Initialize test runner.

        Args:
            verbosity: Output verbosity level
            failfast: Stop on first failure
            profile: Enable performance profiling
            coverage: Enable coverage reporting
            html_report: Generate HTML test report
        """
        self.verbosity = verbosity
        self.failfast = failfast
        self.profile = profile
        self.enable_coverage = coverage
        self.html_report = html_report
        self.coverage = None
        self.test_dir = Path(__file__).parent

    def run(self, pattern: str = "test_*.py") -> StructureTestResult:
        """Run all tests matching pattern.

        Args:
            pattern: Test file pattern to match

        Returns:
            Test result object
        """
        # Start coverage if enabled
        if self.enable_coverage:
            self.coverage = coverage.Coverage()
            self.coverage.start()

        # Discover and load tests
        loader = unittest.TestLoader()
        suite = loader.discover(str(self.test_dir), pattern=pattern)

        # Create result object
        stream = sys.stderr
        descriptions = True
        result = StructureTestResult(stream, descriptions, self.verbosity)

        # Run tests
        if self.profile:
            profiler = cProfile.Profile()
            profiler.enable()

        suite.run(result)

        if self.profile:
            profiler.disable()
            stats = pstats.Stats(profiler)
            stats.sort_stats("cumulative")
            result.add_profile("full_suite", stats)

        # Stop coverage
        if self.enable_coverage:
            self.coverage.stop()
            self.coverage.save()

        # Generate reports
        self._generate_reports(result)

        return result

    def _generate_reports(self, result: StructureTestResult):
        """Generate test reports."""
        # Print test summary
        self._print_summary(result)

        # Generate coverage report if enabled
        if self.enable_coverage:
            self._generate_coverage_report()

        # Generate HTML report if enabled
        if self.html_report:
            self._generate_html_report(result)

    def _print_summary(self, result: StructureTestResult):
        """Print test result summary."""
        logger.info("\nTest Summary:")
        logger.info("-" * 60)

        # Overall stats
        total = result.testsRun
        passed = total - len(result.failures) - len(result.errors) - len(result.skipped)
        logger.info(f"Total tests: {total}")
        logger.info(f"Passed: {passed}")
        logger.info(f"Failed: {len(result.failures)}")
        logger.info(f"Errors: {len(result.errors)}")
        logger.info(f"Skipped: {len(result.skipped)}")

        # Timing information
        if result.test_timings:
            logger.info("\nTest Timings:")
            for name, elapsed in sorted(
                result.test_timings.items(),
                key=lambda x: x[1],
                reverse=True,
            )[:10]:
                logger.info(f"{name}: {elapsed:.3f}s")

        # Profile information
        if result.test_profiles:
            logger.info("\nProfile Data:")
            for name, stats in result.test_profiles.items():
                logger.info(f"\nProfile for {name}:")
                stats.print_stats(10)

        # Failed tests
        if result.failures:
            logger.info("\nFailures:")
            for test, trace in result.failures:
                logger.info(f"\n{test}")
                logger.info(trace)

        # Errors
        if result.errors:
            logger.info("\nErrors:")
            for test, trace in result.errors:
                logger.info(f"\n{test}")
                logger.info(trace)

    def _generate_coverage_report(self):
        """Generate coverage report."""
        if self.coverage:
            # Generate reports
            self.coverage.report()
            self.coverage.html_report(
                directory=str(test_config.get_output_path("coverage"))
            )

    def _generate_html_report(self, result: StructureTestResult):
        """Generate HTML test report."""
        report_dir = test_config.get_output_path("test_report")
        report_dir.mkdir(exist_ok=True)

        # Basic HTML report template
        html = """
        <html>
        <head>
            <title>Structure Tests Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                .summary { margin: 20px 0; }
                .passed { color: green; }
                .failed { color: red; }
                .error { color: darkred; }
                .skipped { color: orange; }
                .timing { margin: 20px 0; }
                pre { background: #f5f5f5; padding: 10px; }
            </style>
        </head>
        <body>
            <h1>Structure Tests Report</h1>
            
            <div class="summary">
                <h2>Summary</h2>
                <p>Total tests: {total}</p>
                <p class="passed">Passed: {passed}</p>
                <p class="failed">Failed: {failed}</p>
                <p class="error">Errors: {errors}</p>
                <p class="skipped">Skipped: {skipped}</p>
            </div>

            <div class="timing">
                <h2>Top 10 Slowest Tests</h2>
                <pre>{timings}</pre>
            </div>

            <div class="failures">
                <h2>Failures</h2>
                <pre>{failures}</pre>
            </div>

            <div class="errors">
                <h2>Errors</h2>
                <pre>{errors_text}</pre>
            </div>
        </body>
        </html>
        """

        # Format timing data
        timings = "\n".join(
            f"{name}: {elapsed:.3f}s"
            for name, elapsed in sorted(
                result.test_timings.items(),
                key=lambda x: x[1],
                reverse=True,
            )[:10]
        )

        # Format failure data
        failures = "\n\n".join(f"{test}\n{trace}" for test, trace in result.failures)

        # Format error data
        errors_text = "\n\n".join(f"{test}\n{trace}" for test, trace in result.errors)

        # Generate report
        total = result.testsRun
        passed = total - len(result.failures) - len(result.errors) - len(result.skipped)
        report = html.format(
            total=total,
            passed=passed,
            failed=len(result.failures),
            errors=len(result.errors),
            skipped=len(result.skipped),
            timings=timings,
            failures=failures,
            errors_text=errors_text,
        )

        # Write report
        report_path = report_dir / "index.html"
        report_path.write_text(report)
        logger.info(f"\nHTML report generated at {report_path}")


# Run tests if executed directly
if __name__ == "__main__":
    runner = StructureTestRunner(
        verbosity=2,
        failfast=False,
        profile=True,
        coverage=True,
        html_report=True,
    )
    result = runner.run()
    sys.exit(not result.wasSuccessful())
