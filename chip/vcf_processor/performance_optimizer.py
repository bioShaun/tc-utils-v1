"""Performance optimization utilities for VCF processor."""

import gc
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import psutil
import os

from loguru import logger


class MemoryMonitor:
    """Monitor memory usage during VCF processing."""
    
    def __init__(self, enable_monitoring: bool = True):
        """Initialize memory monitor.
        
        Args:
            enable_monitoring: Whether to enable memory monitoring
        """
        self.enable_monitoring = enable_monitoring
        self.process = psutil.Process(os.getpid()) if enable_monitoring else None
        self.peak_memory = 0
        self.memory_samples = []
        self.start_memory = 0
    
    def start(self) -> None:
        """Start memory monitoring."""
        if self.enable_monitoring and self.process:
            self.start_memory = self.process.memory_info().rss / (1024 * 1024)  # MB
            self.peak_memory = self.start_memory
            logger.debug(f"Memory monitoring started: {self.start_memory:.1f} MB")
    
    def sample(self, label: str = "") -> float:
        """Take a memory sample.
        
        Args:
            label: Optional label for the sample
            
        Returns:
            Current memory usage in MB
        """
        if not self.enable_monitoring or not self.process:
            return 0.0
        
        current_memory = self.process.memory_info().rss / (1024 * 1024)  # MB
        self.peak_memory = max(self.peak_memory, current_memory)
        
        sample_data = {
            "timestamp": time.time(),
            "memory_mb": current_memory,
            "label": label
        }
        self.memory_samples.append(sample_data)
        
        if label:
            logger.debug(f"Memory sample [{label}]: {current_memory:.1f} MB")
        
        return current_memory
    
    def get_peak_memory(self) -> float:
        """Get peak memory usage in MB."""
        return self.peak_memory
    
    def get_memory_increase(self) -> float:
        """Get memory increase from start in MB."""
        if not self.enable_monitoring:
            return 0.0
        return self.peak_memory - self.start_memory
    
    def force_gc(self) -> None:
        """Force garbage collection and sample memory."""
        if self.enable_monitoring:
            before_gc = self.sample("before_gc")
            gc.collect()
            after_gc = self.sample("after_gc")
            
            if before_gc > 0 and after_gc > 0:
                freed = before_gc - after_gc
                if freed > 1:  # Only log if significant memory was freed
                    logger.debug(f"Garbage collection freed {freed:.1f} MB")


class PerformanceProfiler:
    """Profile performance of VCF processing operations."""
    
    def __init__(self, enable_profiling: bool = True):
        """Initialize performance profiler.
        
        Args:
            enable_profiling: Whether to enable performance profiling
        """
        self.enable_profiling = enable_profiling
        self.timings = {}
        self.counters = {}
        self.memory_monitor = MemoryMonitor(enable_profiling)
    
    def start_timer(self, operation: str) -> None:
        """Start timing an operation.
        
        Args:
            operation: Name of the operation
        """
        if self.enable_profiling:
            self.timings[operation] = {"start": time.time(), "end": None, "duration": None}
    
    def end_timer(self, operation: str) -> float:
        """End timing an operation.
        
        Args:
            operation: Name of the operation
            
        Returns:
            Duration in seconds
        """
        if not self.enable_profiling or operation not in self.timings:
            return 0.0
        
        end_time = time.time()
        self.timings[operation]["end"] = end_time
        duration = end_time - self.timings[operation]["start"]
        self.timings[operation]["duration"] = duration
        
        return duration
    
    def increment_counter(self, counter: str, amount: int = 1) -> None:
        """Increment a performance counter.
        
        Args:
            counter: Name of the counter
            amount: Amount to increment by
        """
        if self.enable_profiling:
            self.counters[counter] = self.counters.get(counter, 0) + amount
    
    def get_timing_summary(self) -> Dict[str, float]:
        """Get summary of all timings.
        
        Returns:
            Dictionary of operation names to durations
        """
        summary = {}
        for operation, timing_data in self.timings.items():
            if timing_data["duration"] is not None:
                summary[operation] = timing_data["duration"]
        return summary
    
    def get_performance_report(self) -> str:
        """Generate a performance report.
        
        Returns:
            Formatted performance report
        """
        if not self.enable_profiling:
            return "Performance profiling disabled"
        
        lines = []
        lines.append("PERFORMANCE REPORT")
        lines.append("=" * 50)
        
        # Timing information
        if self.timings:
            lines.append("\nTIMINGS:")
            for operation, timing_data in self.timings.items():
                if timing_data["duration"] is not None:
                    lines.append(f"  {operation}: {timing_data['duration']:.3f}s")
        
        # Counter information
        if self.counters:
            lines.append("\nCOUNTERS:")
            for counter, value in self.counters.items():
                lines.append(f"  {counter}: {value:,}")
        
        # Memory information
        peak_memory = self.memory_monitor.get_peak_memory()
        memory_increase = self.memory_monitor.get_memory_increase()
        
        if peak_memory > 0:
            lines.append("\nMEMORY:")
            lines.append(f"  Peak usage: {peak_memory:.1f} MB")
            lines.append(f"  Memory increase: {memory_increase:.1f} MB")
        
        # Performance metrics
        if self.timings and self.counters:
            lines.append("\nPERFORMANCE METRICS:")
            
            # Calculate rates
            total_time = sum(t["duration"] for t in self.timings.values() if t["duration"] is not None)
            if total_time > 0:
                for counter, value in self.counters.items():
                    if "variants" in counter.lower():
                        rate = value / total_time
                        lines.append(f"  {counter} rate: {rate:.1f}/sec")
        
        return "\n".join(lines)


class BatchOptimizer:
    """Optimize batch processing parameters based on system resources."""
    
    @staticmethod
    def get_optimal_batch_size(
        file_size_mb: float,
        available_memory_mb: float,
        num_samples: int,
        target_memory_usage: float = 0.5
    ) -> int:
        """Calculate optimal batch size based on system resources.
        
        Args:
            file_size_mb: Size of VCF file in MB
            available_memory_mb: Available system memory in MB
            num_samples: Number of samples in VCF
            target_memory_usage: Target fraction of available memory to use
            
        Returns:
            Optimal batch size
        """
        # Estimate memory per variant (rough heuristic)
        # Assumes each variant uses approximately: samples * 10 bytes + overhead
        estimated_bytes_per_variant = num_samples * 10 + 100
        
        # Calculate target memory for processing
        target_memory_bytes = available_memory_mb * target_memory_usage * 1024 * 1024
        
        # Calculate batch size that fits in target memory
        optimal_batch_size = int(target_memory_bytes / estimated_bytes_per_variant)
        
        # Apply reasonable bounds
        min_batch_size = 100
        max_batch_size = 50000
        
        optimal_batch_size = max(min_batch_size, min(optimal_batch_size, max_batch_size))
        
        logger.debug(f"Calculated optimal batch size: {optimal_batch_size} "
                    f"(file: {file_size_mb:.1f}MB, memory: {available_memory_mb:.1f}MB, "
                    f"samples: {num_samples})")
        
        return optimal_batch_size
    
    @staticmethod
    def get_optimal_thread_count(cpu_count: Optional[int] = None) -> int:
        """Get optimal thread count for processing.
        
        Args:
            cpu_count: Override CPU count (uses system CPU count if None)
            
        Returns:
            Optimal thread count
        """
        if cpu_count is None:
            cpu_count = os.cpu_count() or 1
        
        # For I/O bound operations like VCF processing, use fewer threads
        # to avoid overwhelming the system
        optimal_threads = min(cpu_count, 4)  # Cap at 4 threads
        
        # Ensure at least 1 thread
        optimal_threads = max(1, optimal_threads)
        
        logger.debug(f"Calculated optimal thread count: {optimal_threads} (CPU count: {cpu_count})")
        
        return optimal_threads
    
    @staticmethod
    def get_system_info() -> Dict[str, float]:
        """Get current system resource information.
        
        Returns:
            Dictionary with system resource info
        """
        try:
            # Memory information
            memory = psutil.virtual_memory()
            available_memory_mb = memory.available / (1024 * 1024)
            total_memory_mb = memory.total / (1024 * 1024)
            
            # CPU information
            cpu_count = os.cpu_count() or 1
            cpu_percent = psutil.cpu_percent(interval=1)
            
            # Disk information (for output directory)
            disk_usage = psutil.disk_usage('/')
            available_disk_gb = disk_usage.free / (1024 * 1024 * 1024)
            
            return {
                "available_memory_mb": available_memory_mb,
                "total_memory_mb": total_memory_mb,
                "memory_usage_percent": memory.percent,
                "cpu_count": cpu_count,
                "cpu_usage_percent": cpu_percent,
                "available_disk_gb": available_disk_gb
            }
        
        except Exception as e:
            logger.warning(f"Could not get system info: {e}")
            return {
                "available_memory_mb": 1024,  # Default fallback
                "total_memory_mb": 2048,
                "memory_usage_percent": 50,
                "cpu_count": 1,
                "cpu_usage_percent": 50,
                "available_disk_gb": 10
            }


class ProcessingOptimizer:
    """Main optimizer for VCF processing performance."""
    
    def __init__(self, enable_optimization: bool = True):
        """Initialize processing optimizer.
        
        Args:
            enable_optimization: Whether to enable optimization
        """
        self.enable_optimization = enable_optimization
        self.profiler = PerformanceProfiler(enable_optimization)
        self.memory_monitor = self.profiler.memory_monitor
    
    def optimize_config(self, config, vcf_file_path: Path) -> Dict[str, any]:
        """Optimize processing configuration based on system resources and file size.
        
        Args:
            config: Current processing configuration
            vcf_file_path: Path to VCF file
            
        Returns:
            Dictionary with optimization recommendations
        """
        if not self.enable_optimization:
            return {}
        
        recommendations = {}
        
        try:
            # Get system information
            system_info = BatchOptimizer.get_system_info()
            
            # Get file size
            file_size_mb = vcf_file_path.stat().st_size / (1024 * 1024)
            
            # Estimate number of samples (rough heuristic)
            # This would ideally be done by reading the VCF header
            estimated_samples = getattr(config, 'estimated_samples', 10)
            
            # Optimize batch size
            if hasattr(config, 'batch_size'):
                optimal_batch_size = BatchOptimizer.get_optimal_batch_size(
                    file_size_mb=file_size_mb,
                    available_memory_mb=system_info["available_memory_mb"],
                    num_samples=estimated_samples
                )
                
                if optimal_batch_size != config.batch_size:
                    recommendations["batch_size"] = optimal_batch_size
                    recommendations["batch_size_reason"] = (
                        f"Optimized for {system_info['available_memory_mb']:.0f}MB available memory "
                        f"and {file_size_mb:.1f}MB file size"
                    )
            
            # Optimize thread count
            if hasattr(config, 'threads'):
                optimal_threads = BatchOptimizer.get_optimal_thread_count()
                
                if optimal_threads != config.threads:
                    recommendations["threads"] = optimal_threads
                    recommendations["threads_reason"] = (
                        f"Optimized for {system_info['cpu_count']} CPU cores"
                    )
            
            # Memory usage warnings
            if system_info["memory_usage_percent"] > 80:
                recommendations["memory_warning"] = (
                    f"High memory usage ({system_info['memory_usage_percent']:.1f}%). "
                    "Consider reducing batch size or closing other applications."
                )
            
            # Disk space warnings
            if system_info["available_disk_gb"] < 1:
                recommendations["disk_warning"] = (
                    f"Low disk space ({system_info['available_disk_gb']:.1f}GB available). "
                    "Ensure sufficient space for output files."
                )
            
            # Compression recommendations
            if file_size_mb > 100 and not getattr(config, 'compress_output', False):
                recommendations["compression_suggestion"] = (
                    "Consider enabling output compression for large files to save disk space."
                )
            
        except Exception as e:
            logger.warning(f"Could not optimize configuration: {e}")
        
        return recommendations
    
    def start_profiling(self) -> None:
        """Start performance profiling."""
        if self.enable_optimization:
            self.memory_monitor.start()
            self.profiler.start_timer("total_processing")
    
    def profile_operation(self, operation: str):
        """Context manager for profiling an operation.
        
        Args:
            operation: Name of the operation to profile
        """
        return ProfiledOperation(self.profiler, operation)
    
    def end_profiling(self) -> str:
        """End profiling and return performance report.
        
        Returns:
            Performance report string
        """
        if self.enable_optimization:
            self.profiler.end_timer("total_processing")
            return self.profiler.get_performance_report()
        return "Performance profiling disabled"


class ProfiledOperation:
    """Context manager for profiling individual operations."""
    
    def __init__(self, profiler: PerformanceProfiler, operation: str):
        """Initialize profiled operation.
        
        Args:
            profiler: Performance profiler instance
            operation: Name of the operation
        """
        self.profiler = profiler
        self.operation = operation
    
    def __enter__(self):
        """Start profiling the operation."""
        self.profiler.start_timer(self.operation)
        self.profiler.memory_monitor.sample(f"start_{self.operation}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """End profiling the operation."""
        self.profiler.end_timer(self.operation)
        self.profiler.memory_monitor.sample(f"end_{self.operation}")
        
        # Force garbage collection for memory-intensive operations
        if "process" in self.operation.lower() or "convert" in self.operation.lower():
            self.profiler.memory_monitor.force_gc()


def optimize_for_large_files(file_size_mb: float) -> Dict[str, any]:
    """Get optimization recommendations for large files.
    
    Args:
        file_size_mb: File size in MB
        
    Returns:
        Dictionary with optimization recommendations
    """
    recommendations = {}
    
    if file_size_mb > 1000:  # Files larger than 1GB
        recommendations.update({
            "batch_size": 50000,
            "compress_output": True,
            "memory_monitoring": True,
            "gc_frequency": 1000,  # Force GC every 1000 variants
            "progress_reporting": True
        })
        recommendations["optimization_note"] = "Large file detected: using memory-efficient settings"
    
    elif file_size_mb > 100:  # Files larger than 100MB
        recommendations.update({
            "batch_size": 20000,
            "compress_output": True,
            "progress_reporting": True
        })
        recommendations["optimization_note"] = "Medium file detected: using balanced settings"
    
    else:  # Small files
        recommendations.update({
            "batch_size": 5000,
            "compress_output": False,
            "progress_reporting": False
        })
        recommendations["optimization_note"] = "Small file detected: using fast processing settings"
    
    return recommendations