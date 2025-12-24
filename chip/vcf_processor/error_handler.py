"""Error handling and recovery utilities for VCF processor."""

import traceback
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Type, Union

from loguru import logger

from .config import ProcessingResult


class ErrorSeverity(Enum):
    """Error severity levels."""
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


class ProcessingError(Exception):
    """Base exception for VCF processing errors."""
    
    def __init__(
        self,
        message: str,
        severity: ErrorSeverity = ErrorSeverity.ERROR,
        context: Optional[Dict[str, Any]] = None,
        recoverable: bool = False
    ) -> None:
        """Initialize processing error.
        
        Args:
            message: Error message
            severity: Error severity level
            context: Additional context information
            recoverable: Whether the error is recoverable
        """
        super().__init__(message)
        self.message = message
        self.severity = severity
        self.context = context or {}
        self.recoverable = recoverable


class ValidationError(ProcessingError):
    """Error in data validation."""
    pass


class FileProcessingError(ProcessingError):
    """Error in file processing."""
    pass


class MemoryError(ProcessingError):
    """Error related to memory usage."""
    pass


class ConfigurationError(ProcessingError):
    """Error in configuration."""
    pass


class ErrorHandler:
    """Handles errors and implements recovery strategies.
    
    This class provides comprehensive error handling with different
    recovery strategies based on error type and severity.
    
    Attributes:
        result: Processing result to track errors
        recovery_strategies: Dictionary of recovery strategies
        max_retries: Maximum number of retry attempts
        _error_counts: Count of errors by type
    """
    
    def __init__(
        self,
        result: ProcessingResult,
        max_retries: int = 3
    ) -> None:
        """Initialize error handler.
        
        Args:
            result: Processing result to track errors
            max_retries: Maximum number of retry attempts
        """
        self.result = result
        self.max_retries = max_retries
        self._error_counts: Dict[str, int] = {}
        self.recovery_strategies: Dict[Type[Exception], Callable] = {
            FileNotFoundError: self._handle_file_not_found,
            PermissionError: self._handle_permission_error,
            MemoryError: self._handle_memory_error,
            ValidationError: self._handle_validation_error,
            FileProcessingError: self._handle_file_processing_error,
        }
        
        logger.debug("ErrorHandler initialized")
    
    def handle_error(
        self,
        error: Exception,
        context: Optional[Dict[str, Any]] = None,
        operation: str = "unknown"
    ) -> bool:
        """Handle an error with appropriate recovery strategy.
        
        Args:
            error: The exception that occurred
            context: Additional context information
            operation: Description of the operation that failed
            
        Returns:
            True if error was handled and operation can continue, False otherwise
        """
        error_type = type(error).__name__
        self._error_counts[error_type] = self._error_counts.get(error_type, 0) + 1
        
        # Log the error
        logger.error(
            f"Error in {operation}: {error_type}: {str(error)}"
        )
        
        # Add to result
        error_msg = f"{operation}: {error_type}: {str(error)}"
        self.result.add_error(error_msg)
        
        # Try recovery strategy
        recovery_func = self.recovery_strategies.get(type(error))
        if recovery_func:
            try:
                return recovery_func(error, context, operation)
            except Exception as recovery_error:
                logger.error(f"Recovery strategy failed: {recovery_error}")
                return False
        
        # No recovery strategy available
        logger.warning(f"No recovery strategy for {error_type}")
        return False
    
    def _handle_file_not_found(self, error: FileNotFoundError, context: Dict, operation: str) -> bool:
        """Handle file not found errors.
        
        Args:
            error: The FileNotFoundError
            context: Error context
            operation: Operation description
            
        Returns:
            True if recovered, False otherwise
        """
        logger.info(f"Attempting to recover from file not found: {error}")
        
        # Try to create missing directories
        if context and 'file_path' in context:
            file_path = Path(context['file_path'])
            
            # Create parent directories if they don't exist
            try:
                file_path.parent.mkdir(parents=True, exist_ok=True)
                logger.info(f"Created missing directory: {file_path.parent}")
                
                # If it's an output file, create an empty file
                if 'output' in operation.lower():
                    file_path.touch()
                    logger.info(f"Created empty output file: {file_path}")
                    return True
                    
            except Exception as e:
                logger.error(f"Failed to create directory/file: {e}")
        
        return False
    
    def _handle_permission_error(self, error: PermissionError, context: Dict, operation: str) -> bool:
        """Handle permission errors.
        
        Args:
            error: The PermissionError
            context: Error context
            operation: Operation description
            
        Returns:
            True if recovered, False otherwise
        """
        logger.warning(f"Permission error in {operation}: {error}")
        
        # Try alternative output location
        if context and 'file_path' in context:
            file_path = Path(context['file_path'])
            
            # Try writing to temp directory
            import tempfile
            temp_dir = Path(tempfile.gettempdir())
            alt_path = temp_dir / file_path.name
            
            try:
                # Test if we can write to temp directory
                test_file = alt_path.with_suffix('.test')
                test_file.touch()
                test_file.unlink()
                
                logger.info(f"Alternative path available: {alt_path}")
                context['alternative_path'] = alt_path
                return True
                
            except Exception as e:
                logger.error(f"Alternative path also failed: {e}")
        
        return False
    
    def _handle_memory_error(self, error: MemoryError, context: Dict, operation: str) -> bool:
        """Handle memory errors.
        
        Args:
            error: The MemoryError
            context: Error context
            operation: Operation description
            
        Returns:
            True if recovered, False otherwise
        """
        logger.warning(f"Memory error in {operation}: {error}")
        
        # Suggest reducing batch size
        if context and 'batch_size' in context:
            current_batch_size = context['batch_size']
            new_batch_size = max(100, current_batch_size // 2)
            
            logger.info(f"Reducing batch size from {current_batch_size} to {new_batch_size}")
            context['suggested_batch_size'] = new_batch_size
            return True
        
        # Suggest garbage collection
        try:
            import gc
            gc.collect()
            logger.info("Performed garbage collection")
            return True
        except Exception as e:
            logger.error(f"Garbage collection failed: {e}")
        
        return False
    
    def _handle_validation_error(self, error: ValidationError, context: Dict, operation: str) -> bool:
        """Handle validation errors.
        
        Args:
            error: The ValidationError
            context: Error context
            operation: Operation description
            
        Returns:
            True if recovered, False otherwise
        """
        logger.warning(f"Validation error in {operation}: {error}")
        
        # For validation errors, we might skip the problematic data
        if hasattr(error, 'recoverable') and error.recoverable:
            logger.info("Validation error is recoverable, skipping problematic data")
            return True
        
        return False
    
    def _handle_file_processing_error(self, error: FileProcessingError, context: Dict, operation: str) -> bool:
        """Handle file processing errors.
        
        Args:
            error: The FileProcessingError
            context: Error context
            operation: Operation description
            
        Returns:
            True if recovered, False otherwise
        """
        logger.warning(f"File processing error in {operation}: {error}")
        
        # Try to skip the problematic file/record
        if hasattr(error, 'recoverable') and error.recoverable:
            logger.info("File processing error is recoverable, skipping problematic data")
            return True
        
        return False
    
    def retry_operation(
        self,
        operation: Callable,
        *args,
        max_retries: Optional[int] = None,
        **kwargs
    ) -> Any:
        """Retry an operation with exponential backoff.
        
        Args:
            operation: Function to retry
            *args: Arguments for the operation
            max_retries: Override default max retries
            **kwargs: Keyword arguments for the operation
            
        Returns:
            Result of the operation
            
        Raises:
            Exception: If all retries are exhausted
        """
        import time
        
        max_retries = max_retries or self.max_retries
        last_error = None
        
        for attempt in range(max_retries + 1):
            try:
                return operation(*args, **kwargs)
            except Exception as e:
                last_error = e
                
                if attempt < max_retries:
                    wait_time = 2 ** attempt  # Exponential backoff
                    logger.warning(
                        f"Operation failed (attempt {attempt + 1}/{max_retries + 1}): {e}. "
                        f"Retrying in {wait_time} seconds..."
                    )
                    time.sleep(wait_time)
                else:
                    logger.error(f"Operation failed after {max_retries + 1} attempts")
        
        # All retries exhausted
        raise last_error
    
    def get_error_summary(self) -> Dict[str, Any]:
        """Get summary of errors encountered.
        
        Returns:
            Dictionary with error statistics
        """
        total_errors = sum(self._error_counts.values())
        
        return {
            "total_errors": total_errors,
            "error_counts": self._error_counts.copy(),
            "most_common_error": max(self._error_counts.items(), key=lambda x: x[1])[0] if self._error_counts else None,
            "error_types": list(self._error_counts.keys())
        }
    
    def reset_error_counts(self) -> None:
        """Reset error counters."""
        self._error_counts.clear()
        logger.debug("Error counts reset")


class RecoveryManager:
    """Manages recovery operations and checkpoints."""
    
    def __init__(self, checkpoint_dir: Optional[Path] = None) -> None:
        """Initialize recovery manager.
        
        Args:
            checkpoint_dir: Directory for storing checkpoints
        """
        self.checkpoint_dir = checkpoint_dir or Path.cwd() / ".checkpoints"
        self.checkpoint_dir.mkdir(exist_ok=True)
        
        logger.debug(f"RecoveryManager initialized with checkpoint dir: {self.checkpoint_dir}")
    
    def create_checkpoint(self, name: str, data: Dict[str, Any]) -> Path:
        """Create a recovery checkpoint.
        
        Args:
            name: Checkpoint name
            data: Data to save in checkpoint
            
        Returns:
            Path to checkpoint file
        """
        import json
        import time
        
        timestamp = int(time.time())
        checkpoint_file = self.checkpoint_dir / f"{name}_{timestamp}.json"
        
        try:
            with open(checkpoint_file, 'w') as f:
                json.dump(data, f, indent=2, default=str)
            
            logger.info(f"Checkpoint created: {checkpoint_file}")
            return checkpoint_file
            
        except Exception as e:
            logger.error(f"Failed to create checkpoint: {e}")
            raise
    
    def load_checkpoint(self, name: str) -> Optional[Dict[str, Any]]:
        """Load the most recent checkpoint.
        
        Args:
            name: Checkpoint name pattern
            
        Returns:
            Checkpoint data or None if not found
        """
        import json
        
        # Find most recent checkpoint
        pattern = f"{name}_*.json"
        checkpoints = list(self.checkpoint_dir.glob(pattern))
        
        if not checkpoints:
            logger.info(f"No checkpoints found for pattern: {pattern}")
            return None
        
        # Sort by timestamp (newest first)
        latest_checkpoint = max(checkpoints, key=lambda p: p.stat().st_mtime)
        
        try:
            with open(latest_checkpoint, 'r') as f:
                data = json.load(f)
            
            logger.info(f"Loaded checkpoint: {latest_checkpoint}")
            return data
            
        except Exception as e:
            logger.error(f"Failed to load checkpoint {latest_checkpoint}: {e}")
            return None
    
    def cleanup_old_checkpoints(self, name: str, keep_count: int = 5) -> None:
        """Clean up old checkpoints, keeping only the most recent ones.
        
        Args:
            name: Checkpoint name pattern
            keep_count: Number of checkpoints to keep
        """
        pattern = f"{name}_*.json"
        checkpoints = list(self.checkpoint_dir.glob(pattern))
        
        if len(checkpoints) <= keep_count:
            return
        
        # Sort by timestamp (oldest first)
        checkpoints.sort(key=lambda p: p.stat().st_mtime)
        
        # Remove old checkpoints
        for checkpoint in checkpoints[:-keep_count]:
            try:
                checkpoint.unlink()
                logger.debug(f"Removed old checkpoint: {checkpoint}")
            except Exception as e:
                logger.warning(f"Failed to remove checkpoint {checkpoint}: {e}")


def with_error_handling(
    error_handler: ErrorHandler,
    operation_name: str = "operation",
    context: Optional[Dict[str, Any]] = None
):
    """Decorator for adding error handling to functions.
    
    Args:
        error_handler: ErrorHandler instance
        operation_name: Name of the operation for logging
        context: Additional context for error handling
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                handled = error_handler.handle_error(e, context, operation_name)
                if not handled:
                    raise
                return None
        return wrapper
    return decorator