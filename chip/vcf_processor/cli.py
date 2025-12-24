"""Command-line interface for VCF processor."""

from .config_manager import create_cli_app


def main() -> None:
    """Main CLI entry point."""
    app = create_cli_app()
    app()


if __name__ == "__main__":
    main()