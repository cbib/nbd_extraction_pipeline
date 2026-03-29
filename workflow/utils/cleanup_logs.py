import os
import shutil
from datetime import datetime, timedelta
from pathlib import Path

# Define paths
log_dir = Path("../../.snakemake/log")
archive_dir = log_dir / "log_archive"

# Check if log directory exists
if not log_dir.exists():
    print(f"Log directory does not exist: {log_dir}")
else:
    # Create archive directory if it doesn't exist
    archive_dir.mkdir(parents=True, exist_ok=True)

    # Get the date threshold (one week ago)
    week_ago = datetime.now() - timedelta(days=7)

    # Iterate through log files
    for log_file in log_dir.iterdir():
        # Skip the archive directory itself
        if log_file == archive_dir:
            continue

        # Check if it's a file (not a directory)
        if log_file.is_file():
            # Get file modification time
            file_mtime = datetime.fromtimestamp(log_file.stat().st_mtime)

            # Move file if older than a week
            if file_mtime < week_ago:
                destination = archive_dir / log_file.name
                shutil.move(str(log_file), str(destination))
                print(f"Archived: {log_file.name}")

    print("Log archiving complete.")
