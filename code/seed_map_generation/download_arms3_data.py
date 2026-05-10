#!/usr/bin/env python3
"""
Script to sync S3 data for participants with matched_group == 3
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description='Sync S3 data for participants in matched_group 3'
    )
    parser.add_argument(
        'tsv_file',
        help='Path to TSV file with participant data'
    )
    parser.add_argument(
        '--ses',
        required=True,
        help='Session identifier (e.g., ses-01)'
    )
    parser.add_argument(
        '--bucket',
        default='bucket',
        help='S3 bucket name (default: bucket)'
    )
    parser.add_argument(
        '--s3-path',
        default='path',
        help='Path within S3 bucket (default: path)'
    )
    parser.add_argument(
        '--output-dir',
        default='output/dir',
        help='Local output directory (default: output/dir)'
    )
    parser.add_argument(
        '--s3-prefix',
        default='task-rest_bold_desc-filtered_timeseries.dtseries.nii',
        help='file prefix to download from s3 (default: task-rest_bold_desc-filtered_timeseries.dtseries.nii)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print commands without executing them'
    )
    return parser.parse_args()


def read_matched_participants(tsv_file):
    """Read TSV and return participant_ids where matched_group == 3"""
    participants = []
    
    try:
        with open(tsv_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
                if 'matched_group' not in row or 'participant_id' not in row:
                    print("Error: TSV must contain 'matched_group' and 'participant_id' columns", 
                          file=sys.stderr)
                    sys.exit(1)
                
                if row['matched_group'].strip() == '3':
                    participants.append(row['participant_id'].strip())
        
        return participants
    
    except FileNotFoundError:
        print(f"Error: File '{tsv_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading TSV file: {e}", file=sys.stderr)
        sys.exit(1)


def sync_participant_data(participant_id, ses, bucket, s3_path, output_dir, s3_prefix, dry_run):
    """Execute s3cdm sync for a single participant"""
    s3_source = f"s3://{bucket}/{s3_path}/{participant_id}/{ses}/func/{participant_id}_{ses}_{s3_prefix}"
    local_dest = f"{output_dir}/{participant_id}/{ses}/func/"
    
    # Create local directory if it doesn't exist
    if not dry_run:
        Path(local_dest).mkdir(parents=True, exist_ok=True)
    
    cmd = ['s3cmd', 'sync', s3_source, local_dest, '-v', '-f']
    
    print(f"Syncing {participant_id}...")
    print(f"  Command: {' '.join(cmd)}")
    
    if dry_run:
        print("  [DRY RUN - not executing]")
        # Return the full path to the downloaded file
        local_file = Path(local_dest) / f"{participant_id}_{ses}_{s3_prefix}"
        return True, str(local_file.absolute())
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            print(f"  {result.stdout}")
            # Return the full path to the downloaded file
            local_file = Path(local_dest) / f"{participant_id}_{ses}_{s3_prefix}"
            print(f"  Downloaded file path: {local_file.absolute()}")
            return True, str(local_file.absolute())
        else:
            print(f"  No output from sync command for {participant_id}")
            return False, None
    except subprocess.CalledProcessError as e:
        if "WARNING: Can't set owner/group" in e.stderr:
            print(f"  Warning syncing {participant_id}: {e}", file=sys.stderr)
            if e.stderr:
                print(f"  {e.stderr}", file=sys.stderr)
            # Return the full path to the downloaded file even if there was a warning
            local_file = Path(local_dest) / f"{participant_id}_{ses}_{s3_prefix}"
            print(f"  Downloaded file path (with warning): {local_file.absolute()}")
            return True, str(local_file.absolute())
        
        print(f"  Error syncing {participant_id}: {e}", file=sys.stderr)
        if e.stderr:
            print(f"  {e.stderr}", file=sys.stderr)
        return False, None


def main():
    args = parse_args()
    
    print(f"Reading participants from: {args.tsv_file}")
    participants = read_matched_participants(args.tsv_file)
    
    if not participants:
        print("No participants found with matched_group == 3")
        return
    
    print(f"\nFound {len(participants)} participant(s) with matched_group == 3:")
    for p in participants:
        print(f"  - {p}")
    
    print(f"\nSession: {args.ses}")
    print(f"S3 bucket: {args.bucket}")
    print(f"S3 path: {args.s3_path}")
    print(f"S3 prefix: {args.s3_prefix}")
    print(f"Output directory: {args.output_dir}")
    
    if args.dry_run:
        print("\n*** DRY RUN MODE ***\n")
    
    success_count = 0
    fail_count = 0
    downloaded_files = []  # List to store full paths of downloaded files
    
    print("\nStarting sync operations...\n")
    for participant_id in participants:
        success, file_path = sync_participant_data(participant_id, args.ses, args.bucket, 
                                args.s3_path, args.output_dir, args.s3_prefix, args.dry_run)
        if success:
            success_count += 1
            if file_path:  # Add the full path to downloaded_files list
                downloaded_files.append(file_path)
        else:
            print(f"Failed to sync data for participant: {participant_id}")
            fail_count += 1
        print()
    
    print(f"Summary: {success_count} successful, {fail_count} failed")
    
    # Write .conc file with full paths to downloaded files
    if downloaded_files and not args.dry_run:
        conc_file = Path(args.output_dir) / f"downloaded_dtseries_files_{args.ses}.conc"
        try:
            with open(conc_file, 'w') as f:
                for file_path in downloaded_files:
                    f.write(f"{file_path}\n")
            print(f"\nCreated .conc file: {conc_file.absolute()}")
            print(f"Contains {len(downloaded_files)} file path(s)")
        except Exception as e:
            print(f"Error writing .conc file: {e}", file=sys.stderr)
    
    print(f"Summary: {success_count} successful, {fail_count} failed")


if __name__ == '__main__':
    main()