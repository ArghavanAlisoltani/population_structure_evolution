#!/bin/bash

# Define the input file containing the changes
INPUT_FILE="changed_region.txt"

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

# Read the input file line by line
# IFS=$'\t' sets the delimiter to a tab
while IFS=$'\t' read -r filename string_to_remove string_to_replace; do
    
    # Check if the target .sh file exists before trying to modify it
    if [ -f "$filename" ]; then
        echo "Processing '$filename': Replacing '$string_to_remove' with '$string_to_replace'"
        
        # Use sed to perform the in-place replacement.
        # Double quotes are used to allow variable expansion.
        sed -i "s/$string_to_remove/$string_to_replace/g" "$filename"
    else
        echo "Warning: File '$filename' not found. Skipping."
    fi

done < "$INPUT_FILE"

echo "Script finished."
