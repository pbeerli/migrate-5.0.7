#!/bin/bash 
# used to update Makefile with the last git version available

# 1. Initialize the default version (this is the fallback)
MYVERSION="Distribution-version"

# 2. Check if the 'git' command exists
if [ -x "$(command -v git)" ]; then
    
    # Use a temporary variable to capture the result of the git command
    # We run the command, and use '|| true' to ensure the script continues 
    # even if git fails, while capturing the *real* failure state.
    git_version=$(git describe --abbrev=7 --dirty --always --tag --long 2> /dev/null)
    
    # 3. CRITICAL CHECK: Check the exit status of the last command (git describe)
    # $? is the exit code of the previous command. 0 means success.
    if [ $? -eq 0 ]; then
        # Command succeeded (it found a tag/version)
        MYVERSION="$git_version"
    else
        # Command failed (it probably meant to fall back)
        # In this case, MYVERSION remains the default set at the start.
        : # Do nothing, let it keep the default value
    fi
fi

echo $MYVERSION
