# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
# Give bash history a decent size
export HISTFILESIZE=10000
# Helpful alias for users to quickly see their queue status
alias squ='squeue -u $USER'

module load R/4.0.5
