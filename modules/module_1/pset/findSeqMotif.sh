# BASH script to search for a sequence motif in a FASTA file
# Arguments: 
### 1. FASTA file
### 2. target motif as a regular expression 
# Returns the number of times the motif appears

cat $1 | tail -n +2 |
tr -d "\n" |
grep -o $2 | wc -l