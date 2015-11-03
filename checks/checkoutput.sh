unset GIT_DIR
cmd="$(git rev-parse --show-toplevel)/checks/checkpatch.pl"
git diff origin/master --staged | eval "\"$cmd\"" --no-tree - > style_check_output.txt
