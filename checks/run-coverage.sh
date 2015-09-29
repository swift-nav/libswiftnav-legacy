# requires gcovr and diff-cover
# pip install gcovr diff-cover
set -e
base="$(git rev-parse --show-toplevel)"
echo "running coverage check from ${base}"
gcovr -r "${base}" --object-directory . --xml -o "${base}/.coverage.xml"
gcovr -r "${base}" --object-directory . --html -o coverage.html
( cd "${base}" && diff-cover .coverage.xml --fail-under=85 --html-report .diffreport.html )
mv "${base}/.diffreport.html" ./diffreport.html
